module Poseidon.Generator.CLI.AdmixPops where

import           Poseidon.Generator.Parsers
import           Poseidon.Generator.SampleGeno
import           Poseidon.Generator.Types
import           Poseidon.Generator.Utils

import           Control.Exception             (catch, throwIO)
import           Control.Monad                 (forM, unless, when)
import           Control.Monad.Reader          (ask)
import           Data.Function                 ((&))
import           Data.List
import           Data.Maybe
import           Data.Time                     (getCurrentTime)
import           Lens.Family2                  (view)
import           Pipes
import qualified Pipes.Group                   as PG
import qualified Pipes.Prelude                 as P
import           Pipes.Safe                    (runSafeT)
import           Poseidon.GenotypeData
import           Poseidon.Janno
import           Poseidon.Package
import           Poseidon.Utils
import           SequenceFormats.Eigenstrat
import           SequenceFormats.Plink         (writePlink)
import           System.Directory              (createDirectoryIfMissing)
import           System.FilePath               (takeBaseName, (<.>), (</>))

data AdmixPopsMethodSettings =
    PerSNP {
        _admixMarginalizeMissing :: Bool
    } | InChunks {
        _admixChunkSize :: Int
    }

data AdmixPopsOptions = AdmixPopsOptions {
      _admixGenoSources             :: [GenoDataSource]
    , _admixIndWithAdmixtureSet     :: [InIndAdmixpops]
    , _admixIndWithAdmixtureSetFile :: Maybe FilePath
    , _admixMethodSettings          :: AdmixPopsMethodSettings
    , _admixOutFormat               :: GenotypeFormatSpec
    , _admixOutPath                 :: FilePath
    , _forgeOutPacName              :: Maybe String
    }

pacReadOpts :: PackageReadOptions
pacReadOpts = defaultPackageReadOptions {
      _readOptVerbose          = False
    , _readOptStopOnDuplicates = False
    , _readOptIgnoreChecksums  = True
    , _readOptIgnoreGeno       = False
    , _readOptGenoCheck        = True
    }

runAdmixPops :: AdmixPopsOptions -> PoseidonLogIO ()
runAdmixPops (
    AdmixPopsOptions
        genoSources
        popsWithFracsDirect
        popsWithFracsFile
        methodSetting
        outFormat
        outPath
        maybeOutName
    ) = do
    -- compile individuals
    popsWithFracsFromFile <- case popsWithFracsFile of
        Nothing -> return []
        Just f  -> liftIO $ readIndWithAdmixtureSetFromFile f
    let requestedInds = popsWithFracsDirect ++ popsWithFracsFromFile
    -- validating input
    logInfo "Checking requested, artificial individuals"
    logInfo $ "Individuals: " ++ renderRequestedInds requestedInds
    liftIO $ checkIndsWithAdmixtureSets requestedInds
    -- load Poseidon packages
    properPackages <- readPoseidonPackageCollection pacReadOpts $ [getPacBaseDirs x | x@PacBaseDir {} <- genoSources]
    pseudoPackages <- liftIO $ mapM makePseudoPackageFromGenotypeData $ [getGenoDirect x | x@GenoDirect {} <- genoSources]
    logInfo $ "Unpackaged genotype data files loaded: " ++ show (length pseudoPackages)
    let allPackages = properPackages ++ pseudoPackages
    -- determine relevant packages
    let popNames = concat $ map (map _inPopName) $ map _inPopSet requestedInds
    relevantPackages <- liftIO $ filterPackagesByPops popNames allPackages
    -- gather additional info for requested inds
    preparedInds <- liftIO $ mapM (`gatherInfoForInd` relevantPackages) requestedInds
    -- create new package --
    let outName = case maybeOutName of -- take basename of outPath, if name is not provided
            Just x  -> x
            Nothing -> takeBaseName outPath
    when (outName == "") $ liftIO $ throwIO PoseidonEmptyOutPacNameException
    -- create new directory
    logInfo $ "Writing to directory (will be created if missing): " ++ outPath
    liftIO $ createDirectoryIfMissing True outPath
    -- compile genotype data structure
    let [outInd, outSnp, outGeno] = case outFormat of
            GenotypeFormatEigenstrat -> [outName <.> ".ind", outName <.> ".snp", outName <.> ".geno"]
            GenotypeFormatPlink -> [outName <.> ".fam", outName <.> ".bim", outName <.> ".bed"]
    let [outG, outS, outI] = map (outPath </>) [outGeno, outSnp, outInd]
    let genotypeData = GenotypeDataSpec outFormat outGeno Nothing outSnp Nothing outInd Nothing Nothing
        pac = newMinimalPackageTemplate outPath "admixpops_package" genotypeData
    liftIO $ writePoseidonPackage pac
    -- compile genotype data
    logInfo "Compiling individuals"
    logEnv <- ask
    currentTime <- liftIO getCurrentTime
    liftIO $ catch (
        runSafeT $ do
            (_, eigenstratProd) <- getJointGenotypeData logEnv False relevantPackages Nothing
            let newIndEntries = map (\x -> EigenstratIndEntry (_indName x) Unknown (_groupName x)) preparedInds
            let outConsumer = case outFormat of
                    GenotypeFormatEigenstrat -> writeEigenstrat outG outS outI newIndEntries
                    GenotypeFormatPlink      -> writePlink      outG outS outI newIndEntries
            case methodSetting of
                PerSNP marginalizeMissing -> do
                    runEffect $ eigenstratProd >->
                        printSNPCopyProgress logEnv currentTime >->
                        P.mapM (samplePerSNP marginalizeMissing preparedInds) >->
                        outConsumer
                InChunks chunkSize -> do
                    runEffect $ (
                            eigenstratProd &
                            chunkEigenstratByNrSnps chunkSize &
                            PG.maps (samplePerChunk logEnv preparedInds) &
                            PG.concats
                        ) >->
                        printSNPCopyProgress logEnv currentTime >->
                        outConsumer
        ) (\e -> throwIO $ PoseidonGenotypeExceptionForward e)
    logInfo "Done"
    where
        chunkEigenstratByNrSnps chunkSize = view (PG.chunksOf chunkSize)

checkIndsWithAdmixtureSets :: [InIndAdmixpops] -> IO ()
checkIndsWithAdmixtureSets requestedInds = do
    checkDuplicateIndNames requestedInds
    mapM_ checkPopFracList requestedInds
    where
        checkDuplicateIndNames :: [InIndAdmixpops] -> IO ()
        checkDuplicateIndNames xs =
            let individualsGrouped = filter (\x -> length x > 1) $ group $ sort $ map _inIndName xs
            in unless (null individualsGrouped) $ do
                throwIO $ PoseidonGeneratorCLIParsingException $
                    "Duplicate individual names: " ++ intercalate "," (nub $ concat individualsGrouped)
        checkPopFracList :: InIndAdmixpops -> IO ()
        checkPopFracList x = do
            let xs = _inPopSet x
                fracs = map _inPopFrac xs
            when (sum fracs /= 1) $ do
                throwIO $ PoseidonGeneratorCLIParsingException $
                    "Fractions in " ++ show x ++ " do not to sum to 100%"

filterPackagesByPops :: [String] -> [PoseidonPackage] -> IO [PoseidonPackage]
filterPackagesByPops pops packages = do
    fmap catMaybes . forM packages $ \pac -> do
        inds <- loadIndividuals (posPacBaseDir pac) (posPacGenotypeData pac)
        let groupNamesPac = [groupName | EigenstratIndEntry _ _ groupName <- inds]
        if   not (null (groupNamesPac `intersect` pops))
        then return (Just pac)
        else return Nothing

gatherInfoForInd :: InIndAdmixpops -> [PoseidonPackage] -> IO IndAdmixpops
gatherInfoForInd (InIndAdmixpops name_ group_ set_) pacs = do
    inds <- mapM (`extractIndsPerPop` pacs) set_
    return $ IndAdmixpops name_ group_ inds

extractIndsPerPop :: InPopAdmixpops -> [PoseidonPackage] -> IO PopAdmixpops
extractIndsPerPop (InPopAdmixpops pop_ frac_) relevantPackages = do
    let allPackageNames = map posPacTitle relevantPackages
    allIndEntries <- mapM (\pac -> loadIndividuals (posPacBaseDir pac) (posPacGenotypeData pac)) relevantPackages
    let filterFunc (_,_,EigenstratIndEntry _ _ _group) = _group == pop_
        indNames = map extractIndName $ filter filterFunc (zipGroup allPackageNames allIndEntries)
        indIDs = map extractFirst $ filter filterFunc (zipGroup allPackageNames allIndEntries)
    return (PopAdmixpops pop_ frac_ (zip indNames indIDs))
    where
        extractIndName (_,_,EigenstratIndEntry x _ _) = x
