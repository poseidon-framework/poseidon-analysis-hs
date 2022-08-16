module Poseidon.Generator.CLI.AdmixPops where

import           Poseidon.Generator.Parsers
import           Poseidon.Generator.SampleGeno
import           Poseidon.Generator.Types
import           Poseidon.Generator.Utils

import           Control.Exception             (catch, throwIO)
import           Control.Monad                 (forM, unless, when)
import           Control.Monad.Reader          (ask)
import           Data.List
import           Data.Maybe
import           Data.Ratio                    ((%))
import           Data.Time                     (getCurrentTime)
import           Pipes
import qualified Pipes.Prelude                 as P
import           Pipes.Safe                    (SafeT, runSafeT)
import           Poseidon.GenotypeData
import           Poseidon.Janno
import           Poseidon.Package
import           Poseidon.Utils
import           SequenceFormats.Eigenstrat
import           SequenceFormats.Plink         (writePlink)
import           System.Directory              (createDirectoryIfMissing)
import           System.FilePath               (takeBaseName, (<.>), (</>))
import Lens.Family2 (view)
import qualified Pipes.Group as PG
import Data.Function ((&))

data AdmixPopsOptions = AdmixPopsOptions {
      _admixGenoSources             :: [GenoDataSource]
    , _admixIndWithAdmixtureSet     :: [IndWithAdmixtureSet]
    , _admixIndWithAdmixtureSetFile :: Maybe FilePath
    , _admixMarginalizeMissing      :: Bool
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
runAdmixPops (AdmixPopsOptions genoSources popsWithFracsDirect popsWithFracsFile marginalizeMissing outFormat outPath maybeOutName) = do
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
    -- determine relevant packages and indices
    let popsWithFracs = map (_popFracList . _admixSet) requestedInds
        pops = map (map pop) popsWithFracs
    relevantPackages <- liftIO $ filterPackagesByPops (concat pops) allPackages
    popsFracsInds <- liftIO $  mapM (mapM (`extractIndsPerPop` relevantPackages)) popsWithFracs
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
            let newIndEntries = map (\x -> EigenstratIndEntry (_admixInd x) Unknown (_admixUnit x)) requestedInds
            let outConsumer = case outFormat of
                    GenotypeFormatEigenstrat -> writeEigenstrat outG outS outI newIndEntries
                    GenotypeFormatPlink      -> writePlink      outG outS outI newIndEntries
            if True
            then do 
                runEffect $ eigenstratProd >->
                    printSNPCopyProgress logEnv currentTime >->
                    P.mapM (sampleGenoForMultipleIndWithAdmixtureSet marginalizeMissing popsFracsInds) >->
                    outConsumer
            else do
                runEffect $ (
                        eigenstratProd &
                        chunkEigenstratByNrSnps 5000 &
                        PG.maps sampleChunk &
                        PG.concats
                    ) >->
                    printSNPCopyProgress logEnv currentTime >->
                    outConsumer
        ) (\e -> throwIO $ PoseidonGenotypeExceptionForward e)
    logInfo "Done"
    where
        chunkEigenstratByNrSnps chunkSize = view (PG.chunksOf chunkSize)

sampleChunk :: Producer (EigenstratSnpEntry, GenoLine) (SafeT IO) r ->
               Producer (EigenstratSnpEntry, GenoLine) (SafeT IO) r
sampleChunk prod = for prod handleEntry
  where
    handleEntry :: (EigenstratSnpEntry, GenoLine) -> Producer (EigenstratSnpEntry, GenoLine) (SafeT IO) ()
    handleEntry x = do
        yield x


renderRequestedInds :: [IndWithAdmixtureSet] -> String
renderRequestedInds requestedInds =
    let indString = intercalate ";" $ map show $ take 5 requestedInds
    in if length requestedInds > 5
       then indString ++ "..."
       else indString

checkIndsWithAdmixtureSets :: [IndWithAdmixtureSet] -> IO ()
checkIndsWithAdmixtureSets requestedInds = do
    checkDuplicateIndNames requestedInds
    mapM_ checkPopFracList requestedInds
    where
        checkDuplicateIndNames :: [IndWithAdmixtureSet] -> IO ()
        checkDuplicateIndNames xs =
            let individualsGrouped = filter (\x -> length x > 1) $ group $ sort $ map _admixInd xs
            in unless (null individualsGrouped) $ do
                throwIO $ PoseidonGeneratorCLIParsingException $
                    "Duplicate individual names: " ++ intercalate "," (nub $ concat individualsGrouped)
        checkPopFracList :: IndWithAdmixtureSet -> IO ()
        checkPopFracList x = do
            let xs = (_popFracList . _admixSet) x
                fracs = map frac xs
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

extractIndsPerPop :: PopulationWithFraction -> [PoseidonPackage] -> IO ([Int], Rational)
extractIndsPerPop (PopulationWithFraction _pop _frac) relevantPackages = do
    let allPackageNames = map posPacTitle relevantPackages
    allIndEntries <- mapM (\pac -> loadIndividuals (posPacBaseDir pac) (posPacGenotypeData pac)) relevantPackages
    let filterFunc (_,_,EigenstratIndEntry _ _ _group) = _group == _pop
    return (map extractFirst $ filter filterFunc (zipGroup allPackageNames allIndEntries), _frac)
