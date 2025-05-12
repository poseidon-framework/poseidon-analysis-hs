module Poseidon.Generator.CLI.AdmixPops where

import           Poseidon.Generator.Parsers
import           Poseidon.Generator.SampleGeno
import           Poseidon.Generator.Types
import           Poseidon.Generator.Utils

import           Control.Exception             (catch, throwIO)
import           Control.Monad                 (forM, unless, when)
import qualified Data.ByteString.Char8         as B
import           Data.Function                 ((&))
import           Data.List
import           Data.Maybe
import           Data.Time                     (getCurrentTime)
import           Lens.Family2                  (view)
import           Pipes
import qualified Pipes.Group                   as PG
import qualified Pipes.Prelude                 as P
import           Pipes.Safe                    (runSafeT)
import           Poseidon.EntityTypes          (HasNameAndVersion (..))
import           Poseidon.GenotypeData
import           Poseidon.Janno
import           Poseidon.Package
import           Poseidon.Utils
import           SequenceFormats.Eigenstrat
import           SequenceFormats.Plink         (eigenstratInd2PlinkFam,
                                                writePlink)
import           System.Directory              (createDirectoryIfMissing)
import           System.FilePath               (takeBaseName, (<.>), (</>))

data AdmixPopsMethodSettings =
    PerSNP {
        _admixMarginalizeMissing :: Bool
    } |
    InChunks {
        _admixChunkSize :: Int
    }

data AdmixPopsOptions = AdmixPopsOptions {
      _admixGenoSources             :: [GenoDataSource]
    , _admixIndWithAdmixtureSet     :: [RequestedInd]
    , _admixIndWithAdmixtureSetFile :: Maybe FilePath
    , _admixMethodSettings          :: AdmixPopsMethodSettings
    , _admixOutFormat               :: GenotypeOutFormatSpec
    , _admixOutPath                 :: FilePath
    , _admixOutPacName              :: Maybe String
    , _admixOutputPlinkPopMode      :: PlinkPopNameMode
    , _admixOnlyLatest              :: Bool
    , _admixOutZip                  :: Bool
    }

runAdmixPops :: AdmixPopsOptions -> PoseidonIO ()
runAdmixPops (
    AdmixPopsOptions
        genoSources
        popsWithFracsDirect
        popsWithFracsFile
        methodSetting
        outFormat
        outPath
        maybeOutName
        outPlinkPopMode
        onlyLatest
        outZip
    ) = do
    let pacReadOpts = defaultPackageReadOptions {
      _readOptIgnoreChecksums  = True
    , _readOptIgnoreGeno       = False
    , _readOptGenoCheck        = True
    , _readOptOnlyLatest       = onlyLatest
    }
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
    properPackages <- readPoseidonPackageCollection pacReadOpts $ [getPacBaseDir x | x@PacBaseDir {} <- genoSources]
    pseudoPackages <- mapM makePseudoPackageFromGenotypeData $ [getGenoDirect x | x@GenoDirect {} <- genoSources]
    logInfo $ "Unpackaged genotype data files loaded: " ++ show (length pseudoPackages)
    let allPackages = properPackages ++ pseudoPackages
    -- determine relevant packages
    let popNames = concatMap (map _inPopName . _inPopSet) requestedInds
    relevantPackages <- filterPackagesByPops popNames allPackages
    -- gather additional info for requested inds
    preparedInds <- mapM (`gatherInfoForInd` relevantPackages) requestedInds
    -- create new package
    let outName = case maybeOutName of -- take basename of outPath, if name is not provided
            Just x  -> x
            Nothing -> takeBaseName outPath
    when (outName == "") $ liftIO $ throwIO PoseidonEmptyOutPacNameException
    -- create new directory
    logInfo $ "Writing to directory (will be created if missing): " ++ outPath
    liftIO $ createDirectoryIfMissing True outPath
    -- compile genotype data structure
    let gz = if outZip then "gz" else ""
        gFileSpec = case outFormat of
            GenotypeOutFormatEigenstrat -> GenotypeEigenstrat (outName <.> ".geno" <.> gz) Nothing (outName <.> ".snp" <.> gz) Nothing (outName <.> ".ind") Nothing
            GenotypeOutFormatPlink      -> GenotypePlink (outName <.> ".bed" <.> gz) Nothing (outName <.> ".bim" <.> gz) Nothing (outName <.> ".fam") Nothing
            GenotypeOutFormatVCF        -> GenotypeVCF (outName <.> ".vcf" <.> gz) Nothing
    let genotypeData = GenotypeDataSpec gFileSpec Nothing
    pac <- newMinimalPackageTemplate outPath outName genotypeData
    liftIO $ writePoseidonPackage pac
    -- compile genotype data
    logInfo "Compiling individuals"
    logA <- envLogAction
    currentTime <- liftIO getCurrentTime
    errLength <- envErrorLength
    liftIO $ catch (
        runSafeT $ do
            eigenstratProd <- getJointGenotypeData logA False relevantPackages Nothing
            let (outG, outS, outI) = case gFileSpec of
                    GenotypeEigenstrat outGeno _ outSnp _ outInd _ -> (outPath </> outGeno, outPath </> outSnp, outPath </> outInd)
                    GenotypePlink      outGeno _ outSnp _ outInd _ -> (outPath </> outGeno, outPath </> outSnp, outPath </> outInd)
                    GenotypeVCF        outGeno _                   -> (outPath </> outGeno, "", "")
            let newIndEntries = map (\x -> EigenstratIndEntry (B.pack . _indName $ x) Unknown (B.pack . _groupName $ x)) preparedInds
                jannoRows = getJannoRows $ posPacJanno pac
            let outConsumer = case outFormat of
                    GenotypeOutFormatEigenstrat -> writeEigenstrat outG outS outI newIndEntries
                    GenotypeOutFormatPlink      -> writePlink      outG outS outI (map (eigenstratInd2PlinkFam outPlinkPopMode) newIndEntries)
                    GenotypeOutFormatVCF        -> writeVCF        logA jannoRows outG
            case methodSetting of
                PerSNP marginalizeMissing -> do
                    runEffect $ eigenstratProd >->
                        printSNPCopyProgress logA currentTime >->
                        P.mapM (samplePerSNP marginalizeMissing preparedInds) >->
                        outConsumer
                InChunks chunkSize -> do
                    runEffect $ (
                            eigenstratProd &
                            chunkEigenstratByNrSnps chunkSize &
                            PG.maps (samplePerChunk logA preparedInds) &
                            PG.concats
                        ) >->
                        printSNPCopyProgress logA currentTime >->
                        outConsumer
        ) (throwIO . PoseidonGenotypeExceptionForward errLength)
    logInfo "Done"
    where
        chunkEigenstratByNrSnps chunkSize = view (PG.chunksOf chunkSize)

checkIndsWithAdmixtureSets :: [RequestedInd] -> IO ()
checkIndsWithAdmixtureSets requestedInds = do
    checkDuplicateIndNames requestedInds
    mapM_ checkPopFracList requestedInds
    where
        checkDuplicateIndNames :: [RequestedInd] -> IO ()
        checkDuplicateIndNames xs =
            let individualsGrouped = filter (\x -> length x > 1) $ group $ sort $ map _inIndName xs
            in unless (null individualsGrouped) $ do
                throwIO $ PoseidonGeneratorCLIParsingException $
                    "Duplicate individual names: " ++ intercalate "," (nub $ concat individualsGrouped)
        checkPopFracList :: RequestedInd -> IO ()
        checkPopFracList x = do
            let xs = _inPopSet x
                fracs = map _inPopFrac xs
            when (sum fracs /= 1) $ do
                throwIO $ PoseidonGeneratorCLIParsingException $
                    "Fractions in " ++ show x ++ " do not to sum to 100%"

filterPackagesByPops :: [String] -> [PoseidonPackage] -> PoseidonIO [PoseidonPackage]
filterPackagesByPops pops packages = do
    fmap catMaybes . forM packages $ \pac -> do
        inds <- loadIndividuals (posPacBaseDir pac) (posPacGenotypeData pac)
        let groupNamesPac = [B.unpack groupName | EigenstratIndEntry _ _ groupName <- inds]
        if   not (null (groupNamesPac `intersect` pops))
        then return (Just pac)
        else return Nothing

gatherInfoForInd :: RequestedInd -> [PoseidonPackage] -> PoseidonIO IndConcrete
gatherInfoForInd (RequestedInd name_ group_ set_) pacs = do
    inds <- mapM (`extractIndsPerPop` pacs) set_
    return $ IndConcrete name_ group_ inds

extractIndsPerPop :: PopFrac -> [PoseidonPackage] -> PoseidonIO PopFracConcrete
extractIndsPerPop (PopFrac pop_ frac_) relevantPackages = do
    let allPackageNames = map getPacName relevantPackages
    allIndEntries <- mapM (\pac -> loadIndividuals (posPacBaseDir pac) (posPacGenotypeData pac)) relevantPackages
    let filterFunc (_,_,EigenstratIndEntry _ _ _group) = B.unpack _group == pop_
        indNames = map (B.unpack . extractIndName) $ filter filterFunc (zipGroup allPackageNames allIndEntries)
        indIDs = map extractFirst $ filter filterFunc (zipGroup allPackageNames allIndEntries)
    return (PopFracConcrete pop_ frac_ (zip indNames indIDs))
    where
        extractIndName (_,_,EigenstratIndEntry x _ _) = x
