module Poseidon.Generator.CLI.AdmixPops where

import           Poseidon.Generator.Parsers
import           Poseidon.Generator.SampleGeno
import           Poseidon.Generator.Types
import           Poseidon.Generator.Utils

import           Colog.Message                 (logInfo)
import           Control.Exception             (throwIO)
import           Control.Monad                 (forM, unless, when)
import           Data.List
import           Data.Maybe
import           Data.Ratio                    ((%))
import           Data.Text                     (pack)
import           Pipes
import qualified Pipes.Prelude                 as P
import           Pipes.Safe                    (runSafeT)
import           Poseidon.GenotypeData
import           Poseidon.Janno
import           Poseidon.Package
import           Poseidon.Utils
import           SequenceFormats.Eigenstrat    (EigenstratIndEntry (..),
                                                writeEigenstrat)
import           SequenceFormats.Plink         (writePlink)
import           System.Directory              (createDirectoryIfMissing)
import           System.FilePath               ((</>))
import           System.IO                     (hPrint, hPutStrLn, stderr)

data AdmixPopsOptions = AdmixPopsOptions {
      _admixBaseDirs                :: [FilePath]
    , _admixIndWithAdmixtureSet     :: [IndWithAdmixtureSet]
    , _admixIndWithAdmixtureSetFile :: Maybe FilePath
    , _admixMarginalizeMissing      :: Bool
    , _admixOutFormat               :: GenotypeFormatSpec
    , _admixOutPath                 :: FilePath
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
runAdmixPops (AdmixPopsOptions baseDirs popsWithFracsDirect popsWithFracsFile marginalizeMissing outFormat outDir) = do
    -- compile individuals
    popsWithFracsFromFile <- case popsWithFracsFile of
        Nothing -> return []
        Just f  -> liftIO $ readIndWithAdmixtureSetFromFile f
    let requestedInds = popsWithFracsDirect ++ popsWithFracsFromFile
    -- validating input
    logInfo $ pack "Checking chimeras"
    logInfo $ pack $ renderRequestedInds requestedInds
    liftIO $ checkIndsWithAdmixtureSets requestedInds
    -- load Poseidon packages
    allPackages <- readPoseidonPackageCollection pacReadOpts baseDirs
    -- determine relevant packages and indices
    let popsWithFracs = map (_popFracList . _admixSet) requestedInds
        pops = map (map pop) popsWithFracs
    relevantPackages <- liftIO $ filterPackagesByPops (concat pops) allPackages
    popsFracsInds <- liftIO $  mapM (mapM (`extractIndsPerPop` relevantPackages)) popsWithFracs
    -- compile genotype data structure
    let [outInd, outSnp, outGeno] = case outFormat of
            GenotypeFormatEigenstrat -> ["admixpops_package.ind", "admixpops_package.snp", "admixpops_package.geno"]
            GenotypeFormatPlink -> ["admixpops_package.fam", "admixpops_package.bim", "admixpops_package.bed"]
    -- create output poseidon package
    logInfo $ pack "Creating output Poseidon package"
    liftIO $ createDirectoryIfMissing True outDir
    let genotypeData = GenotypeDataSpec outFormat outGeno Nothing outSnp Nothing outInd Nothing Nothing
        pac = newMinimalPackageTemplate outDir "admixpops_package" genotypeData
    liftIO $ writePoseidonPackage pac
    -- compile genotype data
    logInfo $ pack "Compiling chimeras"
    liftIO $ runSafeT $ do
        (_, eigenstratProd) <- getJointGenotypeData SimpleLog False relevantPackages Nothing
        let [outG, outS, outI] = map (outDir </>) [outGeno, outSnp, outInd]
            newIndEntries = map (\x -> EigenstratIndEntry (_admixInd x) Unknown (_admixUnit x)) requestedInds
        let outConsumer = case outFormat of
                GenotypeFormatEigenstrat -> writeEigenstrat outG outS outI newIndEntries
                GenotypeFormatPlink      -> writePlink      outG outS outI newIndEntries
        runEffect $ eigenstratProd >->
            printSNPCopyProgress >->
            P.mapM (sampleGenoForMultipleIndWithAdmixtureSet marginalizeMissing popsFracsInds) >->
            outConsumer
    logInfo $ pack "Done"

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
        checkDuplicateIndNames requestedInds =
            let individualsGrouped = filter (\x -> length x > 1) $ group $ sort $ map _admixInd requestedInds
            in unless (null individualsGrouped) $ do
                throwIO $ PoseidonGeneratorCLIParsingException $
                    "Duplicate individual names: " ++ intercalate "," (nub $ concat individualsGrouped)
        checkPopFracList :: IndWithAdmixtureSet -> IO ()
        checkPopFracList cur = do
            let xs = (_popFracList . _admixSet) cur
                fracs = map frac xs
            when (sum fracs /= 100) $ do
                throwIO $ PoseidonGeneratorCLIParsingException $
                    "Fractions in " ++ show cur ++ " do not to sum to 100%"

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
    return (map extractFirst $ filter filterFunc (zipGroup allPackageNames allIndEntries), toInteger _frac % 100)
