module Poseidon.Analysis.CLI.PCproject where

import           Poseidon.Analysis.SnpLoadingParser (snpLoadingParser)
import           Poseidon.Analysis.Utils            (JackknifeMode)

import           Control.Monad                      (filterM)
import           Poseidon.EntityTypes               (EntityInput,
                                                     PoseidonEntity (..),
                                                     isLatestInCollection,
                                                     makePacNameAndVersion,
                                                     readEntityInputs)
import           Poseidon.GenotypeData              (GenoDataSource (..))
import           Poseidon.Package                   (PackageReadOptions (..),
                                                     defaultPackageReadOptions,
                                                     makePseudoPackageFromGenotypeData,
                                                     readPoseidonPackageCollection)
import           Poseidon.Utils                     (PoseidonIO, logInfo)

data PCprojectOpts = PCprojectOpts
    { _pcGenoSource      :: [GenoDataSource]
    , _pcSnpLoadingsFile :: FilePath
    , _pcJackknifeMode   :: JackknifeMode
    , _pcEntities        :: [EntityInput PoseidonEntity]
}

pacReadOpts :: PackageReadOptions
pacReadOpts = defaultPackageReadOptions {
      _readOptIgnoreChecksums  = True
    , _readOptIgnoreGeno       = False
    , _readOptGenoCheck        = True
    }

runPCproject :: PCprojectOpts -> PoseidonIO ()
runPCproject (PCprojectOpts genoSources snpLoadingFile jackknifeMode entityInputs) = do

    properPackages <- readPoseidonPackageCollection pacReadOpts $ [getPacBaseDir x | x@PacBaseDir {} <- genoSources]
    pseudoPackages <- mapM makePseudoPackageFromGenotypeData [getGenoDirect x | x@GenoDirect {} <- genoSources]
    logInfo $ "Unpackaged genotype data files loaded: " ++ show (length pseudoPackages)
    let allPackages = properPackages ++ pseudoPackages

    -- compile entities
    entitiesUser <- readEntityInputs entityInputs

    allLatestPackages <- filterM (isLatestInCollection allPackages) allPackages

    entities <- case entitiesUser of
        [] -> do
            logInfo "No requested entities. Implicitly selecting all individuals from all latest packages."
            return $ map (Pac . makePacNameAndVersion) allLatestPackages
        e  -> return e

    relevantPackages <- filterToRelevantPackages entities allPackages
    logInfo $ (show . length $ relevantPackages) ++ " packages contain data for this forging operation"
    when (null relevantPackages) $ liftIO $ throwIO PoseidonEmptyForgeException

    logInfo $ "Reading Snp Loadings from file " ++ snpLoadingFile

    runSafeT $ do
        snpLoadingProd <- readSnpLoadingFile snpLoadingFile
        eigenstratProd <- getJointGenotypeData logA False relevantPackages Nothing

