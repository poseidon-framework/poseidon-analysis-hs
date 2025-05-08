module Poseidon.Analysis.CLI.PCproject where

import           Poseidon.Analysis.Utils (JackknifeMode)

import Control.Monad (filterM)
import           Poseidon.EntityTypes    (EntityInput, PoseidonEntity (..),
                                          readEntityInputs, isLatestInCollection, makePacNameAndVersion)
import           Poseidon.GenotypeData   (GenoDataSource (..))
import           Poseidon.Package        (PackageReadOptions (..),
                                          defaultPackageReadOptions,
                                          makePseudoPackageFromGenotypeData,
                                          readPoseidonPackageCollection)
import           Poseidon.Utils          (PoseidonIO, logInfo)

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
            logInfo "No requested entities. Implicitly forging all packages (latest versions)."
            return $ map (Pac . makePacNameAndVersion) allLatestPackages
        e  -> return e

    logInfo $ "Called pcProject with " ++ show entities
