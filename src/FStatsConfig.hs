module FStatsConfig where

import Utils (GroupDef)
import Poseidon.EntityList (PoseidonEntity)

data FStatConfig = FStatConfig {
    fcGroupDefs :: [GroupDef],
    fcFStatSpecs :: [FStatSpec]
}

-- | A datatype to represent F-Statistics to be computed from genotype data.
data FStatSpec = FStatSpec FStatType [PoseidonEntity] (Maybe AscertainmentSpec) deriving (Eq, Show)

data AscertainmentSpec = AscertainmentSpec {
    ascOutgroupSpec  :: Maybe PoseidonEntity, -- Here, a Nothing denotes that minor allele frequency in ascRefgroup should be used, 
                                         --otherwise use derived allele frequencies in the ref-group, where the ancestral allele is given by the outgroup
    ascRefgroupSpec  :: PoseidonEntity,
    ascLowerFreqSpec :: Double,
    ascUpperFreqSpec :: Double
} deriving (Show, Eq)
