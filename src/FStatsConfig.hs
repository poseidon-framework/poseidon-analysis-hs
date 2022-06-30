{-# LANGUAGE OverloadedStrings #-}
module FStatsConfig where

import Utils (GroupDef)
import Poseidon.EntitiesList (PoseidonEntity(..), SignedEntity(..))
-- import qualified Dhall as D

-- | A datatype to represent different types of F-Statistics
data FStatType = F4 | F3 | F2 | PWM | Het | FST | F3vanilla | F2vanilla | FSTvanilla
    deriving (Show, Read, Eq)

-- | A datatype to represent F-Statistics to be computed from genotype data.
data FStatSpec = FStatSpec FStatType [PoseidonEntity] (Maybe AscertainmentSpec) deriving (Eq, Show)

data AscertainmentSpec = AscertainmentSpec {
    ascOutgroupSpec  :: Maybe PoseidonEntity, -- Here, a Nothing denotes that minor allele frequency in ascRefgroup should be used, 
                                         --otherwise use derived allele frequencies in the ref-group, where the ancestral allele is given by the outgroup
    ascRefgroupSpec  :: PoseidonEntity,
    ascLowerFreqSpec :: Double,
    ascUpperFreqSpec :: Double
} deriving (Show, Eq)

data FStatConfigYaml = FStatConfigYaml {
    fcGroupDefs :: [GroupDef],
    fcFStatSpecs :: [MultiFStatSpec]
} deriving (Show)

data MultiFStatSpec = MultiFStatSpec {
    mFstatType          :: FStatType,
    mFstatSlotA         :: [PoseidonEntity],
    mFstatSlotB         :: [PoseidonEntity],
    mFstatSlotC         :: [PoseidonEntity],
    mFstatSlotD         :: [PoseidonEntity],
    mFstatAscertainment :: Maybe AscertainmentSpec
} deriving (Show)

instance FromJSON FStatConfigYaml where
    parseJSON = withObject "FStatConfigYaml" $ \v -> FStatConfigYaml
        <$> parseGroupDefsFromJSON v
        <*> parseMultiFStatSpecs v

instance FromJSON MultiFStatSpec where
    parseJSON = withObject ""

-- readFstatConfigDhall :: FilePath -> IO FStatConfig
-- readFstatConfigDhall = D.inputFile fStatD

-- fStatD :: D.Decoder FStatConfig
-- fStatD = D.record (
--     FStatConfig <$> D.field "groupDefs" (D.list groupDefD)
--                 <*> D.field "fStatSpecs" (D.list fStatSpecD))

-- groupDefD :: D.Decoder GroupDef
-- groupDefD = D.record (
--     (,) <$> D.field "name" D.string <*> D.field "entities" (D.list signedEntityD))
--   where
--     signedEntityD = (Include . Group) <$> D.string

-- fStatSpecD :: D.Decoder FStatSpec 
-- fStatSpecD = D.record (
--     FStatSpec <$> D.field "type" fTypeD
--               <*> D.field "slots" (D.list entityD)
--               <*> D.field "ascertainment" (D.maybe ascertainmentSpecD))
--   where
--     entityD = Group <$> D.string
--     fTypeD = (const F3) <$> D.string

-- ascertainmentSpecD :: D.Decoder AscertainmentSpec
-- ascertainmentSpecD = D.record (
--     AscertainmentSpec <$> D.field "outgroup" (D.maybe (Group <$> D.string))
--                       <*> D.field "refgroup" (Group <$> D.string)
--                       <*> D.field "lofreq" D.double
--                       <*> D.field "hifreq" D.double)
