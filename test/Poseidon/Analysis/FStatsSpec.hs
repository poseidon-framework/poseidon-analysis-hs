module Poseidon.Analysis.FStatsSpec (spec) where

import Poseidon.EntitiesList (PoseidonEntity(..))
import Poseidon.Analysis.CLI.FStats (collectStatSpecGroups)
import Poseidon.Analysis.FStatsConfig (FStatSpec(..), FStatType(..), AscertainmentSpec(..))

import Test.Hspec

spec :: Spec
spec = do
    testCollectStats

testCollectStats :: Spec
testCollectStats = describe "collectStatSpecGroups" $ do
    it "should correctly collect stats" $ do
        let statSpecs = [
                FStatSpec F3 [Group "French", Group "Spanish", Ind "Chimp.REF"] Nothing,
                FStatSpec F3 [Group "French", Group "Mbuti", Ind "Chimp.REF"] (Just (AscertainmentSpec (Just (Ind "Human.REF")) (Group "CEU") 0.05 0.95))]
            entities = [Group "French", Group "Spanish", Ind "Chimp.REF", Group "Mbuti", Ind "Human.REF", Group "CEU"]
        collectStatSpecGroups statSpecs `shouldMatchList` entities
