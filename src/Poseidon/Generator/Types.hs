module Poseidon.Generator.Types (
    IndWithAdmixtureSet (..),
    AdmixtureSet (..),
    PopulationWithFraction (..),
    IndWithPosition (..),
    SpatialTemporalPosition (..),
    GenoEntry (..)
) where

import           Data.List                  (intercalate)
import           Poseidon.Janno
import           SequenceFormats.Eigenstrat (GenoEntry (..))

data IndWithAdmixtureSet = IndWithAdmixtureSet {
      _admixInd  :: String
    , _admixUnit :: String
    , _admixSet  :: AdmixtureSet
}

instance Show IndWithAdmixtureSet where
    show (IndWithAdmixtureSet _admixInd _admixUnit (AdmixtureSet _popFracList)) =
        "[" ++ _admixInd ++ ":" ++ _admixUnit ++ "]" ++
        "(" ++ intercalate "+" (map show _popFracList) ++ ")"

newtype AdmixtureSet = AdmixtureSet {
    _popFracList :: [PopulationWithFraction]
} deriving (Show)

data PopulationWithFraction = PopulationWithFraction {
      pop  :: String
    , frac :: Rational
}

instance Show PopulationWithFraction where
    show (PopulationWithFraction _pop _frac) =
        _pop ++ "=" ++ show _frac

data IndWithPosition = IndWithPosition {
      spatInd  :: String
    , spatUnit :: String
    , spatPos  :: SpatialTemporalPosition
} deriving (Show)

data SpatialTemporalPosition = SpatialTemporalPosition {
      time :: Int
    , lat  :: Latitude
    , lon  :: Longitude
} deriving (Show)
