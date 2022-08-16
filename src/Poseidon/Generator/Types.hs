module Poseidon.Generator.Types where

import           Data.List      (intercalate)
import           Poseidon.Janno

data IndAdmixpops = IndAdmixpops {
      _indName   :: String
    , _groupName :: String
    , _popSet    :: [PopAdmixpops]
} deriving (Show)

data PopAdmixpops = PopAdmixpops {
      _popName :: String
    , _popFrac :: Rational
    , _popInds :: [(String, Int)]
} deriving (Show)

data InIndAdmixpops = InIndAdmixpops {
      _inIndName   :: String
    , _inGroupName :: String
    , _inPopSet    :: [InPopAdmixpops]
}

instance Show InIndAdmixpops where
    show (InIndAdmixpops _admixInd _admixUnit _popFracList) =
        "[" ++ _admixInd ++ ":" ++ _admixUnit ++ "]" ++
        "(" ++ intercalate "+" (map show _popFracList) ++ ")"

data InPopAdmixpops = InPopAdmixpops {
      _inPopName :: String
    , _inPopFrac :: Rational
}

instance Show InPopAdmixpops where
    show (InPopAdmixpops _pop _frac) =
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
