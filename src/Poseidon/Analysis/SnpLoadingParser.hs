module Poseidon.Analysis.SnpLoadingParser where

import qualified Data.Attoparsec.ByteString.Char8 as A

data SnpLoadingEntry = SnpLoadingEntry
    { snpId         :: B.ByteString
    , snpChrom      :: Chrom
    , snpPos        :: Int
    , snpLoadings   :: [Double]
    }
    deriving (Eq, Show)

snpLoadingParser :: A.Parser SnpLoadingEntry
snpLoadingParser = do
    snpId_ <- A.skipMany A.space >> word
    chrom <- A.skipMany1 A.space >> word
    pos <- A.skipMany1 A.space >> A.decimal
    snpLoadings <- A.double `A.sepBy1` (A.skipMany1 A.space)
    return $ SnpLoadingEntry snpId snpChrom snpPos snpLoadings
  where
    word = A.takeTill isSpace