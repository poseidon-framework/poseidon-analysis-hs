module Utils (JackknifeMode(..), PopSpec(..), GenomPos, parsePopSpecsN) where

import           SequenceFormats.Utils      (Chrom)

import Control.Applicative ((<|>))
import           Data.Char                  (isSpace)
import qualified Text.Parsec as P
import qualified Text.Parsec.String         as P

-- | A datatype representing the two options for how to run the Block-Jackknife
data JackknifeMode = JackknifePerN Int -- ^ Run the Jackknife in blocks of N consecutive SNPs
    | JackknifePerChromosome -- ^ Run the Jackknife in blocks defined as entire chromosomes.
  deriving (Show)

-- | A datatype to represent a group or an individual
data PopSpec = PopSpecGroup String -- ^ Define a group name
    | PopSpecInd String -- ^ Define an individual name
    deriving (Eq)

instance Show PopSpec where
    show (PopSpecGroup n) = n
    show (PopSpecInd   n) = "<" ++ n ++ ">"

-- | A helper type to represent a genomic position.
type GenomPos = (Chrom, Int)

parsePopSpecsN :: Int -> P.Parser [PopSpec]
parsePopSpecsN n = sepByN n parsePopSpec (P.char ',' <* P.spaces)

sepByN :: Int -> P.Parser a -> P.Parser sep -> P.Parser [a]
sepByN 0 _ _ = return []
sepByN 1 p _ = fmap (: []) p
sepByN n p s = do
    x <- p
    _ <- s
    xs <- sepByN (n - 1) p s
    return (x:xs)

parsePopSpec :: P.Parser PopSpec
parsePopSpec = parseIndividualSpec <|> parseGroupSpec
  where
    parseIndividualSpec = PopSpecInd <$> P.between (P.char '<') (P.char '>') parseName
    parseGroupSpec = PopSpecGroup <$> parseName
    parseName = P.many1 (P.satisfy (\c -> not (isSpace c || c `elem` ",<>()")))


