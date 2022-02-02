module Utils (JackknifeMode(..), PopSpec(..), GenomPos, parsePopSpecsN, computeAlleleFreq, computeJackknife) where

import           SequenceFormats.Eigenstrat (GenoEntry (..), GenoLine)
import           SequenceFormats.Utils      (Chrom)

import           Control.Applicative        ((<|>))
import           Data.Char                  (isSpace)
import           Data.Vector                ((!))
import qualified Text.Parsec                as P
import qualified Text.Parsec.String         as P

-- | A datatype representing the two options for how to run the Block-Jackknife
data JackknifeMode = JackknifePerN Int
    | JackknifePerChromosome
    deriving (Show)

-- | A datatype to represent a group or an individual
data PopSpec = PopSpecGroup String
    | PopSpecInd String
    deriving (Eq)

instance Show PopSpec where
    show (PopSpecGroup n) = n
    show (PopSpecInd   n) = "<" ++ n ++ ">"

-- | A helper type to represent a genomic position.
type GenomPos = (Chrom, Int)

data BlockData = BlockData
    { blockStartPos  :: GenomPos
    , blockEndPos    :: GenomPos
    , blockSiteCount :: Int
    , blockVal       :: Double
    }
    deriving (Show)

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


computeAlleleFreq :: GenoLine -> [Int] -> Maybe Double
computeAlleleFreq line indices =
    let nrNonMissing = length . filter (/=Missing) . map (line !) $ indices
        nrDerived = sum $ do
            i <- indices
            case line ! i of
                HomRef  -> return (0 :: Integer)
                Het     -> return 1
                HomAlt  -> return 2
                Missing -> return 0
    in  if nrNonMissing > 0
        then Just (fromIntegral nrDerived / fromIntegral nrNonMissing)
        else Nothing

computeJackknife :: [Int] -> [Double] -> (Double, Double)
computeJackknife weights values =
    let weights'    = map fromIntegral weights
        sumWeights  = sum weights'
        g           = fromIntegral (length weights)
        theta       = sum [mj * val | (mj, val) <- zip weights' values] / sumWeights
        sigmaSquare = sum [mj * (val - theta) ^ (2 :: Int) / (sumWeights - mj) | (mj, val) <- zip weights' values] / g
    in  (theta, sqrt sigmaSquare)
