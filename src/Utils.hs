module Utils (JackknifeMode(..), PopSpec(..), GenomPos, popSpecsNparser, computeAlleleFreq,
computeJackknife, getPopIndices, popSpecParser, P.runParser) where

import           SequenceFormats.Eigenstrat (GenoEntry (..), GenoLine, EigenstratIndEntry(..))
import           SequenceFormats.Utils      (Chrom)

import           Control.Applicative        ((<|>))
import           Data.Char                  (isSpace)
import           Data.Vector                ((!))
import           Poseidon.Utils             (PoseidonException (..))
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

popSpecsNparser :: Int -> P.Parser [PopSpec]
popSpecsNparser n = sepByNparser n popSpecParser (P.char ',' <* P.spaces)

sepByNparser :: Int -> P.Parser a -> P.Parser sep -> P.Parser [a]
sepByNparser 0 _ _ = return []
sepByNparser 1 p _ = fmap (: []) p
sepByNparser n p s = do
    x <- p
    _ <- s
    xs <- sepByNparser (n - 1) p s
    return (x:xs)

popSpecParser :: P.Parser PopSpec
popSpecParser = parseIndividualSpec <|> parseGroupSpec
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

getPopIndices :: [EigenstratIndEntry] -> PopSpec -> Either PoseidonException [Int]
getPopIndices indEntries popSpec =
    let ret = do
            (i, EigenstratIndEntry indName _ popName) <- zip [0..] indEntries
            True <- case popSpec of
                PopSpecGroup name -> return (name == popName)
                PopSpecInd   name -> return (name == indName)
            return i
    in  if null ret
        then case popSpec of
            PopSpecGroup n -> Left $ PoseidonIndSearchException ("Group name " ++ n ++ " not found")
            PopSpecInd   n -> Left $ PoseidonIndSearchException ("Individual name " ++ n ++ " not found")
        else Right ret
