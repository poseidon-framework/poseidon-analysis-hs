module Utils (JackknifeMode(..), PopSpec(..), GenomPos, popSpecsNparser, computeAlleleFreq,
computeJackknife, getPopIndices, popSpecParser, P.runParser, findRelevantPackages) where

import           SequenceFormats.Eigenstrat (EigenstratIndEntry (..),
                                             GenoEntry (..), GenoLine)
import           SequenceFormats.Utils      (Chrom)

import           Control.Applicative        ((<|>))
import           Control.Monad              (forM)
import           Data.Char                  (isSpace)
import           Data.List                  (intersect)
import           Data.Maybe                 (catMaybes)
import           Data.Vector                ((!))
import           Poseidon.Janno             (JannoRow (..), JannoList(..))
import           Poseidon.Package           (PoseidonPackage (..),
                                             getIndividuals)
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

getPopIndices :: [JannoRow] -> PopSpec -> Either PoseidonException [Int]
getPopIndices jointJanno popSpec =
    let ret = do
            (i, jannoRow) <- zip [0..] jointJanno
            True <- case popSpec of
                PopSpecGroup name -> return (name `elem` (getJannoList . jGroupName $ jannoRow))
                PopSpecInd   name -> return (name == jPoseidonID jannoRow)
            return i
    in  if null ret
        then case popSpec of
            PopSpecGroup n -> Left $ PoseidonIndSearchException ("Group name " ++ n ++ " not found")
            PopSpecInd   n -> Left $ PoseidonIndSearchException ("Individual name " ++ n ++ " not found")
        else Right ret

findRelevantPackages :: [PopSpec] -> [PoseidonPackage] -> IO [PoseidonPackage]
findRelevantPackages popSpecs packages = do
    let indNamesStats   = [ind   | PopSpecInd   ind   <- popSpecs]
        groupNamesStats = [group | PopSpecGroup group <- popSpecs]
    fmap catMaybes . forM packages $ \pac -> do
        let janno = posPacJanno pac
        let indNamesPac   = map jPoseidonID janno
            groupNamesPac = concatMap (getJannoList . jGroupName) janno
        if   not (null (indNamesPac `intersect` indNamesStats)) || not (null (groupNamesPac `intersect` groupNamesStats))
        then return (Just pac)
        else return Nothing
