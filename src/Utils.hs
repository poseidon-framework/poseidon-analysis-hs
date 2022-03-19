{-# LANGUAGE OverloadedStrings #-}

module Utils (JackknifeMode(..), GenomPos, popSpecsNparser, computeAlleleFreq, computeAlleleCount,
computeJackknifeAdditive, P.runParser, PopConfig(..), GroupDef, XerxesException(..)) where

import           Control.Applicative        ((<|>))
import           Control.Exception          (Exception)
import           Control.Monad              (forM)
import           Data.Aeson                 (FromJSON, Object, parseJSON,
                                             withObject, (.:), (.:?))
import           Data.Aeson.Types           (Parser)
import           Data.Char                  (isSpace)
import           Data.HashMap.Strict        (toList)
import qualified Data.Vector                as V
import           SequenceFormats.Eigenstrat (GenoEntry (..), GenoLine)
import           SequenceFormats.Utils      (Chrom)

import           Data.Text                  (Text)
import           Data.Vector                ((!))
import           Poseidon.EntitiesList      (EntitiesList, PoseidonEntity (..),
                                             SignedEntitiesList, entitiesListP,
                                             entitySpecParser)
import qualified Text.Parsec                as P
import qualified Text.Parsec.String         as P

-- | A datatype representing the two options for how to run the Block-Jackknife
data JackknifeMode = JackknifePerN Int
    | JackknifePerChromosome
    deriving (Show)

-- | A helper type to represent a genomic position.
type GenomPos = (Chrom, Int)

customEntitySpecParser :: P.Parser PoseidonEntity
customEntitySpecParser = parsePac <|> parseGroup <|> parseInd
    where
    parsePac   = Pac   <$> P.between (P.char '*') (P.char '*') parseName
    parseGroup = Group <$> parseName
    parseInd   = Ind   <$> P.between (P.char '<') (P.char '>') parseName
    charList :: [Char]
    charList = ",<>*()"
    parseName  = P.many1 (P.satisfy (\c -> not (isSpace c || c `elem` charList)))


popSpecsNparser :: Int -> P.Parser [PoseidonEntity]
popSpecsNparser n = sepByNparser n customEntitySpecParser (P.char ',' <* P.spaces)

sepByNparser :: Int -> P.Parser a -> P.Parser sep -> P.Parser [a]
sepByNparser 0 _ _ = return []
sepByNparser 1 p _ = fmap (: []) p
sepByNparser n p s = do
    x <- p
    _ <- s
    xs <- sepByNparser (n - 1) p s
    return (x:xs)

computeAlleleCount :: GenoLine -> [Int] -> (Int, Int)
computeAlleleCount line indices =
    let nrNonMissing = length . filter (/=Missing) . map (line V.!) $ indices
        nrDerived = sum $ do
            i <- indices
            case line V.! i of
                HomRef  -> return (0 :: Int)
                Het     -> return 1
                HomAlt  -> return 2
                Missing -> return 0
    in  (nrDerived, 2 * nrNonMissing)

computeAlleleFreq :: GenoLine -> [Int] -> Maybe Double
computeAlleleFreq line indices =
    let (nrDerivedHaps, nrNonMissingHaps) = computeAlleleCount line indices
    in  if nrNonMissingHaps > 0
        then Just (fromIntegral nrDerivedHaps / fromIntegral nrNonMissingHaps)
        else Nothing

-- ^ This is the original Delete-mj Jackknife based on Busing et al. 1999. The notation is the same.
-- theta_n is the full estimate. m_j is the weights of block j. theta_j is the
-- estimate based on the j'th block removed.
computeJackknifeOriginal :: Double -> [Int] -> [Double] -> (Double, Double)
computeJackknifeOriginal theta_n m_j_list_Int theta_j_list =
    let g               = fromIntegral (length m_j_list_Int)
        m_j_list        = map fromIntegral m_j_list_Int
        n               = sum m_j_list
        theta_Jackknife = g * theta_n - sum [(1.0 - m_j / n) * theta_j | (m_j, theta_j) <- zip m_j_list theta_j_list]
        h_j_list        = [n / m_j | m_j <- m_j_list]
        pseudo_values   = [h_j * theta_n  - (h_j - 1.0) * theta_j | (h_j, theta_j) <- zip h_j_list theta_j_list]
        sigma_square    = 1.0 / g * sum [1.0 / (h_j - 1.0) * (pseudo_value - theta_Jackknife) ^ 2 | (h_j, pseudo_value) <- zip h_j_list pseudo_values]
    in  (theta_Jackknife, sqrt sigma_square)

-- ^ A simplified Jackknife formula based on an additive estimate, in which the full estimate equals the Jackknife estimate
computeJackknifeAdditive :: [Int] -> [Double] -> (Double, Double)
computeJackknifeAdditive weights values =
    let weights'    = map fromIntegral weights
        sumWeights  = sum weights'
        g           = fromIntegral (length weights)
        theta       = sum [mj * val | (mj, val) <- zip weights' values] / sumWeights
        sigmaSquare = sum [mj * (val - theta) ^ (2 :: Int) / (sumWeights - mj) | (mj, val) <- zip weights' values] / g
    in  (theta, sqrt sigmaSquare)



type GroupDef = (String, SignedEntitiesList)

data PopConfig = PopConfigYamlStruct
    { popConfigGroupDef :: [GroupDef]
    , popConfigLefts    :: EntitiesList
    , popConfigRights   :: EntitiesList
    , popConfigOutgroup :: Maybe PoseidonEntity
    }

instance FromJSON PopConfig where
    parseJSON = withObject "PopConfigYamlStruct" $ \v -> PopConfigYamlStruct
        <$> parseGroupDefsFromJSON v
        <*> parsePopSpecsFromJSON v "popLefts"
        <*> parsePopSpecsFromJSON v "popRights"
        <*> parseMaybePopSpecFromJSON v "outgroup"
      where
        parseGroupDefsFromJSON :: Object -> Parser [GroupDef]
        parseGroupDefsFromJSON v = do
            maybeObj <- v .:? "groupDefs"
            case maybeObj of
                Nothing -> return []
                Just obj -> return $ do
                    (key, value) <- toList obj
                    case P.runParser entitiesListP () "" value of
                        Left err -> fail (show err)
                        Right p  -> return (key, p)
        parsePopSpecsFromJSON :: Object -> Text -> Parser [PoseidonEntity]
        parsePopSpecsFromJSON v label = do
            popDefStrings <- v .: label
            forM popDefStrings $ \popDefString -> do
                case P.runParser entitySpecParser () "" popDefString of
                    Left err -> fail (show err)
                    Right p  -> return p
        parseMaybePopSpecFromJSON :: Object -> Text -> Parser (Maybe PoseidonEntity)
        parseMaybePopSpecFromJSON v label = do
            maybePopDefString <- v .:? label
            case maybePopDefString of
                Nothing -> return Nothing
                Just p -> case P.runParser entitySpecParser () "" p of
                    Left err -> fail (show err)
                    Right p' -> return (Just p')

data XerxesException = PopConfigYamlException FilePath String
    | GroupDefException String
    deriving (Show)

instance Exception XerxesException
