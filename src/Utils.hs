{-# LANGUAGE OverloadedStrings #-}

module Utils (JackknifeMode(..), GenomPos, popSpecsNparser, computeAlleleFreq,
computeJackknife, P.runParser, PopConfig(..), GroupDef, XerxesException(..)) where

import           Control.Exception          (Exception)
import           Control.Monad              (forM)
import           Data.Aeson                 (FromJSON, Object, parseJSON,
                                             withObject, (.:), (.:?))
import           Data.Aeson.Types           (Parser)
import           Data.HashMap.Strict        (toList)
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

popSpecsNparser :: Int -> P.Parser [PoseidonEntity]
popSpecsNparser n = sepByNparser n entitySpecParser (P.char ',' <* P.spaces)

sepByNparser :: Int -> P.Parser a -> P.Parser sep -> P.Parser [a]
sepByNparser 0 _ _ = return []
sepByNparser 1 p _ = fmap (: []) p
sepByNparser n p s = do
    x <- p
    _ <- s
    xs <- sepByNparser (n - 1) p s
    return (x:xs)

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
