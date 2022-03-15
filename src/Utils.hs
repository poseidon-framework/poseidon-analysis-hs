{-# LANGUAGE OverloadedStrings #-}

module Utils (JackknifeMode(..), GenomPos, popSpecsNparser, computeAlleleFreq,
computeJackknife, P.runParser, PopConfig(..), GroupDef, XerxesException(..), readPopConfig,
addGroupDefs, computeAlleleCount) where

import           Control.Applicative        ((<|>))
import           Control.Exception          (Exception, throwIO)
import           Control.Monad              (forM)
import           Data.Aeson                 (FromJSON, Object, parseJSON,
                                             withObject, (.:), (.:?))
import           Data.Aeson.Types           (Parser)
import qualified Data.ByteString            as B
import           Data.Char                  (isSpace)
import           Data.HashMap.Strict        (toList)
import           Data.Vector                ((!))
import           Data.Yaml                  (decodeEither')
import           SequenceFormats.Eigenstrat (GenoEntry (..), GenoLine)
import           SequenceFormats.Utils      (Chrom)

import           Data.Text                  (Text)
import           Poseidon.EntitiesList      (EntitiesList, PoseidonEntity (..),
                                             SignedEntitiesList, entitiesListP,
                                             entitySpecParser, indInfoConformsToEntitySpec)
import           Poseidon.SecondaryTypes    (IndividualInfo (..))
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

computeAlleleFreq :: GenoLine -> [Int] -> Maybe Double
computeAlleleFreq line indices =
    let (nrDerived, nrNonMissing) = computeAlleleCount line indices
    in  if nrNonMissing > 0
        then Just (fromIntegral nrDerived / fromIntegral nrNonMissing)
        else Nothing

computeAlleleCount :: GenoLine -> [Int] -> (Int, Int)
computeAlleleCount line indices =
    let nrNonMissing = length . filter (/=Missing) . map (line !) $ indices
        nrDerived = sum $ do
            i <- indices
            case line ! i of
                HomRef  -> return (0 :: Int)
                Het     -> return 1
                HomAlt  -> return 2
                Missing -> return 0
    in  (nrDerived, 2 * nrNonMissing)

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

readPopConfig :: FilePath -> IO PopConfig
readPopConfig fn = do
    bs <- B.readFile fn
    case decodeEither' bs of
        Left err -> throwIO $ PopConfigYamlException fn (show err)
        Right x  -> return x

addGroupDefs :: [GroupDef] -> [IndividualInfo] -> [IndividualInfo]
addGroupDefs groupDefs indInfoRows = do
    indInfo@(IndividualInfo _ groupNames _) <- indInfoRows
    let additionalGroupNames = do
            (groupName, signedEntityList) <- groupDefs
            True <- return $ indInfoConformsToEntitySpec signedEntityList indInfo
            return groupName
    return $ indInfo {indInfoGroups = groupNames ++ additionalGroupNames}
