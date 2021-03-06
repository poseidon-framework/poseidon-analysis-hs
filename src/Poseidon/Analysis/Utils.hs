{-# LANGUAGE OverloadedStrings #-}

module Poseidon.Analysis.Utils where

import           Control.Applicative        ((<|>))
import           Control.Exception          (Exception)
import           Control.Monad              (forM)
import           Data.Aeson                 ((.:))
import           Data.Aeson.Types           (Object, Parser)
import           Data.Char                  (isSpace)
import           Data.HashMap.Strict        (toList)
import           Data.Text                  (unpack)
import qualified Data.Vector                as V
import           Pipes                      (Pipe, cat)
import qualified Pipes.Prelude              as P
import           Poseidon.EntitiesList      (PoseidonEntity (..),
                                             SignedEntitiesList,
                                             indInfoConformsToEntitySpec)
import           Poseidon.SecondaryTypes    (IndividualInfo (..))
import           SequenceFormats.Eigenstrat (EigenstratSnpEntry (..),
                                             GenoEntry (..), GenoLine)
import           SequenceFormats.Utils      (Chrom)
import qualified Text.Parsec                as P
import qualified Text.Parsec.String         as P

-- | A datatype representing the two options for how to run the Block-Jackknife
data JackknifeMode = JackknifePerN Int
    | JackknifePerChromosome
    deriving (Show)

-- | A helper type to represent a genomic position.
type GenomPos = (Chrom, Int)

type GroupDef = (String, SignedEntitiesList)

addGroupDefs :: [GroupDef] -> [IndividualInfo] -> [IndividualInfo]
addGroupDefs groupDefs indInfoRows = do
    indInfo@(IndividualInfo _ groupNames _) <- indInfoRows
    let additionalGroupNames = do
            (groupName, signedEntityList) <- groupDefs
            True <- return $ indInfoConformsToEntitySpec signedEntityList indInfo
            return groupName
    return $ indInfo {indInfoGroups = groupNames ++ additionalGroupNames}

parseGroupDefsFromJSON :: Object -> Parser [GroupDef]
parseGroupDefsFromJSON obj = forM (toList obj) $ \(key, _) -> do
    entities <- obj .: key
    return (unpack key, entities)

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
computeJackknifeOriginal :: Double -> [Double] -> [Double] -> (Double, Double)
computeJackknifeOriginal theta_n m_j_list theta_j_list =
    let g               = fromIntegral (length m_j_list)
        n               = sum m_j_list
        theta_Jackknife = g * theta_n - sum [(1.0 - m_j / n) * theta_j | (m_j, theta_j) <- zip m_j_list theta_j_list]
        h_j_list        = [n / m_j | m_j <- m_j_list]
        pseudo_values   = [h_j * theta_n  - (h_j - 1.0) * theta_j | (h_j, theta_j) <- zip h_j_list theta_j_list]
        sigma_square    = 1.0 / g * sum [1.0 / (h_j - 1.0) * (pseudo_value - theta_Jackknife) ^ (2 :: Int) | (h_j, pseudo_value) <- zip h_j_list pseudo_values]
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

filterTransitions :: (Monad m) => Bool -> Pipe (EigenstratSnpEntry, a) (EigenstratSnpEntry, a) m r
filterTransitions noTransitions = if noTransitions then
        P.filter (\(EigenstratSnpEntry _ _ _ _ ref alt, _) -> isTransversion ref alt)
    else
        cat
    where
    isTransversion ref alt = not $ isTransition ref alt
    isTransition ref alt =
        ((ref == 'A') && (alt == 'G')) ||
        ((ref == 'G') && (alt == 'A')) ||
        ((ref == 'C') && (alt == 'T')) ||
        ((ref == 'T') && (alt == 'C'))

data XerxesException = PopConfigYamlException FilePath String
    | GroupDefException String
    | FStatException String
    deriving (Show)

instance Exception XerxesException

