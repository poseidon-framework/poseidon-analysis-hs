module FStats (
      FStatSpec(..)
    , P.ParseError
    , FstatsOptions(..)
    , JackknifeMode(..)
    , fStatSpecParser
    , runFstats
    , readStatSpecsFromFile
) where

import           Utils                      (GenomPos, JackknifeMode (..),
                                             computeAlleleFreq,
                                             computeJackknife, popSpecsNparser)

import           Control.Applicative        ((<|>))
import           Control.Exception          (throwIO)
import           Control.Foldl              (Fold (..), list, purely)
import           Control.Monad              (forM_)
import           Control.Monad.Catch        (throwM)
import           Data.List                  (intercalate, nub,
                                             transpose)
import           Lens.Family2               (view)
import           Pipes                      ((>->))
import           Pipes.Group                (chunksOf, folds, groupsBy)
import qualified Pipes.Prelude              as P
import           Pipes.Safe                 (runSafeT)
import           Poseidon.EntitiesList      (PoseidonEntity (..),
                                             filterRelevantPackages,
                                             conformingEntityIndices,
                                             findNonExistentEntities)
import           Poseidon.Package           (PackageReadOptions (..),
                                             PoseidonPackage (..),
                                             defaultPackageReadOptions,
                                             getJointGenotypeData,
                                             getJointIndividualInfo,
                                             readPoseidonPackageCollection)
import           Poseidon.SecondaryTypes    (IndividualInfo (..))
import           Poseidon.Utils             (PoseidonException (..))
import           SequenceFormats.Eigenstrat (EigenstratSnpEntry (..), GenoLine)
import           SequenceFormats.Utils      (Chrom)
import           System.IO                  (hPutStrLn, stderr)
import           Text.Layout.Table          (asciiRoundS, column, def, expand,
                                             rowsG, tableString, titlesH)
import qualified Text.Parsec                as P
import qualified Text.Parsec.String         as P

-- | A datatype representing the command line options for the F-Statistics command
data FstatsOptions = FstatsOptions
    { _foBaseDirs        :: [FilePath] -- ^ the list of base directories to search for packages
    -- ^ The way the Jackknife is performed
    , _foJackknifeMode   :: JackknifeMode -- ^ The way the Jackknife is performed
    -- ^ a list of chromosome names to exclude from the computation
    , _foExcludeChroms   :: [Chrom] -- ^ a list of chromosome names to exclude from the computation
    -- ^ A list of F-statistics to compute
    , _foStatSpecsDirect :: [FStatSpec] -- ^ A list of F-statistics to compute
    -- ^ a file listing F-statistics to compute
    , _foStatSpecsFile   :: Maybe FilePath -- ^ a file listing F-statistics to compute
    -- ^ whether to output the result table in raw TSV instead of nicely formatted ASCII table/
    , _foRawOutput       :: Bool -- ^ whether to output the result table in raw TSV instead of nicely formatted ASCII table/
    }

-- | A datatype to represent Summary Statistics to be computed from genotype data.
data FStatSpec = F4Spec PoseidonEntity PoseidonEntity PoseidonEntity PoseidonEntity
    | F3Spec PoseidonEntity PoseidonEntity PoseidonEntity
    | F2Spec PoseidonEntity PoseidonEntity
    | PWMspec PoseidonEntity PoseidonEntity
    deriving (Eq)

instance Show FStatSpec where
    show (F4Spec  a b c d) = "F4("  ++ show a ++ "," ++ show b ++ "," ++ show c ++ "," ++ show d ++ ")"
    show (F3Spec  a b c  ) = "F3("  ++ show a ++ "," ++ show b ++ "," ++ show c ++ ")"
    show (F2Spec  a b    ) = "F2("  ++ show a ++ "," ++ show b ++ ")"
    show (PWMspec a b    ) = "PWM(" ++ show a ++ "," ++ show b ++ ")"

-- | An internal datatype to represent Summary statistics with indices of individuals given as integers
data FStat = F4 [Int] [Int] [Int] [Int]
    | F3 [Int] [Int] [Int]
    | F2 [Int] [Int]
    | PWM [Int] [Int]


data BlockData = BlockData
    { blockStartPos  :: GenomPos
    , blockEndPos    :: GenomPos
    , blockSiteCount :: Int
    , blockVal       :: Double
    }
    deriving (Show)

-- | A parser to parse Summary Statistic specifications.
fStatSpecParser :: P.Parser FStatSpec
fStatSpecParser = P.try f4SpecParser <|> P.try f3SpecParser <|> P.try f2SpecParser <|> pwmSpecParser

-- | A parser to parse F4Stats
f4SpecParser :: P.Parser FStatSpec
f4SpecParser = do
    _ <- P.string "F4"
    [a, b, c, d] <- P.between (P.char '(') (P.char ')') (popSpecsNparser 4)
    return $ F4Spec a b c d


f3SpecParser :: P.Parser FStatSpec
f3SpecParser = do
    _ <- P.string "F3"
    [a, b, c] <- P.between (P.char '(') (P.char ')') (popSpecsNparser 3)
    return $ F3Spec a b c

f2SpecParser :: P.Parser FStatSpec
f2SpecParser = do
    _ <- P.string "F2"
    [a, b] <- P.between (P.char '(') (P.char ')') (popSpecsNparser 2)
    return $ F2Spec a b

pwmSpecParser :: P.Parser FStatSpec
pwmSpecParser = do
    _ <- P.string "PWM"
    [a, b] <- P.between (P.char '(') (P.char ')') (popSpecsNparser 2)
    return $ PWMspec a b

statSpecsFold :: [IndividualInfo] -> [FStatSpec] -> Either PoseidonException (Fold (EigenstratSnpEntry, GenoLine) [BlockData])
statSpecsFold indInfo fStatSpecs = do
    listOfFolds <- mapM (statSpecFold indInfo) fStatSpecs
    return $ sequenceA listOfFolds

statSpecFold :: [IndividualInfo] -> FStatSpec -> Either PoseidonException (Fold (EigenstratSnpEntry, GenoLine) BlockData)
statSpecFold indInfo fStatSpec = do
    let getI = conformingEntityIndices
        i = indInfo
    let fStat = case fStatSpec of
            F4Spec  a b c d -> F4  (getI [a] i) (getI [b] i) (getI [c] i) (getI [d] i)
            F3Spec  a b c   -> F3  (getI [a] i) (getI [b] i) (getI [c] i)
            F2Spec  a b     -> F2  (getI [a] i) (getI [b] i)
            PWMspec a b     -> PWM (getI [a] i) (getI [b] i)
    return $ Fold (step fStat) initialize extract
  where
    step :: FStat -> (Maybe GenomPos, Maybe GenomPos, Int, Double) ->
        (EigenstratSnpEntry, GenoLine) -> (Maybe GenomPos, Maybe GenomPos, Int, Double)
    step fstat (maybeStartPos, _, count, val) (EigenstratSnpEntry c p _ _ _ _, genoLine) =
        let newStartPos = case maybeStartPos of
                Nothing       -> Just (c, p)
                Just (c', p') -> Just (c', p')
            newEndPos = Just (c, p)
        in  case computeFStat fstat genoLine of
                Just v  -> (newStartPos, newEndPos, count + 1, val + v)
                Nothing -> (newStartPos, newEndPos, count + 1, val)
    initialize :: (Maybe GenomPos, Maybe GenomPos, Int, Double)
    initialize = (Nothing, Nothing, 0, 0.0)
    extract :: (Maybe GenomPos, Maybe GenomPos, Int, Double) -> BlockData
    extract (maybeStartPos, maybeEndPos, count, totalVal) =
        let Just startPos = maybeStartPos
            Just endPos = maybeEndPos
        in  BlockData startPos endPos count (totalVal / fromIntegral count)

computeFStat :: FStat -> GenoLine -> Maybe Double
computeFStat fStat gL = case fStat of
    (F4  aI bI cI dI) -> computeF4  <$> computeAlleleFreq gL aI <*> computeAlleleFreq gL bI <*> computeAlleleFreq gL cI <*> computeAlleleFreq gL dI
    (F3  aI bI cI   ) -> computeF3  <$> computeAlleleFreq gL aI <*> computeAlleleFreq gL bI <*> computeAlleleFreq gL cI
    (F2  aI bI      ) -> computeF2  <$> computeAlleleFreq gL aI <*> computeAlleleFreq gL bI
    (PWM aI bI      ) -> computePWM <$> computeAlleleFreq gL aI <*> computeAlleleFreq gL bI
  where
    computeF4  a b c d = (a - b) * (c - d)
    computeF3  a b c   = (c - a) * (c - b)
    computeF2  a b     = (a - b) * (a - b)
    computePWM a b     = a * (1.0 - b) + (1.0 - a) * b

pacReadOpts :: PackageReadOptions
pacReadOpts = defaultPackageReadOptions {
      _readOptVerbose          = False
    , _readOptStopOnDuplicates = True
    , _readOptIgnoreChecksums  = True
    , _readOptIgnoreGeno       = False
    , _readOptGenoCheck        = True
    }

-- | The main function running the FStats command.
runFstats :: FstatsOptions -> IO ()
runFstats (FstatsOptions baseDirs jackknifeMode exclusionList statSpecsDirect maybeStatSpecsFile rawOutput) = do
    -- load packages --
    allPackages <- readPoseidonPackageCollection pacReadOpts baseDirs
    statSpecsFromFile <- case maybeStatSpecsFile of
        Nothing -> return []
        Just f  -> readStatSpecsFromFile f
    let statSpecs = statSpecsFromFile ++ statSpecsDirect
    if null statSpecs then
        hPutStrLn stderr "No statistics to be computed"
    else do
        let collectedStats = collectStatSpecGroups statSpecs
        let missingEntities = findNonExistentEntities collectedStats (getJointIndividualInfo allPackages)
        if not. null $ missingEntities then
            hPutStrLn stderr $ "The following entities couldn't be found: " ++ (intercalate ", " . map show $ missingEntities)
        else do
            let relevantPackages = filterRelevantPackages collectedStats allPackages
            hPutStrLn stderr $ (show . length $ relevantPackages) ++ " relevant packages for chosen statistics identified:"
            forM_ relevantPackages $ \pac -> hPutStrLn stderr (posPacTitle pac)
            hPutStrLn stderr $ "Computing stats " ++ show statSpecs
            blockData <- runSafeT $ do
                (_, eigenstratProd) <- getJointGenotypeData False False relevantPackages Nothing
                let jointIndInfo = getJointIndividualInfo relevantPackages
                let eigenstratProdFiltered = eigenstratProd >-> P.filter chromFilter
                    eigenstratProdInChunks = case jackknifeMode of
                        JackknifePerChromosome  -> chunkEigenstratByChromosome eigenstratProdFiltered
                        JackknifePerN chunkSize -> chunkEigenstratByNrSnps chunkSize eigenstratProdFiltered
                statsFold <- case statSpecsFold jointIndInfo statSpecs of
                    Left e  ->  throwM e
                    Right f -> return f
                let summaryStatsProd = purely folds statsFold eigenstratProdInChunks
                purely P.fold list (summaryStatsProd >-> P.tee (P.map showBlockLogOutput >-> P.toHandle stderr))
            let jackknifeEstimates = [computeJackknife (map blockSiteCount blocks) (map blockVal blocks) | blocks <- transpose blockData]
                colSpecs = replicate 4 (column expand def def def)
                tableH = ["Statistic", "Estimate", "StdErr", "Z score"]
                tableB = do
                    (fstat, result) <- zip statSpecs jackknifeEstimates
                    return [show fstat, show (fst result), show (snd result), show (uncurry (/) result)]
            if   rawOutput
            then do
                putStrLn $ intercalate "\t" tableH
                forM_ tableB $ \row -> putStrLn (intercalate "\t" row)
            else putStrLn $ tableString colSpecs asciiRoundS (titlesH tableH) [rowsG tableB]
  where
    chromFilter (EigenstratSnpEntry chrom _ _ _ _ _, _) = chrom `notElem` exclusionList
    chunkEigenstratByChromosome = view (groupsBy sameChrom)
    sameChrom (EigenstratSnpEntry chrom1 _ _ _ _ _, _) (EigenstratSnpEntry chrom2 _ _ _ _ _, _) =
        chrom1 == chrom2
    chunkEigenstratByNrSnps chunkSize = view (chunksOf chunkSize)
    showBlockLogOutput blocks = "computing chunk range " ++ show (blockStartPos (head blocks)) ++ " - " ++
        show (blockEndPos (head blocks)) ++ ", size " ++ (show . blockSiteCount . head) blocks ++ ", values " ++
        (show . map blockVal) blocks

readStatSpecsFromFile :: FilePath -> IO [FStatSpec]
readStatSpecsFromFile statSpecsFile = do
    let multiFstatSpecParser = fStatSpecParser `P.sepBy1` (P.newline *> P.spaces)
    eitherParseResult <- P.parseFromFile (P.spaces *> multiFstatSpecParser <* P.spaces) statSpecsFile
    case eitherParseResult of
        Left err -> throwIO (PoseidonFStatsFormatException (show err))
        Right r  -> return r

collectStatSpecGroups :: [FStatSpec] -> [PoseidonEntity]
collectStatSpecGroups statSpecs = nub . concat $ do
    stat <- statSpecs
    case stat of
        F4Spec  a b c d -> return [a, b, c, d]
        F3Spec  a b c   -> return [a, b, c]
        F2Spec  a b     -> return [a, b]
        PWMspec a b     -> return [a, b]

