module FStats (
      FStatSpec(..)
    , P.ParseError
    , FstatsOptions(..)
    , JackknifeMode(..)
    , fStatSpecParser
    , runFstats
    , readStatSpecsFromFile
) where

import           Utils                       (GenomPos, JackknifeMode (..),
                                              computeAlleleFreq,
                                              computeJackknife, popSpecsNparser)

import           Control.Applicative         ((<|>))
import           Control.Exception           (throwIO)
import           Control.Foldl               (FoldM (..), impurely, list,
                                              purely)
import           Control.Monad               (forM_)
import           Control.Monad.Catch         (throwM)
import           Control.Monad.IO.Class      (MonadIO, liftIO)
import           Data.List                   (intercalate, nub, transpose)
import qualified Data.Vector.Unboxed         as VU
import qualified Data.Vector.Unboxed.Mutable as VUM
import           Debug.Trace                 (trace)
import           Lens.Family2                (view)
import           Pipes                       (cat, (>->))
import           Pipes.Group                 (chunksOf, foldsM, groupsBy)
import qualified Pipes.Prelude               as P
import           Pipes.Safe                  (runSafeT)
import           Poseidon.EntitiesList       (PoseidonEntity (..),
                                              conformingEntityIndices,
                                              filterRelevantPackages,
                                              findNonExistentEntities)
import           Poseidon.Package            (PackageReadOptions (..),
                                              PoseidonPackage (..),
                                              defaultPackageReadOptions,
                                              getJointGenotypeData,
                                              getJointIndividualInfo,
                                              readPoseidonPackageCollection)
import           Poseidon.SecondaryTypes     (IndividualInfo (..))
import           Poseidon.Utils              (PoseidonException (..))
import           SequenceFormats.Eigenstrat  (EigenstratSnpEntry (..), GenoLine)
import           SequenceFormats.Utils       (Chrom)
import           System.IO                   (hPutStrLn, stderr)
import           Text.Layout.Table           (asciiRoundS, column, def, expand,
                                              rowsG, tableString, titlesH)
import qualified Text.Parsec                 as P
import qualified Text.Parsec.String          as P

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
    , _foMaxSnps         :: Maybe Int
    , _foFullTable         :: Bool
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
    | PWM [Int] [Int] deriving (Show)

data BlockData = BlockData
    { blockStartPos  :: GenomPos
    , blockEndPos    :: GenomPos
    , blockSiteCount :: Int
    , blockStatVal   :: [Double]
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

buildStatSpecsFold :: (MonadIO m) => [IndividualInfo] -> [FStatSpec] -> FoldM m (EigenstratSnpEntry, GenoLine) BlockData
buildStatSpecsFold indInfo fStatSpecs =
    let getI = conformingEntityIndices
        i = indInfo
        fStats = do
            fStatSpec <- fStatSpecs
            case fStatSpec of
                F4Spec  a b c d -> return $ F4  (getI [a] i) (getI [b] i) (getI [c] i) (getI [d] i)
                F3Spec  a b c   -> return $ F3  (getI [a] i) (getI [b] i) (getI [c] i)
                F2Spec  a b     -> return $ F2  (getI [a] i) (getI [b] i)
                PWMspec a b     -> return $ PWM (getI [a] i) (getI [b] i)
    in  FoldM (step fStats) initialize extract
  where
    step :: (MonadIO m) => [FStat] -> (Maybe GenomPos, Maybe GenomPos, Int, VUM.IOVector Int, VUM.IOVector Double) ->
        (EigenstratSnpEntry, GenoLine) -> m (Maybe GenomPos, Maybe GenomPos, Int, VUM.IOVector Int, VUM.IOVector Double)
    step fstats (maybeStartPos, _, count, normVec, valVec) (EigenstratSnpEntry c p _ _ _ _, genoLine) = do
        let newStartPos = case maybeStartPos of
                Nothing       -> Just (c, p)
                Just (c', p') -> Just (c', p')
            newEndPos = Just (c, p)
        forM_ (zip [0..] fstats) $ \(i, fstat) -> do
            case computeFStat fstat genoLine of
                Just v  -> do
                    liftIO $ VUM.modify normVec (+1) i
                    liftIO $ VUM.modify valVec (+v) i
                Nothing -> return ()
        return (newStartPos, newEndPos, count + 1, normVec, valVec)
    initialize :: (MonadIO m) => m (Maybe GenomPos, Maybe GenomPos, Int, VUM.IOVector Int, VUM.IOVector Double)
    initialize = do
        normVec <- liftIO $ VUM.replicate (length fStatSpecs) 0
        valVec <- liftIO $ VUM.replicate (length fStatSpecs) 0.0
        return (Nothing, Nothing, 0, normVec, valVec)
    extract :: (MonadIO m) => (Maybe GenomPos, Maybe GenomPos, Int, VUM.IOVector Int, VUM.IOVector Double) -> m BlockData
    extract (maybeStartPos, maybeEndPos, count, normVec, valVec) = do
        let Just startPos = maybeStartPos
            Just endPos = maybeEndPos
        normVecF <- liftIO $ VU.freeze normVec
        valVecF <- liftIO $ VU.freeze valVec
        let statVals = zipWith (\v n -> v / fromIntegral n) (VU.toList valVecF) (VU.toList normVecF)
        return $ BlockData startPos endPos count statVals

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
runFstats opts = do
    -- load packages --
    allPackages <- readPoseidonPackageCollection pacReadOpts (_foBaseDirs opts)
    statSpecsFromFile <- case _foStatSpecsFile opts of
        Nothing -> return []
        Just f  -> readStatSpecsFromFile f
    let statSpecs = statSpecsFromFile ++ _foStatSpecsDirect opts
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
            let jointIndInfo = getJointIndividualInfo relevantPackages
            let statsFold = buildStatSpecsFold jointIndInfo statSpecs
            blocks <- runSafeT $ do
                (_, eigenstratProd) <- getJointGenotypeData False False relevantPackages Nothing
                let eigenstratProdFiltered = eigenstratProd >-> P.filter chromFilter >-> capNrSnps (_foMaxSnps opts)
                    eigenstratProdInChunks = case _foJackknifeMode opts of
                        JackknifePerChromosome  -> chunkEigenstratByChromosome eigenstratProdFiltered
                        JackknifePerN chunkSize -> chunkEigenstratByNrSnps chunkSize eigenstratProdFiltered
                let summaryStatsProd = impurely foldsM statsFold eigenstratProdInChunks
                purely P.fold list (summaryStatsProd >-> P.tee (P.map showBlockLogOutput >-> P.toHandle stderr))
            let jackknifeWeights = map blockSiteCount blocks
            let jackknifeEstimates = do
                    (i, _) <- zip [0..] statSpecs
                    return $ computeJackknife jackknifeWeights (map ((!!i) . blockStatVal) blocks)
            let (colSpecs, tableH, tableB) =
                    if _foFullTable opts then
                        let c = replicate 5 (column expand def def def)
                            tH = ["Block", "Statistic", "Estimate", "StdErr", "Z score"]
                            tB = concat $ do
                                (i, fstat, result) <- zip3 [0..] statSpecs jackknifeEstimates
                                let perStat = do
                                        (j, block) <- zip [0..] blocks
                                        return [show j, show fstat, show (blockStatVal block !! i), "n/a", "n/a"]
                                return $ perStat ++ [["Full", show fstat, show (fst result), show (snd result), show (uncurry (/) result)]]
                        in  (c, tH, tB)
                    else 
                        let c = replicate 4 (column expand def def def)
                            tH = ["Statistic", "Estimate", "StdErr", "Z score"]
                            tB = do
                                (fstat, result) <- zip statSpecs jackknifeEstimates
                                return [show fstat, show (fst result), show (snd result), show (uncurry (/) result)]
                        in  (c, tH, tB)
            if   _foRawOutput opts
            then do
                putStrLn $ intercalate "\t" tableH
                forM_ tableB $ \row -> putStrLn (intercalate "\t" row)
            else putStrLn $ tableString colSpecs asciiRoundS (titlesH tableH) [rowsG tableB]
  where
    chromFilter (EigenstratSnpEntry chrom _ _ _ _ _, _) = chrom `notElem` _foExcludeChroms opts
    capNrSnps Nothing  = cat
    capNrSnps (Just n) = P.take n
    chunkEigenstratByChromosome = view (groupsBy sameChrom)
    sameChrom (EigenstratSnpEntry chrom1 _ _ _ _ _, _) (EigenstratSnpEntry chrom2 _ _ _ _ _, _) =
        chrom1 == chrom2
    chunkEigenstratByNrSnps chunkSize = view (chunksOf chunkSize)
    showBlockLogOutput block = "computing chunk range " ++ show (blockStartPos block) ++ " - " ++
        show (blockEndPos block) ++ ", size " ++ (show . blockSiteCount) block ++ " SNPs"

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

