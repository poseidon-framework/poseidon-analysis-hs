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
                                              computeAlleleCount,
                                              computeAlleleFreq,
                                              computeJackknifeOriginal,
                                              popSpecsNparser)

import           Control.Applicative         ((<|>))
import           Control.Exception           (throwIO)
import           Control.Foldl               (FoldM (..), impurely, list,
                                              purely)
import           Control.Monad               (forM_)
import           Control.Monad.IO.Class      (MonadIO, liftIO)
import           Data.List                   (elemIndex, intercalate, nub)
import qualified Data.Vector.Unboxed         as VU
import qualified Data.Vector.Unboxed.Mutable as VUM
-- import           Debug.Trace                 (trace)
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
import           System.IO                   (IOMode (..), hPutStrLn, stderr,
                                              withFile)
import           Text.Layout.Table           (asciiRoundS, column, def, expand,
                                              rowsG, tableString, titlesH)
import qualified Text.Parsec                 as P
import qualified Text.Parsec.String          as P
import           Text.Printf                 (printf)

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
    , _foMaxSnps         :: Maybe Int
    , _foTableOut        :: Maybe FilePath
    }

data FStatType = F4 | F3 | F2 | PWM | Het | FST | F3vanilla | F2vanilla | FSTvanilla
    deriving (Show, Read)

-- | A datatype to represent Summary Statistics to be computed from genotype data.
data FStatSpec = FStatSpec FStatType [PoseidonEntity] deriving (Eq)

instance Show FStatSpec where
    show (FStatSpec type slots) = show type ++ "("  ++ intercalate "," (map show slots) ++ ")"

-- | An internal datatype to represent Summary statistics with indices of individuals given as integers
data FStat = FStat FStatType [[Int]] deriving (Eq)

data BlockData = BlockData
    { blockStartPos  :: GenomPos
    , blockEndPos    :: GenomPos
    , blockSiteCount :: Int
    , blockStatVal   :: [Double]
    }
    deriving (Show)

-- | A parser to parse Summary Statistic specifications.
fStatSpecParser :: P.Parser FStatSpec
fStatSpecParser = do
        typeStr <- P.satisfy (\c -> c /= '(')
        type <- case readMaybe typeStr of
            Just t -> return t
            Nothing -> fail $ "Cannot parse Statistic type " ++ typeStr ++
                ". Must be one of F4, F3, F2, FST, PWM, Het, F3vanilla, FSTvanilla"
        [a, b, c, d] <- P.between (P.char '(') (P.char ')') (popSpecsNparser 4)
        return $ F4Spec a b c d
    f3SpecParser = do
        _ <- P.string "F3"
        [a, b, c] <- P.between (P.char '(') (P.char ')') (popSpecsNparser 3)
        return $ F3Spec a b c
    f2SpecParser = do
        _ <- P.string "F2"
        [a, b] <- P.between (P.char '(') (P.char ')') (popSpecsNparser 2)
        return $ F2Spec a b
    pwmSpecParser = do
        _ <- P.string "PWM"
        [a, b] <- P.between (P.char '(') (P.char ')') (popSpecsNparser 2)
        return $ PWMspec a b
    hetSpecParser = do
        _ <- P.string "Het"
        [a] <- P.between (P.char '(') (P.char ')') (popSpecsNparser 1)
        return $ HetSpec a
    fstSpecParser = do
        _ <- P.string "FST"
        [a, b] <- P.between (P.char '(') (P.char ')') (popSpecsNparser 2)
        return $ FSTspec a b
    f3VanillaSpecParser = do
        _ <- P.string "F3vanilla"
        [a, b, c] <- P.between (P.char '(') (P.char ')') (popSpecsNparser 3)
        return $ F3vanillaSpec a b c
    f2VanillaSpecParser = do
        _ <- P.string "F2vanilla"
        [a, b] <- P.between (P.char '(') (P.char ')') (popSpecsNparser 2)
        return $ F2vanillaSpec a b
    fstVanillaSpecParser = do
        _ <- P.string "FSTvanilla"
        [a, b] <- P.between (P.char '(') (P.char ')') (popSpecsNparser 2)
        return $ FSTvanillaSpec a b

buildStatSpecsFold :: (MonadIO m) => [IndividualInfo] -> [FStatSpec] -> ([FStat], FoldM m (EigenstratSnpEntry, GenoLine) BlockData)
buildStatSpecsFold indInfo fStatSpecs =
    let getI = conformingEntityIndices
        i = indInfo
        fStatsBasic = do
            fStatSpec <- fStatSpecs
            case fStatSpec of
                F4Spec         a b c d -> return $ F4         (getI [a] i) (getI [b] i) (getI [c] i) (getI [d] i)
                F3Spec         a b c   -> return $ F3         (getI [a] i) (getI [b] i) (getI [c] i)
                F2Spec         a b     -> return $ F2         (getI [a] i) (getI [b] i)
                PWMspec        a b     -> return $ PWM        (getI [a] i) (getI [b] i)
                HetSpec        a       -> return $ Het        (getI [a] i)
                FSTspec        a b     -> return $ FST        (getI [a] i) (getI [b] i)
                F3vanillaSpec  a b c   -> return $ F3vanilla  (getI [a] i) (getI [b] i) (getI [c] i)
                F2vanillaSpec  a b     -> return $ F2vanilla  (getI [a] i) (getI [b] i)
                FSTvanillaSpec a b     -> return $ FST        (getI [a] i) (getI [b] i)
        additionalHetStats = [Het cI | F3 _ _ cI <- fStatsBasic]
        fStatsAll = fStatsBasic ++ additionalHetStats
        n = length fStatsAll
    in  (fStatsAll, FoldM (step fStatsAll) (initialize n) extract)
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
    initialize :: (MonadIO m) => Int -> m (Maybe GenomPos, Maybe GenomPos, Int, VUM.IOVector Int, VUM.IOVector Double)
    initialize n = do
        normVec <- liftIO $ VUM.replicate n 0
        valVec <- liftIO $ VUM.replicate n 0.0
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
computeFStat fStat gL =
    let caf = computeAlleleFreq
        cac gL' i = case computeAlleleCount gL' i of
            (_, 0) -> Nothing
            x      -> Just x
    in  case fStat of
            F4         aI bI cI dI -> computeF4         <$> caf gL aI <*> caf gL bI <*> caf gL cI <*> caf gL dI
            F3vanilla  aI bI cI    -> computeF3vanilla  <$> caf gL aI <*> caf gL bI <*> caf gL cI
            F2vanilla  aI bI       -> computeF2vanilla  <$> caf gL aI <*> caf gL bI
            FSTvanilla aI bI       -> computeFSTvanilla <$> caf gL aI <*> caf gL bI
            PWM        aI bI       -> computePWM        <$> caf gL aI <*> caf gL bI
            F3         aI bI cI    -> computeF3         <$> caf gL aI <*> caf gL bI <*> cac gL cI
            Het        aI          -> computeHet        <$> cac gL aI
            F2         aI bI       -> computeF2         <$> cac gL aI <*> cac gL bI
            FST        aI bI       -> computeFST        (cac gL aI) (cac gL bI)
  where
    -- these formulas are mostly taken from Patterson et al. 2012 Appendix A (page 25 in the PDF)
    computeF4         a b c d = (a - b) * (c - d)
    computeF3vanilla  a b c   = (c - a) * (c - b)
    computeF2vanilla  a b     = (a - b) * (a - b)
    computeFSTvanilla a b     = (a - b)^(2 :: Int) / (a * (1 - b) + b * (1 - a))
    computePWM        a b     = a * (1.0 - b) + (1.0 - a) * b
    computeHet (na, sa) = 2.0 * fromIntegral (na * (sa - na)) / fromIntegral (sa * (sa - 1))
    computeF3 a b (nc, sc) =
        let c = computeFreq nc sc
            corrFac = 0.5 * computeHet (nc, sc) / fromIntegral sc
        in  computeF3vanilla a b c - corrFac
    computeF2 (na, sa) (nb, sb) =
        let a = computeFreq na sa
            b = computeFreq nb sb
            corrFac = 0.5 * computeHet (na, sa) / fromIntegral sa + 0.5 * computeHet (nb, sb) / fromIntegral sb
        in  computeF2vanilla a b - corrFac
    computeFST maybeNaSa maybeNbSb = case (maybeNaSa, maybeNbSb) of
        (Just (na, sa), Just (nb, sb)) ->
            let num = computeF2 (na, sa) (nb, sb)
                denom = computeF2 (na, sa) (nb, sb) + 0.5 * computeHet (na, sa) + 0.5 * computeHet (nb, sb)
            in  if denom > 0 then Just (num / denom) else Nothing
        _ -> Nothing
    computeFreq na sa = fromIntegral na / fromIntegral sa

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
            let (fStats, statsFold) = buildStatSpecsFold jointIndInfo statSpecs
            blocks <- runSafeT $ do
                (_, eigenstratProd) <- getJointGenotypeData False False relevantPackages Nothing
                let eigenstratProdFiltered = eigenstratProd >-> P.filter chromFilter >-> capNrSnps (_foMaxSnps opts)
                    eigenstratProdInChunks = case _foJackknifeMode opts of
                        JackknifePerChromosome  -> chunkEigenstratByChromosome eigenstratProdFiltered
                        JackknifePerN chunkSize -> chunkEigenstratByNrSnps chunkSize eigenstratProdFiltered
                let summaryStatsProd = impurely foldsM statsFold eigenstratProdInChunks
                purely P.fold list (summaryStatsProd >-> P.tee (P.map showBlockLogOutput >-> P.toHandle stderr))
            let jackknifeEstimates = processBlocks statSpecs fStats blocks
            let colSpecs = replicate 4 (column expand def def def)
                tableH = ["Statistic", "Estimate", "StdErr", "Z score"]
                tableB = do
                    (fstat, result) <- zip statSpecs jackknifeEstimates
                    return [show fstat, printf "%.4g" (fst result), printf "%.4g" (snd result), show (uncurry (/) result)]
            putStrLn $ tableString colSpecs asciiRoundS (titlesH tableH) [rowsG tableB]
            case _foTableOut opts of
                Nothing -> return ()
                Just outFn -> do
                    withFile outFn WriteMode $ \h -> do
                        hPutStrLn h . intercalate "\t" $ ["Statistic", "a", "b", "c", "d", "Estimate", "StdErr", "Z score"]
                        forM_ (zip statSpecs jackknifeEstimates) $ \(fstat, result) -> do
                            let abcd = collectStatSpecGroups [fstat]
                            let abcdStr = take 4 (map show abcd ++ repeat "")
                            hPutStrLn h . intercalate "\t" $ [show fstat] ++ abcdStr ++ [show (fst result), show (snd result), show (uncurry (/) result)]
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
        F4Spec         a b c d -> return [a, b, c, d]
        F3Spec         a b c   -> return [a, b, c]
        F2Spec         a b     -> return [a, b]
        PWMspec        a b     -> return [a, b]
        HetSpec        a       -> return [a]
        FSTspec        a b     -> return [a, b]
        F3vanillaSpec  a b c   -> return [a, b, c]
        F2vanillaSpec  a b     -> return [a, b]
        FSTvanillaSpec a b     -> return [a, b]

processBlocks :: [FStatSpec] -> [FStat] -> [BlockData] -> [(Double, Double)]
processBlocks statSpecs stats blocks = do
    (i, statSpec, stat) <- zip3 [0..] statSpecs stats
    case statSpec of
        F3Spec {} ->
            let F3 _ _ cI = stat
                relatedHetIndex = case elemIndex (Het cI) stats of
                    Nothing -> error "should never happen, cannot find related het-statistics for F3 stat"
                    Just j -> j
                numerator_values = map ((!!i) . blockStatVal) blocks
                denominator_values = map ((!!relatedHetIndex) . blockStatVal) blocks
                block_weights = map (fromIntegral . blockSiteCount) blocks
                num_full = sum [m * v / sum block_weights | (m, v) <- zip block_weights numerator_values]
                denom_full = sum [m * v / sum block_weights | (m, v) <- zip block_weights denominator_values]
                full_estimate = num_full / denom_full
                partial_estimates = do
                    j <- [0..(length block_weights - 1)]
                    let weight_norm = sum [m | (k, m) <- zip [0..] block_weights, k /= j]
                        num = sum [m * v / weight_norm | (k, m, v) <- zip3 [0..] block_weights numerator_values, k /= j]
                        denom = sum [m * v / weight_norm | (k, m, v) <- zip3 [0..] block_weights denominator_values, k /= j]
                    return $ num / denom
            in  return $ computeJackknifeOriginal full_estimate block_weights partial_estimates
        _ ->
            let values = map ((!!i) . blockStatVal) blocks
                block_weights = map (fromIntegral . blockSiteCount) blocks
                full_estimate = sum [m * v / sum block_weights | (m, v) <- zip block_weights values]
                partial_estimates = do
                    j <- [0..(length block_weights - 1)]
                    let weight_norm = sum [m | (k, m) <- zip [0..] block_weights, k /= j]
                    return $ sum [m * v / weight_norm | (k, m, v) <- zip3 [0..] block_weights values, k /= j]
            in  return $ computeJackknifeOriginal full_estimate block_weights partial_estimates


