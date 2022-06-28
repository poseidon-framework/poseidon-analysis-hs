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
                                              computeJackknifeOriginal)

import           Control.Applicative         ((<|>))
import           Control.Exception           (throwIO)
import           Control.Foldl               (FoldM (..), impurely, list,
                                              purely)
import           Control.Monad               (forM_, when, forM)
import           Control.Monad.IO.Class      (MonadIO, liftIO)
import           Data.Char                   (isSpace)
import           Data.IORef                  (IORef, newIORef, writeIORef,
                                              modifyIORef', readIORef)
import           Data.List                   (intercalate, nub)
import qualified Data.Vector                 as V
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
import           Text.Read                   (readMaybe)

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

-- | A datatype to represent different types of F-Statistics
data FStatType = F4 | F3 | F2 | PWM | Het | FST | F3vanilla | F2vanilla | FSTvanilla
    deriving (Show, Read, Eq)

-- | A datatype to represent F-Statistics to be computed from genotype data.
data FStatSpec = FStatSpec FStatType [PoseidonEntity] (Maybe AscertainmentSpec) deriving (Eq, Show)

-- | An internal datatype to represent Summary statistics with indices of individuals given as integers, and ascertainment information
data FStat = FStat {
    _fType :: FStatType,
    _slotI :: [[Int]],
    _ascOutgroupI :: [Int], -- an empty list here means that minor allele frequencies are to be taken as ascertaining frequency
    _ascRefI :: [Int], -- empty list means no ascertainment
    _ascLower :: Double,
    _ascUpper :: Double
} deriving (Eq)

data AscertainmentSpec = AscertainmentSpec {
    ascOutgroupSpec  :: Maybe PoseidonEntity, -- Here, a Nothing denotes that minor allele frequency in ascRefgroup should be used, 
                                         --otherwise use derived allele frequencies in the ref-group, where the ancestral allele is given by the outgroup
    ascRefgroupSpec  :: PoseidonEntity,
    ascLowerFreqSpec :: Double,
    ascUpperFreqSpec :: Double
} deriving (Show, Eq)

data BlockData = BlockData
    { blockStartPos  :: GenomPos
    , blockEndPos    :: GenomPos
    , blockSiteCount :: Int
    , blockStatVal   :: [[Double]] -- multiple per-block-accumulators per statistics, can be used however the specific statistic needs to.
    -- For example, most stats will use two numbers, one for the accumulating statistic, and one for the normalisation. Some use more, like F3, which 
    -- accumulates not only (c-a)(c-b), but also the heterozygosity for the denominator, both of which also need normalisation, so F3 has four accumulators.
    }
    deriving (Show)

data BlockAccumulator = BlockAccumulator {
    accMaybeStartPos :: IORef (Maybe GenomPos),
    accMaybeEndPos :: IORef (Maybe GenomPos),
    accCount :: IORef Int,
    accValues :: V.Vector (VUM.IOVector Double) -- this is the key value accumulator: for each statistic there is a list of accumulators if needed, see above.
        -- the outer vector is boxed and immutable, the inner vector is unboxed and mutable, so that values can be updated as we loop through the data.
}

summaryPrintFstats :: FStatSpec -> String
summaryPrintFstats (FStatSpec fType slots maybeAsc) =
    let ascString = case maybeAsc of
            Nothing -> ""
            Just _ -> "_ascertained"
    in  show fType ++ "("  ++ intercalate "," (map show slots) ++ ")" ++ ascString

-- | A parser to parse Summary Statistic specifications from the simple text file input. Every line is one statistics. No ascertainment can be given with this interface.
fStatSpecParser :: P.Parser FStatSpec
fStatSpecParser = do
    typeStr <- P.many1 (P.satisfy (\c -> c /= '('))
    fStatType <- case readMaybe typeStr of
        Just t -> return t
        Nothing -> fail $ "Cannot parse Statistic type " ++ typeStr ++
            ". Must be one of F4, F3, F2, FST, PWM, Het, F3vanilla, FSTvanilla"
    slots <- P.between (P.char '(') (P.char ')') parseEntities
    when (not (checkFstatSlotLength fStatType slots)) $
        fail $ "Not the right number of arguments to Statistic " ++ show fStatType
    return $ FStatSpec fStatType slots Nothing -- no ascertainment can be specified using this input interface.
  where
    parseEntities = P.sepBy1 customEntitySpecParser (P.char ',' <* P.spaces)
    customEntitySpecParser = parsePac <|> parseGroup <|> parseInd
      where
        parsePac   = Pac   <$> P.between (P.char '*') (P.char '*') parseName
        parseGroup = Group <$> parseName
        parseInd   = Ind   <$> P.between (P.char '<') (P.char '>') parseName
        parseName  = P.many1 (P.satisfy (\c -> not (isSpace c || c `elem` charList)))
        charList = ",<>*()" :: [Char] -- we use a custom parser here, because we cannot tolerate bracket openings, which is something that is not constrained in 
            -- Poseidon.EntityList.
    checkFstatSlotLength fStatType slots =
        let n = length slots
        in  case fStatType of
                F4         -> n == 4
                F3         -> n == 3
                F2         -> n == 2
                FST        -> n == 2
                PWM        -> n == 2
                Het        -> n == 1
                F3vanilla  -> n == 2
                FSTvanilla -> n == 2
                F2vanilla  -> n == 2

-- This functioin builds the central Fold that is run over each block of sites of the input data. The return is a tuple of the internal FStats datatypes and the fold.
buildStatSpecsFold :: (MonadIO m) => [IndividualInfo] -> [FStatSpec] -> m ([FStat], FoldM m (EigenstratSnpEntry, GenoLine) BlockData)
buildStatSpecsFold indInfo fStatSpecs = do
    let getI = (flip conformingEntityIndices) indInfo . return
        fStats = do
            FStatSpec t slots maybeAsc <- fStatSpecs
            let (ascOG, ascRef, lo, up) = case maybeAsc of
                    Nothing -> ([], [], 0.0, 1.0)
                    Just (AscertainmentSpec ogE refE lo' up') -> (maybe [] getI ogE, getI refE, lo', up')
            return $ FStat t (map getI slots) ascOG ascRef lo up
    blockAccum <- initialize fStats
    return $ (fStats, FoldM (step fStats blockAccum) (return ()) (extract blockAccum))
  where
    step :: (MonadIO m) => [FStat] -> BlockAccumulator -> () -> (EigenstratSnpEntry, GenoLine) -> m ()
    step fStats blockAccum _ (EigenstratSnpEntry c p _ _ _ _, genoLine) = do
        -- this function is called for every SNP.
        startPos <- liftIO $ readIORef (accMaybeStartPos blockAccum)
        case startPos of
            Nothing       -> liftIO $ writeIORef (accMaybeStartPos blockAccum) (Just (c, p))
            Just _ -> return ()
        liftIO $ writeIORef (accMaybeEndPos blockAccum) (Just (c, p))
        forM_ (zip [0..] fStats) $ \(i, fstat) -> do
            -- loop over all statistics
            let maybeAccValues = computeFStatAccumulators fstat genoLine -- compute the accumulating values for that site.
            forM_ (zip [0..] maybeAccValues) $ \(j, maybeAccVal) -> do
                case maybeAccVal of
                    Nothing -> return ()
                    Just x -> liftIO $ VUM.modify (accValues blockAccum V.! i) (+x) j -- add the value to the respective accumulator.
        -- counts the number of SNPs in a block, and ignores missing data. Missing data is considered within each accumulator.
        liftIO $ modifyIORef' (accCount blockAccum) (+1)
        return ()
    initialize :: (MonadIO m) => [FStat] -> m BlockAccumulator
    initialize fStats = do
        listOfInnerVectors <- forM fStats $ \(FStat fType _ _ _ _ _) -> do
            case fType of
                F3 -> liftIO $ VUM.replicate 4 0.0 --only F3 has four accumulators: one numerator, one denominator, and one normaliser for each of the two.
                _  -> liftIO $ VUM.replicate 2 0.0 -- all other statistics have just one value and one normaliser.      
        liftIO $ BlockAccumulator <$> newIORef Nothing <*> newIORef Nothing <*> newIORef 0 <*> pure (V.fromList listOfInnerVectors)
    extract :: (MonadIO m) => BlockAccumulator -> () -> m BlockData
    extract (BlockAccumulator maybeStartPosRef maybeEndPosRef countRef valVec) _ = do
        maybeStartPos <- liftIO $ readIORef maybeStartPosRef
        maybeEndPos <- liftIO $ readIORef maybeEndPosRef
        count <- liftIO $ readIORef countRef
        statVals <- liftIO $ mapM (fmap VU.toList . VU.freeze) (V.toList valVec) -- this unfolds the vector of mutable vectors into a list of lists.
        case (maybeStartPos, maybeEndPos) of
            (Just startPos, Just endPos) -> return $ BlockData startPos endPos count statVals
            _ -> error "should never happen"

-- TODO: Currently ignores ascertainment!
computeFStatAccumulators :: FStat -> GenoLine -> [Maybe Double] -- returns a number of accumulated variables, in most cases a value and a normalising count,
-- but in case of F3, for example, also a second accumulator capturing the heterozygosity
computeFStatAccumulators (FStat fType indices ascOgI ascRefI ascLo ascHi) gL =
    let caf = computeAlleleFreq gL -- this returns Nothing if missing data
        cac i = case computeAlleleCount gL i of -- this also returns Nothing if missing data.
            (_, 0) -> Nothing
            x      -> Just x
        ascCond = 
            if null ascRefI then
                return True
            else do -- Maybe Monad - if any of the required allele frequencies are missing, this returns a Nothing
                let ascRefX = caf ascRefI
                ascFreq <- if null ascOgI then do
                        x <- ascRefX
                        if x > 0.5 then return (1.0 - x) else return x -- use minor allele frequency if no outgroup is given
                    else do -- Maybe Monad
                        ogX <- caf ascOgI
                        x <- ascRefX
                        if (ogX < 0.5) then return x else return (1.0 - x)
                return $ ascFreq >= ascLo && ascFreq <= ascHi
        applyAsc x = do -- Maybe Monad
            cond <- ascCond
            if cond then return x else return 0.0
    in  case (fType, indices) of
            (F4,         [aI, bI, cI, dI]) -> retWithNormAcc $ (computeF4         <$> caf aI <*> caf bI <*> caf cI <*> caf dI) >>= applyAsc
            (F3vanilla,  [aI, bI, cI])     -> retWithNormAcc $ (computeF3vanilla  <$> caf aI <*> caf bI <*> caf cI)            >>= applyAsc
            (F2vanilla,  [aI, bI])         -> retWithNormAcc $ (computeF2vanilla  <$> caf aI <*> caf bI)                       >>= applyAsc
            (PWM,        [aI, bI])         -> retWithNormAcc $ (computePWM        <$> caf aI <*> caf bI)                       >>= applyAsc
            (Het,        [aI])             -> retWithNormAcc $ (computeHet        <$> cac aI)                                  >>= applyAsc
            (F2,         [aI, bI])         -> retWithNormAcc $ (computeF2         <$> cac aI <*> cac bI)                       >>= applyAsc
            (FSTvanilla, [aI, bI])         -> retWithNormAcc $ (computeFSTvanilla    (caf aI)   (caf bI))                      >>= applyAsc
            (FST,        [aI, bI])         -> retWithNormAcc $ (computeFST           (cac aI)   (cac bI))                      >>= applyAsc
            (F3,         [aI, bI, cI])     ->
                retWithNormAcc ((computeF3noNorm   <$> caf aI <*> caf bI <*> cac cI) >>= applyAsc) ++
                retWithNormAcc ((computeHet <$> cac cI) >>= applyAsc)
            _ -> error "should never happen"
  where
    retWithNormAcc (Just x) = [Just x, Just 1.0]
    retWithNormAcc Nothing = [Nothing, Nothing]
    -- these formulas are mostly taken from Patterson et al. 2012 Appendix A (page 25 in the PDF)
    computeF4         a b c d = (a - b) * (c - d)
    computeF3vanilla  a b c   = (c - a) * (c - b)
    computeF2vanilla  a b     = (a - b) * (a - b)
    computePWM        a b     = a * (1.0 - b) + (1.0 - a) * b
    computeHet (na, sa) = 2.0 * fromIntegral (na * (sa - na)) / fromIntegral (sa * (sa - 1))
    computeF3noNorm a b (nc, sc) =
        let c = computeFreq nc sc
            corrFac = 0.5 * computeHet (nc, sc) / fromIntegral sc
        in  computeF3vanilla a b c - corrFac
    computeF2 (na, sa) (nb, sb) =
        let a = computeFreq na sa
            b = computeFreq nb sb
            corrFac = 0.5 * computeHet (na, sa) / fromIntegral sa + 0.5 * computeHet (nb, sb) / fromIntegral sb
        in  computeF2vanilla a b - corrFac
    computeFSTvanilla maybeA maybeB =
        case (maybeA, maybeB) of
            (Just a, Just b) ->
                let num = (a - b)^(2 :: Int)
                    denom = a * (1 - b) + b * (1 - a)
                in  if denom > 0 then Just (num / denom) else Nothing
            _ -> Nothing
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
            hPutStrLn stderr $ "Computing stats:"
            mapM_ (hPutStrLn stderr . summaryPrintFstats) statSpecs
            let jointIndInfo = getJointIndividualInfo relevantPackages
            (fStats, blocks) <- runSafeT $ do
                (fStats, statsFold) <- buildStatSpecsFold jointIndInfo statSpecs
                (_, eigenstratProd) <- getJointGenotypeData False False relevantPackages Nothing
                let eigenstratProdFiltered = eigenstratProd >-> P.filter chromFilter >-> capNrSnps (_foMaxSnps opts)
                    eigenstratProdInChunks = case _foJackknifeMode opts of
                        JackknifePerChromosome  -> chunkEigenstratByChromosome eigenstratProdFiltered
                        JackknifePerN chunkSize -> chunkEigenstratByNrSnps chunkSize eigenstratProdFiltered
                let summaryStatsProd = impurely foldsM statsFold eigenstratProdInChunks
                blocks <- purely P.fold list (summaryStatsProd >-> P.tee (P.map showBlockLogOutput >-> P.toHandle stderr))
                return (fStats, blocks)
            let jackknifeEstimates = processBlocks fStats blocks
            let colSpecs = replicate 10 (column expand def def def)
                tableH = ["Statistic", "a", "b", "c", "d", "Asc (Og, Ref)", "Asc (Lo, Up)", "Estimate", "StdErr", "Z score"]
                tableB = do
                    (fstat, result) <- zip statSpecs jackknifeEstimates
                    let FStatSpec fType slots maybeAsc = fstat
                        abcdStr = take 4 (map show slots ++ repeat "")
                        (asc1, asc2) = case maybeAsc of
                            Just (AscertainmentSpec (Just og) ref lo up) -> (show (og, ref), show (lo, up))
                            Just (AscertainmentSpec Nothing ref lo up) -> (show ("n/a", ref), show (lo, up))
                            _ -> ("n/a", "n/a")
                    return $ [show fType] ++ abcdStr ++ [asc1, asc2] ++ [printf "%.4g" (fst result), printf "%.4g" (snd result), show (uncurry (/) result)]
            putStrLn $ tableString colSpecs asciiRoundS (titlesH tableH) [rowsG tableB]
            case _foTableOut opts of
                Nothing -> return ()
                Just outFn -> withFile outFn WriteMode $ \h -> mapM_ (hPutStrLn h . intercalate "\t") (tableH : tableB)
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
collectStatSpecGroups statSpecs = nub $ concat [slots | FStatSpec _ slots _ <- statSpecs]

processBlocks :: [FStat] -> [BlockData] -> [(Double, Double)]
processBlocks stats blocks = do
    let block_weights = map (fromIntegral . blockSiteCount) blocks
    (i, stat) <- zip [0..] stats
    case _fType stat of
        F3 ->
            let numerator_values = map ((!!0) . (!!i) . blockStatVal) blocks
                numerator_norm = map ((!!1) . (!!i) . blockStatVal) blocks
                denominator_values = map ((!!2) . (!!i) . blockStatVal) blocks
                denominator_norm = map ((!!3) . (!!i) . blockStatVal) blocks
                num_full = sum numerator_values / sum numerator_norm
                denom_full = sum denominator_values / sum denominator_norm
                full_estimate = num_full / denom_full
                partial_estimates = do
                    j <- [0..(length block_weights - 1)]
                    let num = sum [v | (k, v) <- zip [0..] numerator_values, k /= j]
                        num_norm = sum [v | (k, v) <- zip [0..] numerator_norm, k /= j]
                        denom = sum [v | (k, v) <- zip [0..] denominator_values, k /= j]
                        denom_norm = sum [v | (k, v) <- zip [0..] denominator_norm, k /= j]
                    return $ (num / num_norm) / (denom / denom_norm)
            in  return $ computeJackknifeOriginal full_estimate block_weights partial_estimates
        _ ->
            let values = map ((!!0) . (!!i) . blockStatVal) blocks
                norm = map ((!!1) . (!!i) . blockStatVal) blocks
                full_estimate = sum values / sum norm
                partial_estimates = do
                    j <- [0..(length block_weights - 1)]
                    let val' = sum [v | (k, v) <- zip [0..] values, k /= j]
                        norm' = sum [v | (k, v) <- zip [0..] norm, k /= j]
                    return $ val' / norm'
            in  return $ computeJackknifeOriginal full_estimate block_weights partial_estimates


