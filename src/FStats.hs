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
import           Control.Monad               (forM_, when)
import           Control.Monad.IO.Class      (MonadIO, liftIO)
import           Data.Char                   (isSpace)
import           Data.List                   (elemIndex, intercalate, nub)
import qualified Data.Vector                 as V
import qualified Data.Vector.Unboxed         as VU
import qualified Data.Vector.Unboxed.Mutable as VUM
import qualified Data.Vector.Mutable         as VM
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
    _fType :: FStatType
    _slotI :: [[Int]]
    _ascOutgroupI :: [Int] -- an empty list here means that minor allele frequencies are to be taken as ascertaining frequency
    _ascRefI :: [Int] -- empty list means no ascertainment
    _ascLower :: Double
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
    accMaybeStartPos :: Maybe GenomPos
    accMaybeEndPos :: Maybe GenomPos
    accCount :: Int
    accValues :: V.IOVector (VUM.IOVector Double) -- this is the key value accumulator: for each statistic there is a list of accumulators if needed, see above.
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
buildStatSpecsFold :: (MonadIO m) => [IndividualInfo] -> [FStatSpec] -> ([FStat], FoldM m (EigenstratSnpEntry, GenoLine) BlockData)
buildStatSpecsFold indInfo fStatSpecs =
    let getI = (flip conformingEntityIndices) indInfo . return
        fStats = do
            FStatSpec t slots maybeAsc <- fStatSpecs
            let (ascOG, ascRef, lo, up) = case maybeAsc of
                    Nothing -> ([], [], 0.0, 1.0)
                    Just (AscertainmentSpec ogE refE lo' up') -> (maybe [] getI ogE, getI refE, lo', up')
            return $ FStat t (map getI slots) ascOG ascRef lo up
        n = length fStats
    in  (fStats, FoldM (step fStats) (initialize fStats) (extract fStats))
  where
    step :: (MonadIO m) => [FStat] -> BlockAccumulator -> (EigenstratSnpEntry, GenoLine) -> m BlockAccumulator
    step fStats blockAccum (EigenstratSnpEntry c p _ _ _ _, genoLine) = do
        -- this function is called for every SNP.
        let newStartPos = case accMaybeStartPos blockAccum of
                Nothing       -> Just (c, p)
                Just (c', p') -> Just (c', p')
            newEndPos = Just (c, p)
        forM_ (zip [0..] fStats) $ \(i, fstat) -> do
            -- loop over all statistics
            let maybeAccValues = computeFStatAccumulators fstat genoLine -- compute the accumulating values for that site.
            forM_ (zip [0..] maybeAccValues) $ \(j, maybeAccVal) -> do
                case maybeAccVal of
                    Nothing -> return ()
                    Just x -> liftIO $ VUM.modify (accValues blockAccum V.! i) (+x) j -- add the value to the respective accumulator.
                Nothing -> return ()
        return $ blockAccum {
            accMaybeStartPos = newStartPos,
            accMaybeEndPos   = newEndPos,
            accCount         = accCount blockAccum + 1} -- also increment the overall block size by the number of SNPs in that blockl. Note that this simply 
            -- counts the number of SNPs in a block, and ignores missing data. Missing data is considered within each accumulator.
    initialize :: (MonadIO m) => [FStat] -> m BlockAccumulator
    initialize fStats = do
        listOfInnerVectors <- forM fStats $ \(FStat fType _ _ _ _ _) -> do
            case fType of
                F3 -> VUM.replicate 4 0.0 --only F3 has four accumulators: one numerator, one denominator, and one normaliser for each of the two.
                _  -> VUM.replicate 2 0.0 -- all other statistics have just one value and one normaliser.      
        return $ BlockAccumulator Nothing Nothing 0 (V.fromList listOfInnerVectors)
    extract :: (MonadIO m) => BlockAccumulator -> m BlockData
    extract (BlockAccumulator maybeStartPos maybeEndPos count valVec) = do
        let Just startPos = maybeStartPos
            Just endPos = maybeEndPos
        statVals <- mapM (fmap VU.toList . VU.freeze) (V.toList valVec) -- this unfolds the vector of mutable vectors into a list of lists.
        return $ BlockData startPos endPos count statVals

-- TODO: Currently ignores ascertainment!
computeFStatAccumulators :: FStat -> GenoLine -> [Maybe Double] -- returns a number of accumulated variables, in most cases a value and a normalising count,
-- but in case of F3, for example, also a second accumulator capturing the heterozygosity
computeFStatAccumulators (FStat fType indices) gL =
    let caf = computeAlleleFreq -- this returns Nothing if missing data
        cac gL' i = case computeAlleleCount gL' i of -- this also returns Nothing if missing data.
            (_, 0) -> Nothing
            x      -> Just x
    in  case (fType, indices) of
            (F4,         [aI, bI, cI, dI]) -> computeF4         (caf gL aI) (caf gL bI) (caf gL cI) (caf gL dI)
            (F3vanilla,  [aI, bI, cI])     -> computeF3vanilla  (caf gL aI) (caf gL bI) (caf gL cI)
            (F2vanilla,  [aI, bI])         -> computeF2vanilla  (caf gL aI) (caf gL bI)
            (FSTvanilla, [aI, bI])         -> computeFSTvanilla (caf gL aI) (caf gL bI)
            (PWM,        [aI, bI])         -> computePWM        (caf gL aI) (caf gL bI)
            (F3,         [aI, bI, cI])     -> computeF3         (caf gL aI) (caf gL bI) (cac gL cI)
            (Het,        [aI])             -> computeHet        (cac gL aI)
            (F2,         [aI, bI])         -> computeF2         (cac gL aI) (cac gL bI)
            (FST,        [aI, bI])         -> computeFST        (cac gL aI) (cac gL bI)
            _ -> error "should never happen"
  where
    -- these formulas are mostly taken from Patterson et al. 2012 Appendix A (page 25 in the PDF)
    -- all of these functions return a list of Maybes, with Nothing denoting missing data for this specific
    -- accumulator.
    computeF4 (Just a) (Just b) (Just c) (Just d) = [Just (a - b) * (c - d), Just 1.0]
    computeF4 _                                   = [Nothing, Nothing]
    
    computeF3vanilla (Just a) (Just b) (Just c) = [Just (c - a) * (c - b), Just 1.0]
    computeF3vanilla _                          = [Nothing, Nothing]
    
    computeF2vanilla (Just a) (Just b) = [Just (a - b) * (a - b), Just 1.0]
    computeF2vanilla _                 = [Nothing, Nothing]
    
    computeFSTvanilla (Just a) (Just b) = [Just (a - b)^(2 :: Int) / (a * (1 - b) + b * (1 - a)), Just 1.0]
    computeFSTvanilla _                 = [Nothing, Nothing]
    
    computePWM (Just a) (Just b) = [Just a * (1.0 - b) + (1.0 - a) * b, Just 1.0]
    computeWPM _                 = [Nothing, Nothing]

    computeHet (Just (na, sa)) = [Just 2.0 * fromIntegral (na * (sa - na)) / fromIntegral (sa * (sa - 1)), Just 1.0]
    computeHet Nothing         = [Nothing, Nothing]
    
    -- this is the only statistic currently returning four accumulators.
    computeF3 maybeA maybeB maybeCcounts =
        let (numAcc, numNorm) = case (maybeA, maybeB) of
                (Just a, Just b) -> 
                    let c = computeFreq nc sc
                        corrFac = 0.5 * computeHet (nc, sc) / fromIntegral sc
                    in  (Just (computeF3vanilla a b c - corrFac), Just 1.0)
                _ -> (Nothing, Nothing)
            (denomAcc, denomNorm) = case maybeCcounts of
                Just (nc, sc) -> (Just (computeHet (nc, sc)), Just 1.0)
                _ -> (Nothing, Nothing)
        in  [numAcc, numNorm, denomAcc, denomNorm]
    
    computeF2 (Just (na, sa)) (Just (nb, sb)) =
        let a = computeFreq na sa
            b = computeFreq nb sb
            corrFac = 0.5 * computeHet (na, sa) / fromIntegral sa + 0.5 * computeHet (nb, sb) / fromIntegral sb
        in  [Just (computeF2vanilla a b - corrFac), Just 1.0]
    computeF2 _ = [Nothing, Nothing]

    computeFST (Just (na, sa)) (Just (nb, sb)) =
        let num = computeF2 (na, sa) (nb, sb)
            denom = computeF2 (na, sa) (nb, sb) + 0.5 * computeHet (na, sa) + 0.5 * computeHet (nb, sb)
        in  if denom > 0 then [Just (num / denom), Just 1.0] else [Nothing, Nothing]
    computeFST _ = [Nothing, Nothing]

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
            let (fStats, statsFold) = buildStatSpecsFold jointIndInfo statSpecs
            blocks <- runSafeT $ do
                (_, eigenstratProd) <- getJointGenotypeData False False relevantPackages Nothing
                let eigenstratProdFiltered = eigenstratProd >-> P.filter chromFilter >-> capNrSnps (_foMaxSnps opts)
                    eigenstratProdInChunks = case _foJackknifeMode opts of
                        JackknifePerChromosome  -> chunkEigenstratByChromosome eigenstratProdFiltered
                        JackknifePerN chunkSize -> chunkEigenstratByNrSnps chunkSize eigenstratProdFiltered
                let summaryStatsProd = impurely foldsM statsFold eigenstratProdInChunks
                purely P.fold list (summaryStatsProd >-> P.tee (P.map showBlockLogOutput >-> P.toHandle stderr))
            let jackknifeEstimates = processBlocks fStats blocks
            let colSpecs = replicate 4 (column expand def def def)
                tableH = ["Statistic", "a", "b", "c", "d", "Asc (Og, Ref)", "Asc (Lo, Up)", "Estimate", "StdErr", "Z score"]
                tableB = do
                    (fstat, result) <- zip statSpecs jackknifeEstimates
                    let FStatSpec fType slots maybeAsc = fstat
                        abcdStr = take 4 (map show slots ++ repeat "")
                        (asc1, asc2) = case maybeAsc of
                            Just (AscertainmentSpec (Just og) ref lo up) -> (show (og, ref), show (lo, up))
                            Just (AscertainmentSpec Nothing ref lo up) -> (show ("n/a", ref), show (lo, up))
                            _ -> ("n/a", "n/a")
                    return $ [show fType] ++ abcdStr ++ [asc1, asc2] ++ [( printf "%.4g" (fst result), printf "%.4g" (snd result), show (uncurry (/) result)]
            putStrLn $ tableString colSpecs asciiRoundS (titlesH tableH) [rowsG tableB]
            case _foTableOut opts of
                Nothing -> return ()
                Just outFn -> withFile outFn WriteMode $ \h -> mapM_ (hPutStrLn h . intercalate "\t") tableB
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
    (i, stat) <- zip [0..] stats
    case _fType stat of
        F3 ->
            let numerator_values = map ((!!0) . (!!i) . blockStatVal) blocks
                numerator_norm = map ((!!1) . (!!i) . blockStatVal) blocks
                denominator_values = map ((!!2) . (!!i) . blockStatVal) blocks
                denominator_norm = map ((!!3) . (!!i) . blockStatVal) blocks
                block_weights = map (fromIntegral . blockSiteCount) blocks
                num_full = sum numerator_values / sum numerator_norm
                denom_full = sum denominator_values / sum denominator_norm
                full_estimate = num_full / denom_full
                partial_estimates = do
                    j <- [0..(length block_weights - 1)]
                    let weight_norm = sum [m | (k, m) <- zip [0..] block_weights, k /= j]
                        num = sum [v | (k, v) <- zip3 [0..] block_weights numerator_values, k /= j]
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


