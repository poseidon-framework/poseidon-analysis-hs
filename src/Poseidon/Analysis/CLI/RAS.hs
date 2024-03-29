{-# LANGUAGE OverloadedStrings #-}

module Poseidon.Analysis.CLI.RAS where

import           Poseidon.Analysis.RASconfig (PopConfig (..))
import           Poseidon.Analysis.Utils     (GenomPos, JackknifeMode (..),
                                              PloidyVec, XerxesException (..),
                                              addGroupDefs, computeAlleleCount,
                                              computeAlleleFreq,
                                              computeJackknifeAdditive,
                                              computeJackknifeOriginal,
                                              filterTransitions, makePloidyVec)

import           Control.Exception           (catch, throwIO)
import           Control.Foldl               (FoldM (..), impurely, list,
                                              purely)
import           Control.Monad               (forM_, unless, when)
import           Control.Monad.IO.Class      (MonadIO, liftIO)
import qualified Data.ByteString             as B
import           Data.List                   (intercalate, nub, (\\))
import qualified Data.Vector                 as V
import qualified Data.Vector.Unboxed         as VU
import qualified Data.Vector.Unboxed.Mutable as VUM
import           Data.Yaml                   (decodeEither')
-- import           Debug.Trace                 (trace)
import           Lens.Family2                (view)

import           Pipes                       (cat, (>->))
import           Pipes.Group                 (chunksOf, foldsM, groupsBy)
import qualified Pipes.Prelude               as P
import           Pipes.Safe                  (runSafeT)
import           Poseidon.EntityTypes        (EntitiesList, IndividualInfo (..),
                                              PoseidonEntity (..),
                                              checkIfAllEntitiesExist,
                                              determineRelevantPackages,
                                              resolveUniqueEntityIndices,
                                              underlyingEntity)
import           Poseidon.Janno              (JannoGenotypePloidy (..))
import           Poseidon.Package            (PackageReadOptions (..),
                                              PoseidonPackage (..),
                                              defaultPackageReadOptions,
                                              getJointGenotypeData,
                                              getJointIndividualInfo,
                                              getJointJanno,
                                              readPoseidonPackageCollection)
import           Poseidon.Utils              (PoseidonException (..),
                                              PoseidonIO, envInputPlinkMode,
                                              envLogAction, logInfo, logWithEnv)
import           SequenceFormats.Bed         (filterThroughBed, readBedFile)
import           SequenceFormats.Eigenstrat  (EigenstratSnpEntry (..),
                                              GenoEntry (..), GenoLine)
import           SequenceFormats.Genomic     (genomicPosition)
import           SequenceFormats.Utils       (Chrom (..))
import           System.IO                   (IOMode (..), hPutStrLn, stderr,
                                              withFile)

data FreqSpec = FreqNone
    | FreqK Int
    | FreqX Double
    deriving (Show)

data RASOptions = RASOptions
    { _rasBaseDirs       :: [FilePath]
    , _rasJackknifeMode  :: JackknifeMode
    , _rasExcludeChroms  :: [Chrom]
    , _rasPopConfig      :: FilePath
    , _rasMinFreq        :: FreqSpec
    , _rasMaxFreq        :: FreqSpec
    , _rasMaxMissingness :: Double
    , _rasBlockTableFile :: Maybe FilePath
    , _rasF4tableOutFile :: Maybe FilePath
    , _rasMaxSnps        :: Maybe Int
    , _rasNoTransitions  :: Bool
    , _rasBedFile        :: Maybe FilePath
    }
    deriving (Show)

data BlockData = BlockData
    { blockStartPos  :: GenomPos
    , blockEndPos    :: GenomPos
    , blockSiteCount :: [Int]
    , blockVals      :: [[Double]]
    }
    deriving (Show)

data BlockAccumulator = BlockAccumulator
    { blockAccStart     :: Maybe GenomPos
    , blockAccEnd       :: Maybe GenomPos
    , blockAccSiteCount :: VUM.IOVector Int
    , blockAccVals      :: VUM.IOVector Double
    }

runRAS :: RASOptions -> PoseidonIO ()
runRAS rasOpts = do
    -- reading in the configuration file
    PopConfigYamlStruct groupDefs popLefts popRights maybeOutgroup <- readPopConfig (_rasPopConfig rasOpts)
    unless (null groupDefs) . logInfo $ "Found group definitions: " ++ show groupDefs
    logInfo $ "Found left populations: " ++ show popLefts
    logInfo $ "Found right populations: " ++ show popRights
    case maybeOutgroup of
        Nothing -> return ()
        Just o  -> logInfo $ "Found outgroup: " ++ show o

    -- reading in Poseidon packages
    let pacReadOpts = defaultPackageReadOptions {_readOptIgnoreChecksums = True}
    allPackages <- readPoseidonPackageCollection pacReadOpts (_rasBaseDirs rasOpts)
    logInfo $ "Loaded " ++ show (length allPackages) ++ " packages"

    -- if no outgroup is given, set it as empty list
    let outgroupSpec = case maybeOutgroup of
            Nothing -> []
            Just o  -> [o]

    -- check whether all individuals that are needed for the statistics are there, including individuals needed for the adhoc-group definitions in the config file
    let newGroups = map (Group . fst) groupDefs
    let allEntities = nub (concatMap (map underlyingEntity . snd) groupDefs ++
            popLefts ++ popRights ++ outgroupSpec) \\ newGroups
    checkIfAllEntitiesExist allEntities =<< getJointIndividualInfo allPackages

    -- annotate all individuals in all packages with the new adhoc-group definitions where necessary
    let packagesWithNewGroups = addGroupDefs groupDefs allPackages

    -- select only the packages needed for the statistics to be computed
    relevantPackageNames <- determineRelevantPackages (popLefts ++ popRights ++ outgroupSpec) =<< getJointIndividualInfo packagesWithNewGroups
    let relevantPackages = filter (flip elem relevantPackageNames . posPacNameAndVersion) packagesWithNewGroups
    logInfo $ (show . length $ relevantPackages) ++ " relevant packages for chosen statistics identified:"
    mapM_ (logInfo . show . posPacNameAndVersion) relevantPackages

    -- build the main fold, i.e. the thing that does the actual work with the genotype data and computes the RAS statistics (see buildRasFold)
    rasFold <- buildRasFold relevantPackages (_rasMinFreq rasOpts) (_rasMaxFreq rasOpts)
        (_rasMaxMissingness rasOpts) maybeOutgroup popLefts popRights

    -- build a bed-filter if needed
    let bedFilterFunc = case _rasBedFile rasOpts of
            Nothing -> id
            Just fn -> filterThroughBed (readBedFile fn) (genomicPosition . fst)

    -- run the fold and retrieve the block data needed for RAS computations and output
    logA <- envLogAction
    inPlinkPopMode <- envInputPlinkMode
    blockData <- liftIO $ catch (
        runSafeT $ do
            (_, eigenstratProd) <- getJointGenotypeData logA False inPlinkPopMode relevantPackages Nothing
            let eigenstratProdFiltered =
                    bedFilterFunc (eigenstratProd >->
                                    P.filter (chromFilter (_rasExcludeChroms rasOpts)) >->
                                    capNrSnps (_rasMaxSnps rasOpts) >->
                                    filterTransitions (_rasNoTransitions rasOpts))
                eigenstratProdInChunks = case _rasJackknifeMode rasOpts of
                    JackknifePerChromosome  -> chunkEigenstratByChromosome eigenstratProdFiltered
                    JackknifePerN chunkSize -> chunkEigenstratByNrSnps chunkSize eigenstratProdFiltered
            let summaryStatsProd = impurely foldsM rasFold eigenstratProdInChunks
            logWithEnv logA . logInfo $ "performing counts"
            purely P.fold list (summaryStatsProd >-> P.tee (P.map showBlockLogOutput >-> P.toHandle stderr))
        ) (throwIO . PoseidonGenotypeExceptionForward)

    -- outputting and computing results
    logInfo "collating results"

    -- Output for the standard output (a simple table with RAS estimates)
    liftIO . putStrLn . intercalate "\t" $ ["Left", "Right", "Norm", "RAS", "StdErr"]
    forM_ (zip [0..] popLefts) $ \(i, popLeft) ->
        forM_ (zip [0..] popRights) $ \(j, popRight) -> do
            -- get the raw counts for the Jackknife computation
            let counts = [blockSiteCount bd !! i | bd <- blockData]
                vals = [(blockVals bd !! i) !! j | bd <- blockData]
            -- compute jackknife estimate and standard error (see Utils.hs for implementation of the Jackknife)
            let (val, err) = computeJackknifeAdditive counts vals
            liftIO . putStrLn . intercalate "\t" $ [show popLeft, show popRight, show (sum counts), show val, show err]

    -- Compute and output F4 as the pairwise difference of F3. It's only complicated because of the Jackknife
    case _rasF4tableOutFile rasOpts of
        Nothing -> return ()
        Just outFn -> do
            liftIO . withFile outFn WriteMode $ \h -> do
                hPutStrLn h $ intercalate "\t" ["Left1", "Left2", "Right", "Norm1", "Norm2", "RAS-F4", "StdErr", "Z score"]
                -- loop over all possible left1, left2 and rights
                forM_ (zip [0..] popLefts) $ \(l1, popLeft1) ->
                    forM_ (zip [0..] popLefts) $ \(l2, popLeft2) ->
                        when (l1 /= l2) $
                            forM_ (zip [0..] popRights) $ \(r, popRight) -> do
                                let ras1_vals = [(blockVals bd !! l1) !! r | bd <- blockData]
                                    ras2_vals = [(blockVals bd !! l2) !! r | bd <- blockData]
                                    ras1_norms = [blockSiteCount bd !! l1 | bd <- blockData]
                                    ras2_norms = [blockSiteCount bd !! l2 | bd <- blockData]
                                    ras1_full_estimate = weightedAverage ras1_vals (map fromIntegral ras1_norms)
                                    ras2_full_estimate = weightedAverage ras2_vals (map fromIntegral ras2_norms)
                                    rasf4_full_estimate = ras1_full_estimate - ras2_full_estimate
                                    -- compute rasf4 estimates based on one block removed, in turn:
                                let rasf4_minus1_estimates = do
                                        --loop over block-indices to remove each block in turn
                                        removeIndex <- [0.. (length blockData)]
                                        let blockData_minus1 = [bd | (i, bd) <- zip [0..] blockData, i /= removeIndex]
                                        let ras1_vals_minus1 = [(blockVals bd !! l1) !! r | bd <- blockData_minus1]
                                            ras2_vals_minus1 = [(blockVals bd !! l2) !! r | bd <- blockData_minus1]
                                            ras1_norms_minus1 = [blockSiteCount bd !! l1 | bd <- blockData_minus1]
                                            ras2_norms_minus1 = [blockSiteCount bd !! l2 | bd <- blockData_minus1]
                                            ras1_estimate_minus1 = weightedAverage ras1_vals_minus1 (map fromIntegral ras1_norms_minus1)
                                            ras2_estimate_minus1 = weightedAverage ras2_vals_minus1 (map fromIntegral ras2_norms_minus1)
                                        return $ ras1_estimate_minus1 - ras2_estimate_minus1
                                -- compute the block weights needed for Jackknife. Absolute values don't matter, only relative, so we'll just go with the one with more data to assign block weights.
                                let jackknife_block_weights = if sum ras1_norms > sum ras2_norms then ras1_norms else ras2_norms

                                -- compute the jackknife estimate and stderr
                                let (rasf4_jackknife_estimate, rasf4_jackknife_stderr) =
                                        computeJackknifeOriginal rasf4_full_estimate (map fromIntegral jackknife_block_weights) rasf4_minus1_estimates

                                -- print
                                hPutStrLn h . intercalate "\t" $ [show popLeft1, show popLeft2, show popRight, show (sum ras1_norms),
                                    show (sum ras2_norms), show rasf4_jackknife_estimate, show rasf4_jackknife_stderr,
                                    show (rasf4_jackknife_estimate / rasf4_jackknife_stderr)]
    -- optionally output the block file
    case _rasBlockTableFile rasOpts of
        Nothing -> return ()
        Just fn -> liftIO . withFile fn WriteMode $ \h -> do
            hPutStrLn h . intercalate "\t" $ ["Left", "Right", "BlockNr", "StartChrom", "StartPos", "EndChrom", "EndPos", "Norm", "RAS"]
            forM_ (zip [0..] popLefts) $ \(i, popLeft) ->
                forM_ (zip [0..] popRights) $ \(j, popRight) ->
                    forM_ (zip [(0 :: Int)..] blockData) $ \(k, block) ->
                        liftIO . hPutStrLn h . intercalate "\t" $ [show popLeft, show popRight, show k, show (fst . blockStartPos $ block),
                            show (snd . blockStartPos $ block), show (fst . blockEndPos $ block),
                            show (snd . blockEndPos $ block), show (blockSiteCount block !! i), show ((blockVals block !! i) !! j)]
  where
    chromFilter exclusionList (EigenstratSnpEntry chrom _ _ _ _ _, _) = chrom `notElem` exclusionList
    capNrSnps Nothing  = cat
    capNrSnps (Just n) = P.take n
    chunkEigenstratByChromosome = view (groupsBy sameChrom)
    sameChrom (EigenstratSnpEntry chrom1 _ _ _ _ _, _) (EigenstratSnpEntry chrom2 _ _ _ _ _, _) =
        chrom1 == chrom2
    chunkEigenstratByNrSnps chunkSize = view (chunksOf chunkSize)
    showBlockLogOutput block = "computing chunk range " ++ show (blockStartPos block) ++ " - " ++
        show (blockEndPos block) ++ ", size " ++ (show . blockSiteCount) block

weightedAverage :: [Double] -> [Double] -> Double
weightedAverage vals weights =
    let num = sum $ zipWith (*) vals weights
        denom = sum weights
    in  num / denom

readPopConfig :: FilePath -> PoseidonIO PopConfig
readPopConfig fn = do
    bs <- liftIO $ B.readFile fn
    case decodeEither' bs of
        Left err -> liftIO . throwIO $ PopConfigYamlException fn (show err)
        Right x  -> return x

buildRasFold :: (MonadIO m) => [PoseidonPackage] -> FreqSpec -> FreqSpec -> Double -> Maybe PoseidonEntity -> EntitiesList -> EntitiesList -> PoseidonIO (FoldM m (EigenstratSnpEntry, GenoLine) BlockData)
buildRasFold packages minFreq maxFreq maxM maybeOutgroup popLefts popRights = do
    indInfoCollection <- getJointIndividualInfo packages
    ploidyVec <- makePloidyVec . getJointJanno $ packages
    outgroupI <- case maybeOutgroup of
            Nothing -> return []
            Just o  -> resolveUniqueEntityIndices [o] indInfoCollection
    leftI <- sequence [resolveUniqueEntityIndices [l] indInfoCollection | l <- popLefts]
    rightI <- sequence [resolveUniqueEntityIndices [r] indInfoCollection | r <- popRights]
    let nL = length popLefts
        nR = length popRights
    let indivNames = map indInfoName (fst indInfoCollection)
    return $ FoldM (step ploidyVec indivNames outgroupI leftI rightI) (initialise nL nR) extract
  where
    step :: (MonadIO m) => PloidyVec -> [String] -> [Int] -> [[Int]] -> [[Int]] -> BlockAccumulator ->
        (EigenstratSnpEntry, GenoLine) -> m BlockAccumulator
    step ploidyVec indivNames outgroupI leftI rightI (BlockAccumulator maybeStartPos _ counts vals) (EigenstratSnpEntry c p _ _ _ _, genoLine) = do
        let newStartPos = case maybeStartPos of
                Nothing       -> Just (c, p)
                Just (c', p') -> Just (c', p')
        let newEndPos = Just (c, p)
            rAlleleCountPairs = map (computeAlleleCount genoLine ploidyVec indivNames) rightI
            rTotalDerived = sum . map fst $ rAlleleCountPairs
            rTotalHaps = sum $ do
                i <- rightI
                i2 <- i
                case ploidyVec V.! i2 of
                    Haploid -> return 1
                    Diploid -> return 2
            rTotalFreq = (fromIntegral rTotalDerived / fromIntegral rTotalHaps :: Double)
            rTotalNonMissing = sum . map snd $ rAlleleCountPairs
            rMissingness = fromIntegral (rTotalHaps - rTotalNonMissing) / fromIntegral rTotalHaps
        let maybePolarise =
                if null outgroupI
                then Just (rTotalFreq > 0.5) -- repolarise to minor freq if no ougroup is set
                else
                    case computeAlleleCount genoLine ploidyVec indivNames outgroupI of
                        (_, 0) -> Nothing -- outgroup missing
                        (ogNonRef, ogNonMissing) ->
                            if ogNonRef == 0 then Just False else -- outgroup is fixed at reference allele
                                if ogNonRef == ogNonMissing then Just True else -- outgroup is fixed at alternative allele
                                    Nothing -- outgroup is polymorphic
        when (rMissingness <= maxM) $ do
            case maybePolarise of
                Nothing -> return ()
                Just polarise -> do
                    -- update counts
                    forM_ (zip [0..] leftI) $ \(i1, i2s) ->
                        when (any (/= Missing) [genoLine V.! j | j <- i2s]) . liftIO $ VUM.modify counts (+1) i1
                    let rDirectedTotalCount = if polarise then rTotalHaps - rTotalDerived else rTotalDerived
                        rDirectedTotalFreq = fromIntegral rDirectedTotalCount / fromIntegral rTotalHaps
                    let conditionMin = case minFreq of
                            FreqNone -> True
                            FreqK k  -> rDirectedTotalCount >= k
                            FreqX x  -> rDirectedTotalFreq >= x
                        conditionMax = case maxFreq of
                            FreqNone -> True
                            FreqK k  -> rDirectedTotalCount <= k
                            FreqX x  -> rDirectedTotalFreq <= x
                    when (conditionMin && conditionMax) $ do
                        let nR = length popRights
                            leftFreqsNonRef = map (computeAlleleFreq genoLine ploidyVec indivNames) leftI
                            leftFreqs = if polarise then map (fmap (1.0 - )) leftFreqsNonRef else leftFreqsNonRef
                            rightFreqs = do
                                r <- rightI
                                let nrDerived = fst $ computeAlleleCount genoLine ploidyVec indivNames r
                                let n = sum $ do -- we consider missing = reference allele(s)
                                        i <- r
                                        case ploidyVec V.! i of
                                            Haploid -> return 1.0
                                            Diploid -> return 2.0
                                if polarise then
                                    return $ 1.0 - (fromIntegral nrDerived / n)
                                else
                                    return (fromIntegral nrDerived / n)
                            relevantLeftFreqs  = [(i, x) | (i, Just x) <- zip [0..] leftFreqs,  x > 0.0]
                            relevantRightFreqs = [(i, x) | (i,      x) <- zip [0..] rightFreqs, x > 0.0]
                        forM_ relevantLeftFreqs $ \(i, x) ->
                            forM_ relevantRightFreqs $ \(j, y) -> do
                                let index = i * nR + j
                                liftIO $ VUM.modify vals (+(x * y)) index
        return $ BlockAccumulator newStartPos newEndPos counts vals
    initialise :: (MonadIO m) => Int -> Int -> m BlockAccumulator
    initialise nL nR = do
        countVec <- liftIO $ VUM.replicate nL 0
        valVec <- liftIO $ VUM.replicate (nL * nR) 0.0
        return $ BlockAccumulator Nothing Nothing countVec valVec
    extract :: (MonadIO m) => BlockAccumulator -> m BlockData
    extract (BlockAccumulator maybeStartVec maybeEndVec counts vals) = case (maybeStartVec, maybeEndVec) of
        (Just startPos, Just endPos) -> do
            let nLefts = VUM.length counts
                nRights = VUM.length vals `div` nLefts
            countsF <- liftIO $ VU.freeze counts
            valsF <- liftIO $ VU.freeze vals
            let normalisedVals = do
                    i <- [0 .. (nLefts - 1)]
                    return $ do
                        j <- [0 .. (nRights - 1)]
                        let jointIndex = i * nRights + j
                        let val = valsF VU.! jointIndex
                        return $ val / fromIntegral (countsF VU.! i)
            return $ BlockData startPos endPos (VU.toList countsF) normalisedVals
        _ -> error "this should never happen"
