{-# LANGUAGE OverloadedStrings #-}

module RAS where

import           Utils                       (GenomPos, JackknifeMode (..),
                                              computeAlleleFreq, computeAlleleCount,
                                              computeJackknifeAdditive, PopConfig(..), GroupDef, XerxesException(..))

import           Control.Exception           (throwIO)
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
import           Poseidon.EntitiesList       (EntitiesList, PoseidonEntity (..),
                                              conformingEntityIndices,
                                              findNonExistentEntities,
                                              indInfoConformsToEntitySpec,
                                              indInfoFindRelevantPackageNames,
                                              underlyingEntity)
import           Poseidon.Package            (PackageReadOptions (..),
                                              PoseidonPackage (..),
                                              defaultPackageReadOptions,
                                              getJointGenotypeData,
                                              getJointIndividualInfo,
                                              readPoseidonPackageCollection)
import           Poseidon.SecondaryTypes     (IndividualInfo (..))
import           SequenceFormats.Eigenstrat  (EigenstratSnpEntry (..),
                                              GenoEntry (..), GenoLine)
import           SequenceFormats.Utils       (Chrom (..))
import           System.IO                   (IOMode (..), hPutStrLn, stderr,
                                              withFile)
import           Text.Layout.Table           (asciiRoundS, column, def, expand,
                                              rowsG, tableString, titlesH)

data FreqSpec = FreqNone | FreqK Int | FreqX Double deriving (Show)

data RASOptions = RASOptions
    { _rasBaseDirs       :: [FilePath]
    , _rasJackknifeMode  :: JackknifeMode
    , _rasExcludeChroms  :: [Chrom]
    , _rasPopConfig      :: FilePath
    , _rasMinFreq        :: FreqSpec
    , _rasMaxFreq        :: FreqSpec
    , _rasMaxMissingness :: Double
    , _rasBlockTableFile :: Maybe FilePath
    , _rasTableOutFile   :: Maybe FilePath
    , _rasMaxSnps        :: Maybe Int
    }
    deriving (Show)

data BlockData = BlockData
    { blockStartPos  :: GenomPos
    , blockEndPos    :: GenomPos
    , blockSiteCount :: [Int]
    , blockVals      :: [[Double]]
    }
    deriving (Show)

runRAS :: RASOptions -> IO ()
runRAS rasOpts = do
    let pacReadOpts = defaultPackageReadOptions {_readOptStopOnDuplicates = True, _readOptIgnoreChecksums = True}

    PopConfigYamlStruct groupDefs popLefts popRights maybeOutgroup <- readPopConfig (_rasPopConfig rasOpts)
    unless (null groupDefs) $ hPutStrLn stderr $ "Found group definitions: " ++ show groupDefs
    hPutStrLn stderr $ "Found left populations: " ++ show popLefts
    hPutStrLn stderr $ "Found right populations: " ++ show popRights
    case maybeOutgroup of
        Nothing -> return ()
        Just o  -> hPutStrLn stderr $ "Found outgroup: " ++ show o

    allPackages <- readPoseidonPackageCollection pacReadOpts (_rasBaseDirs rasOpts)
    hPutStrLn stderr ("Loaded " ++ show (length allPackages) ++ " packages")

    let outgroupSpec = case maybeOutgroup of
            Nothing -> []
            Just o  -> [o]
    let newGroups = map (Group . fst) groupDefs
    let allEntities = nub (concatMap (map underlyingEntity . snd) groupDefs ++
            popLefts ++ popRights ++ outgroupSpec) \\ newGroups

    let jointIndInfoAll = getJointIndividualInfo allPackages
    let missingEntities = findNonExistentEntities allEntities jointIndInfoAll
    if not. null $ missingEntities then
        hPutStrLn stderr $ "The following entities couldn't be found: " ++ (intercalate ", " . map show $ missingEntities)
    else do
        let jointIndInfoWithNewGroups = addGroupDefs groupDefs jointIndInfoAll
            relevantPackageNames = indInfoFindRelevantPackageNames (popLefts ++ popRights ++ outgroupSpec) jointIndInfoWithNewGroups
        let relevantPackages = filter (flip elem relevantPackageNames . posPacTitle) allPackages
        hPutStrLn stderr $ (show . length $ relevantPackages) ++ " relevant packages for chosen statistics identified:"
        mapM_ (hPutStrLn stderr . posPacTitle) relevantPackages

        let jointIndInfo = addGroupDefs groupDefs . getJointIndividualInfo $ relevantPackages
        let rasFold = buildRasFold jointIndInfo (_rasMinFreq rasOpts) (_rasMaxFreq rasOpts)
                (_rasMaxMissingness rasOpts) maybeOutgroup popLefts popRights
        blockData <- runSafeT $ do
            (_, eigenstratProd) <- getJointGenotypeData False False relevantPackages Nothing
            let eigenstratProdFiltered = eigenstratProd >-> P.filter (chromFilter (_rasExcludeChroms rasOpts))
                    >-> capNrSnps (_rasMaxSnps rasOpts)
                eigenstratProdInChunks = case _rasJackknifeMode rasOpts of
                    JackknifePerChromosome  -> chunkEigenstratByChromosome eigenstratProdFiltered
                    JackknifePerN chunkSize -> chunkEigenstratByNrSnps chunkSize eigenstratProdFiltered
            let summaryStatsProd = impurely foldsM rasFold eigenstratProdInChunks
            liftIO $ hPutStrLn stderr "performing counts"
            purely P.fold list (summaryStatsProd >-> P.tee (P.map showBlockLogOutput >-> P.toHandle stderr))
        liftIO $ hPutStrLn stderr "collating results"
        let tableH = ["Left", "Right", "Norm", "RAS", "StdErr"]
        let tableB = do
                (i, popLeft) <- zip [0..] popLefts
                (j, popRight) <- zip [0..] popRights
                let counts = [blockSiteCount bd !! i | bd <- blockData]
                    vals = [(blockVals bd !! i) !! j | bd <- blockData]
                let (val, err) = computeJackknifeAdditive counts vals
                return [show popLeft, show popRight, show (sum counts), show val, show err]
        let colSpecs = replicate 7 (column expand def def def)
        putStrLn $ tableString colSpecs asciiRoundS (titlesH tableH) [rowsG tableB]
        case _rasTableOutFile rasOpts of
            Nothing -> return ()
            Just outFn -> do
                withFile outFn WriteMode $ \h -> do
                    hPutStrLn h $ intercalate "\t" tableH
                    forM_ tableB $ \row -> hPutStrLn h (intercalate "\t" row)
        case _rasBlockTableFile rasOpts of
            Nothing -> return ()
            Just fn -> withFile fn WriteMode $ \h -> do
                hPutStrLn h $ intercalate "\t" ["Left", "Right", "BlockNr", "StartChrom", "StartPos", "EndChrom", "EndPos", "Norm", "RAS"]
                let rows = do
                        (i, popLeft) <- zip [0..] popLefts
                        (j, popRight) <- zip [0..] popRights
                        (k, block) <- zip [(0 :: Int)..] blockData
                        return [show popLeft, show popRight, show k, show (fst . blockStartPos $ block), 
                            show (snd . blockStartPos $ block), show (fst . blockEndPos $ block), 
                            show (snd . blockEndPos $ block), show (blockSiteCount block !! i), show ((blockVals block !! i) !! j)]
                forM_ rows $ \row -> hPutStrLn h (intercalate "\t" row)
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

buildRasFold :: (MonadIO m) => [IndividualInfo] -> FreqSpec -> FreqSpec -> Double -> Maybe PoseidonEntity -> EntitiesList -> EntitiesList -> FoldM m (EigenstratSnpEntry, GenoLine) BlockData
buildRasFold indInfo minFreq maxFreq maxM maybeOutgroup popLefts popRights =
    let outgroupI = case maybeOutgroup of
            Nothing -> []
            Just o  -> conformingEntityIndices [o] indInfo
        leftI = [conformingEntityIndices [l] indInfo | l <- popLefts]
        rightI = [conformingEntityIndices [r] indInfo | r <- popRights]
        nL = length popLefts
        nR = length popRights
    in  FoldM (step outgroupI leftI rightI) (initialise nL nR) extract
  where
    step :: (MonadIO m) => [Int] -> [[Int]] -> [[Int]] -> (Maybe GenomPos, Maybe GenomPos, VUM.IOVector Int, VUM.IOVector Double) ->
        (EigenstratSnpEntry, GenoLine) -> m (Maybe GenomPos, Maybe GenomPos, VUM.IOVector Int, VUM.IOVector Double)
    step outgroupI leftI rightI (maybeStartPos, _, counts, vals) (EigenstratSnpEntry c p _ _ _ _, genoLine) = do
        let newStartPos = case maybeStartPos of
                Nothing       -> Just (c, p)
                Just (c', p') -> Just (c', p')
        let newEndPos = Just (c, p)
            alleleCountPairs = map (computeAlleleCount genoLine) rightI
            totalDerived = sum . map fst $ alleleCountPairs
            totalNonMissing = sum . map snd $ alleleCountPairs
            totalHaps = 2 * sum (map length rightI)
            missingness = fromIntegral (totalHaps - totalNonMissing) / fromIntegral totalHaps
        -- liftIO $ hPrint stderr (totalDerived, totalNonMissing, totalHaps, missingness)
        when (missingness <= maxM) $ do
            let outgroupFreq = if null outgroupI then Just 0.0 else computeAlleleFreq genoLine outgroupI
            case outgroupFreq of
                Nothing -> return ()
                Just oFreq -> do
                    -- update counts
                    forM_ (zip [0..] leftI) $ \(i1, i2s) ->
                        when (any (/= Missing) [genoLine V.! j | j <- i2s]) . liftIO $ VUM.modify counts (+1) i1
                    let directedTotalCount = if oFreq < 0.5 then totalDerived else totalHaps - totalDerived
                        directedTotalFreq = fromIntegral directedTotalCount / fromIntegral totalHaps
                    -- liftIO $ hPrint stderr (directedTotalCount, totalDerived, totalNonMissing, totalHaps, missingness)
                    let conditionMin = case minFreq of
                            FreqNone -> True                        
                            FreqK k -> directedTotalCount >= k
                            FreqX x -> directedTotalFreq >= x
                        conditionMax = case maxFreq of
                            FreqNone -> True                        
                            FreqK k -> directedTotalCount <= k
                            FreqX x -> directedTotalFreq <= x
                    when (conditionMin && conditionMax) $ do
                        -- liftIO $ hPrint stderr (directedTotalCount, totalDerived, totalNonMissing, totalHaps, missingness)
                        -- main loop
                        let nR = length popRights
                            leftFreqs = map (computeAlleleFreq genoLine) leftI
                            rightFreqs = do
                                r <- rightI
                                let nrDerived = fst $ computeAlleleCount genoLine r
                                let n = 2 * length r
                                return (fromIntegral nrDerived / fromIntegral n)
                            relevantLeftFreqs  = [(i, x) | (i, Just x) <- zip [0..] leftFreqs,  x > 0.0]
                            relevantRightFreqs = [(i, x) | (i,      x) <- zip [0..] rightFreqs, x > 0.0]
                        forM_ relevantLeftFreqs $ \(i, x) ->
                            forM_ relevantRightFreqs $ \(j, y) -> do
                                let index = i * nR + j
                                liftIO $ VUM.modify vals (+(x * y)) index
        return (newStartPos, newEndPos, counts, vals)
    initialise :: (MonadIO m) => Int -> Int -> m (Maybe GenomPos, Maybe GenomPos, VUM.IOVector Int, VUM.IOVector Double)
    initialise nL nR = do
        countVec <- liftIO $ VUM.replicate nL 0
        valVec <- liftIO $ VUM.replicate (nL * nR) 0.0
        return (Nothing, Nothing, countVec, valVec)
    extract :: (MonadIO m) => (Maybe GenomPos, Maybe GenomPos, VUM.IOVector Int, VUM.IOVector Double) ->
            m BlockData
    extract (maybeStartVec, maybeEndVec, counts, vals) = case (maybeStartVec, maybeEndVec) of
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
