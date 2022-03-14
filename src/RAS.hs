{-# LANGUAGE OverloadedStrings #-}

module RAS where

import           Utils                       (GenomPos, JackknifeMode (..),
                                              computeAlleleFreq,
                                              computeJackknife, PopConfig(..), GroupDef, XerxesException(..))

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
                                              withFile, hPrint)
import           Text.Layout.Table           (asciiRoundS, column, def, expand,
                                              rowsG, tableString, titlesH)

data RASOptions = RASOptions
    { _optBaseDirs       :: [FilePath]
    , _optJackknifeMode  :: JackknifeMode
    , _optExcludeChroms  :: [Chrom]
    , _optPopConfig      :: FilePath
    , _optMaxCutoff      :: Int
    , _optMaxMissingness :: Double
    , _optTableOutFile   :: FilePath
    , _optMaxSnps        :: Maybe Int
    }
    deriving (Show)

data BlockData = BlockData
    { blockStartPos  :: GenomPos
    , blockEndPos    :: GenomPos
    , blockSiteCount :: [Int]
    , blockVals      :: [[[Double]]]
    }
    deriving (Show)

runRAS :: RASOptions -> IO ()
runRAS rasOpts = do
    let pacReadOpts = defaultPackageReadOptions {_readOptStopOnDuplicates = True, _readOptIgnoreChecksums = True}

    PopConfigYamlStruct groupDefs popLefts popRights maybeOutgroup <- readPopConfig (_optPopConfig rasOpts)
    unless (null groupDefs) $ hPutStrLn stderr $ "Found group definitions: " ++ show groupDefs
    hPutStrLn stderr $ "Found left populations: " ++ show popLefts
    hPutStrLn stderr $ "Found right populations: " ++ show popRights
    case maybeOutgroup of
        Nothing -> return ()
        Just o  -> hPutStrLn stderr $ "Found outgroup: " ++ show o

    allPackages <- readPoseidonPackageCollection pacReadOpts (_optBaseDirs rasOpts)
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

        blockData <- runSafeT $ do
            (_, eigenstratProd) <- getJointGenotypeData False False relevantPackages Nothing
            let eigenstratProdFiltered = eigenstratProd >-> P.filter (chromFilter (_optExcludeChroms rasOpts))
                    >-> capNrSnps (_optMaxSnps rasOpts)
                eigenstratProdInChunks = case _optJackknifeMode rasOpts of
                    JackknifePerChromosome  -> chunkEigenstratByChromosome eigenstratProdFiltered
                    JackknifePerN chunkSize -> chunkEigenstratByNrSnps chunkSize eigenstratProdFiltered
            let rasFold = buildRasFold jointIndInfo (_optMaxCutoff rasOpts) (_optMaxMissingness rasOpts) maybeOutgroup popLefts popRights
            let summaryStatsProd = impurely foldsM rasFold eigenstratProdInChunks
            purely P.fold list (summaryStatsProd >-> P.tee (P.map showBlockLogOutput >-> P.toHandle stderr))
        let maxK = _optMaxCutoff rasOpts
        let tableH = ["Left", "Right", "k", "Cumulative", "Norm", "RAS", "StdErr"]
        let tableB = do
                cumul <- [False, True]
                (i, popLeft) <- zip [0..] popLefts
                (j, popRight) <- zip [0..] popRights
                k <- [0 .. (maxK - 2)]
                let counts = [blockSiteCount bd !! i | bd <- blockData]
                    vals = if cumul then
                            [sum [((blockVals bd !! k') !! i) !! j | k' <- [0 .. k]] | bd <- blockData]
                        else
                            [((blockVals bd !! k) !! i) !! j | bd <- blockData]
                let (val, err) = computeJackknife counts vals
                return [show popLeft, show popRight, show (k + 2), show cumul, show (sum counts), show val, show err]
        let colSpecs = replicate 7 (column expand def def def)
        putStrLn $ tableString colSpecs asciiRoundS (titlesH tableH) [rowsG tableB]
        withFile (_optTableOutFile rasOpts) WriteMode $ \h -> do
            hPutStrLn h $ intercalate "\t" tableH
            forM_ tableB $ \row -> hPutStrLn h (intercalate "\t" row)
        return ()
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

buildRasFold :: (MonadIO m) => [IndividualInfo] -> Int -> Double -> Maybe PoseidonEntity -> EntitiesList -> EntitiesList -> FoldM m (EigenstratSnpEntry, GenoLine) BlockData
buildRasFold indInfo maxK maxM maybeOutgroup popLefts popRights =
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
                    -- liftIO $ hPrint stderr (directedTotalCount, totalDerived, totalNonMissing, totalHaps, missingness)
                    when (directedTotalCount >= 2 && directedTotalCount <= maxK) $ do
                        liftIO $ hPrint stderr (directedTotalCount, totalDerived, totalNonMissing, totalHaps, missingness)
                        -- main loop
                        let nL = length popLefts
                            nR = length popRights
                            k = directedTotalCount - 2
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
                                let index = k * nL * nR + i * nR + j
                                liftIO $ VUM.modify vals (+(x * y)) index
        return (newStartPos, newEndPos, counts, vals)
    initialise :: (MonadIO m) => Int -> Int -> m (Maybe GenomPos, Maybe GenomPos, VUM.IOVector Int, VUM.IOVector Double)
    initialise nL nR = do
        countVec <- liftIO $ VUM.replicate nL 0
        valVec <- liftIO $ VUM.replicate ((maxK - 1) * nL * nR) 0.0
        return (Nothing, Nothing, countVec, valVec)
    extract :: (MonadIO m) => (Maybe GenomPos, Maybe GenomPos, VUM.IOVector Int, VUM.IOVector Double) ->
            m BlockData
    extract (maybeStartVec, maybeEndVec, counts, vals) = case (maybeStartVec, maybeEndVec) of
        (Just startPos, Just endPos) -> do
            let nLefts = VUM.length counts
                nRights = VUM.length vals `div` (nLefts * (maxK - 1))
            countsF <- liftIO $ VU.freeze counts
            valsF <- liftIO $ VU.freeze vals
            let normalisedVals = do
                    k <- [0 .. (maxK - 2)]
                    return $ do
                        i <- [0 .. (nLefts - 1)]
                        return $ do
                            j <- [0 .. (nRights - 1)]
                            let jointIndex = k * nLefts * nRights + i * nRights + j
                            let val = valsF VU.! jointIndex
                            return $ val / fromIntegral (countsF VU.! i)
            return $ BlockData startPos endPos (VU.toList countsF) normalisedVals
        _ -> error "this should never happen"

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
