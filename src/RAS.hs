{-# LANGUAGE OverloadedStrings #-}

module RAS where

import           Utils                       (GenomPos, JackknifeMode (..),
                                              computeAlleleFreq, computeJackknife, PopSpec(..),
                                              popSpecParser, getPopIndices)

import           Control.Exception           (Exception, throwIO)
import           Control.Foldl               (FoldM (..), impurely, list,
                                              purely)
import           Control.Monad               (forM, forM_, when)
import           Control.Monad.Catch         (throwM)
import           Control.Monad.IO.Class      (MonadIO, liftIO)
import           Data.Aeson                  (FromJSON, Object, parseJSON,
                                              withObject, (.:))
import           Data.Aeson.Types            (Parser)
import qualified Data.ByteString             as B
import           Data.List                   (intersect, intercalate)
import           Data.Maybe                  (catMaybes)
import           Data.Text                   (Text)
import qualified Data.Vector                 as V
import qualified Data.Vector.Unboxed         as VU
import qualified Data.Vector.Unboxed.Mutable as VUM

import           Data.Yaml                   (decodeEither')
import           Lens.Family2                (view)

import           Pipes                       ((>->), cat)
import           Pipes.Group                 (chunksOf, foldsM, groupsBy)
import qualified Pipes.Prelude               as P
import           Pipes.Safe                  (runSafeT)
import           Poseidon.Package            (PackageReadOptions (..),
                                              PoseidonPackage (..),
                                              defaultPackageReadOptions,
                                              getIndividuals,
                                              getJointGenotypeData,
                                              readPoseidonPackageCollection)
import           Poseidon.Utils              (PoseidonException (..))
import           SequenceFormats.Eigenstrat  (EigenstratIndEntry (..),
                                              EigenstratSnpEntry (..),
                                              GenoEntry (..), GenoLine)
import           SequenceFormats.Utils       (Chrom (..))
import           System.IO                   (hPutStrLn, stderr, withFile, IOMode(..))
import           Text.Layout.Table          (asciiRoundS, column, def, expand,
                                             rowsG, tableString, titlesH)
import qualified Text.Parsec                 as P

data RASOptions = RASOptions
    { _optBaseDirs       :: [FilePath]
    , _optJackknifeMode  :: JackknifeMode
    , _optExcludeChroms  :: [Chrom]
    , _optPopConfig      :: PopConfig
    , _optMaxCutoff      :: Int
    , _optMaxMissingness :: Double
    , _optTableOutFile   :: FilePath
    , _optMaxSnps        :: Maybe Int
    }
    deriving (Show)

data PopConfig = PopConfigDirect [PopSpec] [PopSpec] PopSpec
    | PopConfigFile FilePath
    deriving (Show)

data PopConfigYamlStruct = PopConfigYamlStruct
    { popConfigLefts    :: [PopSpec]
    , popConfigRights   :: [PopSpec]
    , popConfigOutgroup :: PopSpec
    }

instance FromJSON PopConfigYamlStruct where
    parseJSON = withObject "PopConfigYamlStruct" $ \v -> PopConfigYamlStruct
        <$> parsePopSpecsFromJSON v "popLefts"
        <*> parsePopSpecsFromJSON v "popRights"
        <*> parsePopSpecFromJSON v "outgroup"
      where
        parsePopSpecFromJSON :: Object -> Text -> Parser PopSpec
        parsePopSpecFromJSON v label = do
            popDefString <- v .: label
            case P.runParser popSpecParser () "" popDefString of
                Left err -> fail (show err)
                Right p  -> return p
        parsePopSpecsFromJSON :: Object -> Text -> Parser [PopSpec]
        parsePopSpecsFromJSON v label = do
            popDefStrings <- v .: label
            forM popDefStrings $ \popDefString -> do
                case P.runParser popSpecParser () "" popDefString of
                    Left err -> fail (show err)
                    Right p  -> return p

data RascalException = PopConfigYamlException FilePath String
    deriving (Show)

instance Exception RascalException

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
    allPackages <- readPoseidonPackageCollection pacReadOpts (_optBaseDirs rasOpts)

    hPutStrLn stderr ("Loaded " ++ show (length allPackages) ++ " packages")
    (popLefts, popRights, outgroup) <- case _optPopConfig rasOpts of
        PopConfigDirect pl pr og -> return (pl, pr, og)
        PopConfigFile f          -> readPopConfig f
    hPutStrLn stderr $ "Found left populations: " ++ show popLefts
    hPutStrLn stderr $ "Found right populations: " ++ show popRights
    hPutStrLn stderr $ "Found outgroup: " ++ show outgroup

    let collectedStats = popLefts ++ popRights ++ [outgroup]
    relevantPackages <- findRelevantPackages collectedStats allPackages
    hPutStrLn stderr $ (show . length $ relevantPackages) ++ " relevant packages for chosen statistics identified:"
    forM_ relevantPackages $ \pac -> hPutStrLn stderr (posPacTitle pac)

    blockData <- runSafeT $ do
        (eigenstratIndEntries, eigenstratProd) <- getJointGenotypeData False False relevantPackages Nothing
        let eigenstratProdFiltered = eigenstratProd >-> P.filter (chromFilter (_optExcludeChroms rasOpts))
                >-> capNrSnps (_optMaxSnps rasOpts)
            eigenstratProdInChunks = case _optJackknifeMode rasOpts of
                JackknifePerChromosome  -> chunkEigenstratByChromosome eigenstratProdFiltered
                JackknifePerN chunkSize -> chunkEigenstratByNrSnps chunkSize eigenstratProdFiltered
        rasFold <- case buildRasFold eigenstratIndEntries (_optMaxCutoff rasOpts) (_optMaxMissingness rasOpts) outgroup popLefts popRights of
            Left e  ->  throwM e
            Right f -> return f
        let summaryStatsProd = impurely foldsM rasFold eigenstratProdInChunks
        purely P.fold list (summaryStatsProd >-> P.tee (P.map showBlockLogOutput >-> P.toHandle stderr))
    let maxK = _optMaxCutoff rasOpts
        nLefts = length popLefts
        nRights = length popRights
    let jackknifeEstimates = do
            k <- [0 .. (maxK - 2)]
            return $ do
                i <- [0 .. nLefts]
                return $ do
                    j <- [0 .. nRights]
                    let counts    = [blockSiteCount bd !! i | bd <- blockData]
                    let vals      = [((blockVals bd !! k) !! i) !! j | bd <- blockData]
                        cumulVals = do
                            bd <- blockData
                            return $ sum [((blockVals bd !! k') !! i) !! j | k' <- [0 .. k]]
                    return (computeJackknife counts vals, computeJackknife counts cumulVals)
    let tableH = ["Left", "Right", "k", "Cumulative", "RAS", "StdErr"]
        tableB = do
            cumul <- [False, True]
            (i, popLeft) <- zip [0..] popLefts
            (j, popRight) <- zip [0..] popRights
            k <- [0 .. (maxK - 2)]
            let jne = ((jackknifeEstimates !! k) !! i) !! j
            let (val, err) = if cumul then snd jne else fst jne
            return [show popLeft, show popRight, show (k + 2), show cumul, show val, show err]
    let colSpecs = replicate 6 (column expand def def def)
    putStrLn $ tableString colSpecs asciiRoundS (titlesH tableH) [rowsG tableB]
    withFile (_optTableOutFile rasOpts) WriteMode $ \h -> do
        hPutStrLn h $ intercalate "\t" tableH
        forM_ tableB $ \row -> hPutStrLn h (intercalate "\t" row)
    return ()
  where
    chromFilter exclusionList (EigenstratSnpEntry chrom _ _ _ _ _, _) = chrom `notElem` exclusionList
    capNrSnps Nothing = cat
    capNrSnps (Just n) = P.take n
    chunkEigenstratByChromosome = view (groupsBy sameChrom)
    sameChrom (EigenstratSnpEntry chrom1 _ _ _ _ _, _) (EigenstratSnpEntry chrom2 _ _ _ _ _, _) =
        chrom1 == chrom2
    chunkEigenstratByNrSnps chunkSize = view (chunksOf chunkSize)
    showBlockLogOutput block = "computing chunk range " ++ show (blockStartPos block) ++ " - " ++
        show (blockEndPos block) ++ ", size " ++ (show . blockSiteCount) block

readPopConfig :: FilePath -> IO ([PopSpec], [PopSpec], PopSpec)
readPopConfig fn = do
    bs <- B.readFile fn
    PopConfigYamlStruct pl pr og <- case decodeEither' bs of
        Left err -> throwIO $ PopConfigYamlException fn (show err)
        Right x  -> return x
    return (pl, pr, og)

findRelevantPackages :: [PopSpec] -> [PoseidonPackage] -> IO [PoseidonPackage]
findRelevantPackages popSpecs packages = do
    let indNamesStats   = [ind   | PopSpecInd   ind   <- popSpecs]
        groupNamesStats = [group | PopSpecGroup group <- popSpecs]
    fmap catMaybes . forM packages $ \pac -> do
        inds <- getIndividuals pac
        let indNamesPac   = [ind   | EigenstratIndEntry ind _ _     <- inds]
            groupNamesPac = [group | EigenstratIndEntry _   _ group <- inds]
        if   not (null (indNamesPac `intersect` indNamesStats)) || not (null (groupNamesPac `intersect` groupNamesStats))
        then return (Just pac)
        else return Nothing

buildRasFold :: (MonadIO m) => [EigenstratIndEntry] -> Int -> Double -> PopSpec -> [PopSpec] -> [PopSpec] -> Either PoseidonException (FoldM m (EigenstratSnpEntry, GenoLine) BlockData)
buildRasFold indEntries maxK maxM outgroup popLefts popRights = do
    outgroupI <- getPopIndices indEntries outgroup
    leftI <- mapM (getPopIndices indEntries) popLefts
    rightI <- mapM (getPopIndices indEntries) popRights
    let nL = length popLefts
        nR = length popRights
    return $ FoldM (step outgroupI leftI rightI) (initialise nL nR) extract
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
        -- liftIO $ hPutStrLn stderr (show (totalDerived, totalNonMissing, totalHaps, missingness))
        when (missingness <= maxM) $ do
            let outgroupFreq = if null outgroupI then Just 0.0 else computeAlleleFreq genoLine outgroupI
            case outgroupFreq of
                Nothing -> return ()
                Just oFreq -> do
                    -- update counts
                    forM_ (zip [0..] leftI) $ \(i1, i2s) ->
                        when (any (/= Missing) [genoLine V.! j | j <- i2s]) . liftIO $ VUM.modify counts (+1) i1
                    let directedTotalCount = if oFreq < 0.5 then totalDerived else totalHaps - totalDerived
                    when (directedTotalCount >= 2 && directedTotalCount <= maxK) $ do
                        -- main loop
                        let nL = length popLefts
                            nR = length popRights
                            kIndexOffset = (directedTotalCount - 2) * nL * nR
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
                                let index = kIndexOffset + i * j
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
                        i <- [0 .. nLefts]
                        return $ do
                            j <- [0 .. nRights]
                            let jointIndex = k * nLefts * nRights + i * j
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
