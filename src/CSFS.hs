module CSFS where

import           Utils                       (PopConfig (..), addGroupDefs,
                                              computeAlleleCount,
                                              computeAlleleFreq, readPopConfig)

import           Control.Foldl               (FoldM (..), impurely)
import           Control.Monad               (unless)
import           Control.Monad.IO.Class      (MonadIO, liftIO)
import qualified Data.ByteString             as B
import           Data.List                   (intercalate, nub, (\\))
import qualified Data.Vector                 as V
import qualified Data.Vector.Unboxed         as VU
import qualified Data.Vector.Unboxed.Mutable as VUM
import           Pipes                       (cat, (>->))
import qualified Pipes.Prelude               as P
import           Pipes.Safe                  (runSafeT)
import           Poseidon.EntitiesList       (EntitiesList, PoseidonEntity (..),
                                              conformingEntityIndices,
                                              findNonExistentEntities,
                                              indInfoFindRelevantPackageNames,
                                              underlyingEntity)
import           Poseidon.Package            (PackageReadOptions (..),
                                              PoseidonPackage (..),
                                              defaultPackageReadOptions,
                                              getJointGenotypeData,
                                              getJointIndividualInfo,
                                              readPoseidonPackageCollection)
import           Poseidon.SecondaryTypes     (IndividualInfo (..))
import           SequenceFormats.Eigenstrat  (EigenstratSnpEntry (..), GenoLine)
import           SequenceFormats.Utils       (Chrom (..))
import           System.IO                   (hPutStrLn, stderr)

data CSFSOptions = CSFSOptions
    { _csfsBaseDirs       :: [FilePath]
    , _csfsExcludeChroms  :: [Chrom]
    , _csfsPopConfig      :: FilePath
    , _csfsMaxMissingness :: Double
    , _csfsTableOutFile   :: FilePath
    , _csfsMaxSnps        :: Maybe Int
    }
    deriving (Show)

runCSFS :: CSFSOptions -> IO ()
runCSFS csfsOpts = do
    let pacReadOpts = defaultPackageReadOptions {_readOptStopOnDuplicates = True, _readOptIgnoreChecksums = True}

    PopConfigYamlStruct groupDefs popLefts popRights maybeOutgroup <- readPopConfig (_csfsPopConfig csfsOpts)
    unless (null groupDefs) $ hPutStrLn stderr $ "Found group definitions: " ++ show groupDefs
    hPutStrLn stderr $ "Found left populations: " ++ show popLefts
    hPutStrLn stderr $ "Found right populations: " ++ show popRights
    case maybeOutgroup of
        Nothing -> return ()
        Just o  -> hPutStrLn stderr $ "Found outgroup: " ++ show o

    allPackages <- readPoseidonPackageCollection pacReadOpts (_csfsBaseDirs csfsOpts)
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
        let csfsFold = buildCSFSfold jointIndInfo (_csfsMaxMissingness csfsOpts) maybeOutgroup popLefts popRights

        csfsResults <- runSafeT $ do
            (_, eigenstratProd) <- getJointGenotypeData False False relevantPackages Nothing
            let eigenstratProdFiltered = eigenstratProd >-> P.filter (chromFilter (_csfsExcludeChroms csfsOpts))
                    >-> capNrSnps (_csfsMaxSnps csfsOpts)
            impurely P.foldM csfsFold eigenstratProdFiltered
        print csfsResults
        -- let tableH = ["Left", "Right", "k", "Cumulative", "Norm", "RAS", "StdErr"]
        -- let tableB = do
        --         cumul <- [False, True]
        --         (i, popLeft) <- zip [0..] popLefts
        --         (j, popRight) <- zip [0..] popRights
        --         k <- [0 .. (maxK - 2)]
        --         let counts = [blockSiteCount bd !! i | bd <- blockData]
        --             vals = if cumul then
        --                     [sum [((blockVals bd !! k') !! i) !! j | k' <- [0 .. k]] | bd <- blockData]
        --                 else
        --                     [((blockVals bd !! k) !! i) !! j | bd <- blockData]
        --         let (val, err) = computeJackknife counts vals
        --         return [show popLeft, show popRight, show (k + 2), show cumul, show (sum counts), show val, show err]
        -- let colSpecs = replicate 7 (column expand def def def)
        -- putStrLn $ tableString colSpecs asciiRoundS (titlesH tableH) [rowsG tableB]
        -- withFile (_optTableOutFile rasOpts) WriteMode $ \h -> do
        --     hPutStrLn h $ intercalate "\t" tableH
        --     forM_ tableB $ \row -> hPutStrLn h (intercalate "\t" row)
        -- return ()
  where
    chromFilter exclusionList (EigenstratSnpEntry chrom _ _ _ _ _, _) = chrom `notElem` exclusionList
    capNrSnps Nothing  = cat
    capNrSnps (Just n) = P.take n

buildCSFSfold :: (MonadIO m) => [IndividualInfo] -> Double -> Maybe PoseidonEntity -> EntitiesList -> EntitiesList -> FoldM m (EigenstratSnpEntry, GenoLine) (VU.Vector Int)
buildCSFSfold indInfo maxM maybeOutgroup popLefts popRights =
    let outgroupI = case maybeOutgroup of
            Nothing -> []
            Just o  -> conformingEntityIndices [o] indInfo
        leftI = [conformingEntityIndices [l] indInfo | l <- popLefts]
        rightI = [conformingEntityIndices [r] indInfo | r <- popRights]
        nL = length popLefts
        nR = length popRights
    in  FoldM (step outgroupI leftI rightI) (initialise leftI rightI) extract
  where
    step :: (MonadIO m) => [Int] -> [[Int]] -> [[Int]] -> VUM.IOVector Int -> (EigenstratSnpEntry, GenoLine) -> m (VUM.IOVector Int)
    step outgroupI leftI rightI vals (EigenstratSnpEntry c p _ _ _ _, genoLine) = do
        let alleleCountPairs = map (computeAlleleCount genoLine) rightI
            totalDerived = sum . map fst $ alleleCountPairs
            totalNonMissing = sum . map snd $ alleleCountPairs
            totalHaps = 2 * sum (map length rightI)
            missingness = fromIntegral (totalHaps - totalNonMissing) / fromIntegral totalHaps
        let outgroupFreq = if null outgroupI then Just 0.0 else computeAlleleFreq genoLine outgroupI
        case outgroupFreq of
            Nothing -> return vals
            Just oFreq -> do
                let directedTotalCount = if oFreq < 0.5 then totalDerived else totalHaps - totalDerived
                liftIO $ VUM.modify vals (+1) directedTotalCount
                return vals
    initialise :: (MonadIO m) => [[Int]] -> [[Int]] -> m (VUM.IOVector Int)
    initialise leftI rightI = liftIO $ VUM.replicate (2 * sum (map length rightI) + 1) 0
    extract :: (MonadIO m) => VUM.IOVector Int -> m (VU.Vector Int)
    extract vals = liftIO $ VU.freeze vals
