module CSFS where

import Utils (PopConfig(..))

import           SequenceFormats.Utils (Chrom (..))
import qualified Pipes.Prelude               as P
import           Poseidon.Package            (PackageReadOptions (..))

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
        let csfsFold = buildCSFSfold jointIndInfo (_optMaxCutoff rasOpts) (_optMaxMissingness rasOpts) maybeOutgroup popLefts popRights

        csfsResults <- runSafeT $ do
            (_, eigenstratProd) <- getJointGenotypeData False False relevantPackages Nothing
            let eigenstratProdFiltered = eigenstratProd >-> P.filter (chromFilter (_optExcludeChroms rasOpts))
                    >-> capNrSnps (_optMaxSnps rasOpts)
            impurely foldM csfsFold eigenstratProdFiltered
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

buildCSFSfold :: (MonadIO m) => [IndividualInfo] -> Double -> Maybe PoseidonEntity -> EntitiesList -> EntitiesList -> FoldM m (EigenstratSnpEntry, GenoLine) BlockData
buildCSFSfold indInfo maxM maybeOutgroup popLefts popRights =
