{-# LANGUAGE OverloadedStrings #-}

import           Control.Applicative   (some, (<|>))
import           Control.Exception     (Exception, throwIO)
import           Control.Monad         (forM)
import           Data.Aeson            (FromJSON, parseJSON, withObject, (.:))
import           Data.Aeson.Types      (Parser, Value)
import qualified Data.ByteString       as B
import           Data.ByteString.Char8 (pack, splitWith)
import           Data.Char             (isSpace)
import           Data.Text             (Text)
import           Data.Version          (showVersion)
import           Data.Yaml             (decodeEither')
import qualified Options.Applicative   as OP
import           Paths_RAScal_hs       (version)
import           Poseidon.Package      (PackageReadOptions (..),
                                        defaultPackageReadOptions,
                                        readPoseidonPackageCollection)
import           SequenceFormats.Utils (Chrom (..))
import           System.IO             (hPutStrLn, stderr)
import qualified Text.Parsec           as P
import qualified Text.Parsec.String    as PS
import           Text.Read             (readEither)

data CLIoptions = CLIoptions
    { _optBaseDirs      :: [FilePath]
    , _optJackknifeMode :: JackknifeMode
    , _optExcludeChroms :: [Chrom]
    , _optPopConfig     :: PopConfig
    , _optMinCutoff     :: Int
    , _optMaxCutoff     :: Int
    }
    deriving (Show)

data JackknifeMode = JackknifePerN Int
    | JackknifePerChromosome
    deriving (Show)

data PopConfig = PopConfigDirect [PopDef] [PopDef]
    | PopConfigFile FilePath
    deriving (Show)

type PopDef = [PopComponent]

data PopComponent = PopComponentAdd EntitySpec
    | PopComponentSubtract EntitySpec
    deriving (Show)

-- | A datatype to represent a group or an individual
data EntitySpec = EntitySpecGroup String
    | EntitySpecInd String
    deriving (Eq)

instance Show EntitySpec where
    show (EntitySpecGroup n) = n
    show (EntitySpecInd   n) = "<" ++ n ++ ">"

data PopConfigYamlStruct = PopConfigYamlStruct
    { popConfigLefts  :: [PopDef]
    , popConfigRights :: [PopDef]
    }

instance FromJSON PopConfigYamlStruct where
    parseJSON = withObject "PopConfigYamlStruct" $ \v -> PopConfigYamlStruct
        <$> parsePopDefsFromJSON v "popLefts"
        <*> parsePopDefsFromJSON v "popRights"
      where
        parsePopDefsFromJSON :: Value -> Text -> Parser [PopDef]
        parsePopDefsFromJSON v label = do
            popDefStrings <- v .: label
            forM popDefStrings $ \popDefString -> do
                case parsePopDef popDefString of
                    Left err -> fail err
                    Right p  -> return p

data RascalException = PopConfigYamlException FilePath String
    deriving (Show)

instance Exception RascalException

-- | A helper type to represent a genomic position.
-- type GenomPos = (Chrom, Int)

main :: IO ()
main = do
    cliOpts <- OP.customExecParser p optParserInfo
    let pacReadOpts = defaultPackageReadOptions {_readOptStopOnDuplicates = True, _readOptIgnoreChecksums = True}
    allPackages <- readPoseidonPackageCollection pacReadOpts (_optBaseDirs cliOpts)
    hPutStrLn stderr ("Loaded " ++ show (length allPackages) ++ " packages")
    (popLefts, popRights) <- case _optPopConfig cliOpts of
        PopConfigDirect pl pr -> return (pl, pr)
        PopConfigFile f       -> readPopConfig f
    return ()
  where
    p = OP.prefs OP.showHelpOnEmpty
    optParserInfo = OP.info (OP.helper <*> versionOption <*> optParser) (
        OP.briefDesc <>
        OP.progDesc "rascal computes RAS statistics for Poseidon-packaged genotype data")
    versionOption = OP.infoOption (showVersion version) (OP.long "version" <> OP.help "Show version")

optParser :: OP.Parser CLIoptions
optParser = CLIoptions <$> parseBasePaths
                       <*> parseJackknife
                       <*> parseExcludeChroms
                       <*> parsePopConfig
                       <*> parseMinCutoff
                       <*> parseMaxCutoff

parseBasePaths :: OP.Parser [FilePath]
parseBasePaths = OP.some (OP.strOption (OP.long "baseDir" <>
    OP.short 'd' <>
    OP.metavar "DIR" <>
    OP.help "a base directory to search for Poseidon Packages (could be a Poseidon repository)"))

parseJackknife :: OP.Parser JackknifeMode
parseJackknife = OP.option (OP.eitherReader readJackknifeString) (OP.long "jackknife" <> OP.short 'j' <>
    OP.help "Jackknife setting. If given an integer number, this defines the block size in SNPs. \
        \Set to \"CHR\" if you want jackknife blocks defined as entire chromosomes. The default is at 5000 SNPs" <> OP.value JackknifePerChromosome)
  where
    readJackknifeString :: String -> Either String JackknifeMode
    readJackknifeString s = case s of
        "CHR"  -> Right JackknifePerChromosome
        numStr -> let num = readEither numStr
                  in  case num of
                        Left e  -> Left e
                        Right n -> Right (JackknifePerN n)

parseExcludeChroms :: OP.Parser [Chrom]
parseExcludeChroms = OP.option (map Chrom . splitWith (==',') . pack <$> OP.str)
    (OP.long "excludeChroms" <> OP.short 'e' <>
    OP.help "List of chromosome names to exclude chromosomes, given as comma-separated \
        \list. Defaults to X, Y, MT, chrX, chrY, chrMT, 23,24,90" <> OP.value [Chrom "X", Chrom "Y", Chrom "MT",
        Chrom "chrX", Chrom "chrY", Chrom "chrMT", Chrom "23", Chrom "24", Chrom "90"])

parsePopConfig :: OP.Parser PopConfig
parsePopConfig = parsePopConfigDirect <|> parsePopConfigFile
  where
    parsePopConfigDirect = PopConfigDirect <$> some parseLeftPop <*> some parseRightPop
    parsePopConfigFile = PopConfigFile <$> OP.option OP.str (OP.long "popConfigFile" <>
        OP.help "a file containing the population configuration")

parseLeftPop :: OP.Parser PopDef
parseLeftPop = OP.option (OP.eitherReader parsePopDef) (OP.long "popLeft" <> OP.short 'l' <>
    OP.help "Define a left population. can be given multiple times. A single population can be defined in their simplest form my just entering a group label, such as \"French\". \
    \A single individual can be entered within angular brackets, such as \"<I123>\". More complex group definitions can \
    \involve multiple groups or individuals that are added or subtracted, using a comma-separated list of entities \
    \(groups or individuals), and using the \"!\" symbol to mark an entity to exclude from the definition. \
    \Example: \"French,!<I1234>,!<I1235>,Spanish\". Here, French is added as group, then two individuals are removed, \
    \and then Spanish is added. These operations are always executed in the order they appear in the definition. \
    \Note it is also possible to define completely new groups by adding up specific \
    \individuals, such as \"<I123>,<I124>,<I125>\". Note: In bash or zsh, you need to surround group definitions \
    \using single quotes!")

parseRightPop :: OP.Parser PopDef
parseRightPop = OP.option (OP.eitherReader parsePopDef) (OP.long "popRight" <> OP.short 'r' <>
    OP.help "Define a right population. can be given multiple times. The same rules for complex compositions \
    \apply as with --popLeft, see above.")

popDefHelpStr :: String
popDefHelpStr = ""

parsePopDef :: String -> Either String PopDef
parsePopDef s = case P.runParser popDefParser () "" s of
    Left p  -> Left (show p)
    Right x -> Right x

popDefParser :: PS.Parser PopDef
popDefParser = (componentParserSubtract <|> componentParserAdd) `P.sepBy` (P.char ',')
  where
    componentParserSubtract = P.char '!' *> (PopComponentSubtract <$> componentEntityParser)
    componentParserAdd = PopComponentAdd <$> componentEntityParser
    componentEntityParser = entityIndParser <|> entityGroupParser
    entityIndParser = EntitySpecInd <$> P.between (P.char '<') (P.char '>') parseName
    entityGroupParser = EntitySpecGroup <$> parseName
    parseName = P.many1 (P.satisfy (\c -> not (isSpace c || c `elem` [',', '<', '>'])))

parseMinCutoff :: OP.Parser Int
parseMinCutoff = OP.option OP.auto (OP.long "minCutoff" <>
    OP.help "define a minimal allele-count cutoff for the RAS statistics. " <>
    OP.value 0 <> OP.showDefault)

parseMaxCutoff :: OP.Parser Int
parseMaxCutoff = OP.option OP.auto (OP.long "maxCutoff" <>
    OP.help "define a maximal allele-count cutoff for the RAS statistics. " <>
    OP.value 10 <> OP.showDefault)

readPopConfig :: FilePath -> IO ([PopDef], [PopDef])
readPopConfig fn = do
    bs <- B.readFile fn
    PopConfigYamlStruct pl pr <- case decodeEither' bs of
        Left err -> throwIO $ PopConfigYamlException fn (show err)
        Right x  -> return x
    return (pl, pr)

-- data BlockData = BlockData
--     { blockStartPos  :: GenomPos
--     , blockEndPos    :: GenomPos
--     , blockSiteCount :: Int
--     , blockVal       :: Double
--     }
--     deriving (Show)

-- statSpecsFold :: [EigenstratIndEntry] -> [FStatSpec] -> Either PoseidonException (Fold (EigenstratSnpEntry, GenoLine) [BlockData])
-- statSpecsFold indEntries fStatSpecs = do
--     listOfFolds <- mapM (statSpecFold indEntries) fStatSpecs
--     return $ sequenceA listOfFolds

-- statSpecFold :: [EigenstratIndEntry] -> FStatSpec -> Either PoseidonException (Fold (EigenstratSnpEntry, GenoLine) BlockData)
-- statSpecFold iE fStatSpec = do
--     fStat <- case fStatSpec of
--         F4Spec  a b c d -> F4  <$> getPopIndices iE a <*> getPopIndices iE b <*> getPopIndices iE c <*> getPopIndices iE d
--         F3Spec  a b c   -> F3  <$> getPopIndices iE a <*> getPopIndices iE b <*> getPopIndices iE c
--         F2Spec  a b     -> F2  <$> getPopIndices iE a <*> getPopIndices iE b
--         PWMspec a b     -> PWM <$> getPopIndices iE a <*> getPopIndices iE b
--     return $ Fold (step fStat) initialize extract
--   where
--     step :: FStat -> (Maybe GenomPos, Maybe GenomPos, Int, Double) ->
--         (EigenstratSnpEntry, GenoLine) -> (Maybe GenomPos, Maybe GenomPos, Int, Double)
--     step fstat (maybeStartPos, _, count, val) (EigenstratSnpEntry c p _ _ _ _, genoLine) =
--         let newStartPos = case maybeStartPos of
--                 Nothing       -> Just (c, p)
--                 Just (c', p') -> Just (c', p')
--             newEndPos = Just (c, p)
--         in  case computeFStat fstat genoLine of
--                 Just v  -> (newStartPos, newEndPos, count + 1, val + v)
--                 Nothing -> (newStartPos, newEndPos, count + 1, val)
--     initialize :: (Maybe GenomPos, Maybe GenomPos, Int, Double)
--     initialize = (Nothing, Nothing, 0, 0.0)
--     extract :: (Maybe GenomPos, Maybe GenomPos, Int, Double) -> BlockData
--     extract (maybeStartPos, maybeEndPos, count, totalVal) =
--         let Just startPos = maybeStartPos
--             Just endPos = maybeEndPos
--         in  BlockData startPos endPos count (totalVal / fromIntegral count)

-- computeFStat :: FStat -> GenoLine -> Maybe Double
-- computeFStat fStat gL = case fStat of
--     (F4  aI bI cI dI) -> computeF4  <$> computeFreq gL aI <*> computeFreq gL bI <*> computeFreq gL cI <*> computeFreq gL dI
--     (F3  aI bI cI   ) -> computeF3  <$> computeFreq gL aI <*> computeFreq gL bI <*> computeFreq gL cI
--     (F2  aI bI      ) -> computeF2  <$> computeFreq gL aI <*> computeFreq gL bI
--     (PWM aI bI      ) -> computePWM <$> computeFreq gL aI <*> computeFreq gL bI
--   where
--     computeF4  a b c d = (a - b) * (c - d)
--     computeF3  a b c   = (c - a) * (c - b)
--     computeF2  a b     = (a - b) * (a - b)
--     computePWM a b     = a * (1.0 - b) + (1.0 - a) * b

-- computeFreq :: GenoLine -> [Int] -> Maybe Double
-- computeFreq line indices =
--     let nrNonMissing = length . filter (/=Missing) . map (line !) $ indices
--         nrDerived = sum $ do
--             i <- indices
--             case line ! i of
--                 HomRef  -> return (0 :: Integer)
--                 Het     -> return 1
--                 HomAlt  -> return 2
--                 Missing -> return 0
--     in  if nrNonMissing > 0
--         then Just (fromIntegral nrDerived / fromIntegral nrNonMissing)
--         else Nothing

-- getPopIndices :: [EigenstratIndEntry] -> PopSpec -> Either PoseidonException [Int]
-- getPopIndices indEntries popSpec =
--     let ret = do
--             (i, EigenstratIndEntry indName _ popName) <- zip [0..] indEntries
--             True <- case popSpec of
--                 PopSpecGroup name -> return (name == popName)
--                 PopSpecInd   name -> return (name == indName)
--             return i
--     in  if (null ret)
--         then case popSpec of
--             PopSpecGroup n -> Left $ PoseidonIndSearchException ("Group name " ++ n ++ " not found")
--             PopSpecInd   n -> Left $ PoseidonIndSearchException ("Individual name " ++ n ++ " not found")
--         else Right ret

-- | The main function running the FStats command.
-- runFstats :: FstatsOptions -> IO ()
-- runFstats (FstatsOptions baseDirs jackknifeMode exclusionList statSpecsDirect maybeStatSpecsFile rawOutput) = do
--     -- load packages --
--     allPackages <- readPoseidonPackageCollection pacReadOpts baseDirs
--     statSpecsFromFile <- case maybeStatSpecsFile of
--         Nothing -> return []
--         Just f  -> readStatSpecsFromFile f
--     let statSpecs = statSpecsFromFile ++ statSpecsDirect
--     if null statSpecs then
--         hPutStrLn stderr $ "No statistics to be computed"
--     else do
--         let collectedStats = collectStatSpecGroups statSpecs
--         relevantPackages <- findRelevantPackages collectedStats allPackages
--         hPutStrLn stderr $ (show . length $ relevantPackages) ++ " relevant packages for chosen statistics identified:"
--         forM_ relevantPackages $ \pac -> hPutStrLn stderr (posPacTitle pac)
--         hPutStrLn stderr $ "Computing stats " ++ show statSpecs
--         blockData <- runSafeT $ do
--             (eigenstratIndEntries, eigenstratProd) <- getJointGenotypeData False False relevantPackages
--             let eigenstratProdFiltered = eigenstratProd >-> P.filter chromFilter
--                 eigenstratProdInChunks = case jackknifeMode of
--                     JackknifePerChromosome  -> chunkEigenstratByChromosome eigenstratProdFiltered
--                     JackknifePerN chunkSize -> chunkEigenstratByNrSnps chunkSize eigenstratProdFiltered
--             statsFold <- case statSpecsFold eigenstratIndEntries statSpecs of
--                 Left e  ->  throwM e
--                 Right f -> return f
--             let summaryStatsProd = purely folds statsFold eigenstratProdInChunks
--             purely P.fold list (summaryStatsProd >-> P.tee (P.map showBlockLogOutput >-> P.toHandle stderr))
--         let jackknifeEstimates = [computeJackknife (map blockSiteCount blocks) (map blockVal blocks) | blocks <- transpose blockData]
--             colSpecs = replicate 4 (column expand def def def)
--             tableH = ["Statistic", "Estimate", "StdErr", "Z score"]
--             tableB = do
--                 (fstat, result) <- zip statSpecs jackknifeEstimates
--                 return [show fstat, show (fst result), show (snd result), show (fst result / snd result)]
--         if   rawOutput
--         then do
--             putStrLn $ intercalate "\t" tableH
--             forM_ tableB $ \row -> putStrLn (intercalate "\t" row)
--         else putStrLn $ tableString colSpecs asciiRoundS (titlesH tableH) [rowsG tableB]
--   where
--     chromFilter (EigenstratSnpEntry chrom _ _ _ _ _, _) = chrom `notElem` exclusionList
--     chunkEigenstratByChromosome = view (groupsBy sameChrom)
--     sameChrom (EigenstratSnpEntry chrom1 _ _ _ _ _, _) (EigenstratSnpEntry chrom2 _ _ _ _ _, _) =
--         chrom1 == chrom2
--     chunkEigenstratByNrSnps chunkSize = view (chunksOf chunkSize)
--     showBlockLogOutput blocks = "computing chunk range " ++ show (blockStartPos (head blocks)) ++ " - " ++
--         show (blockEndPos (head blocks)) ++ ", size " ++ (show . blockSiteCount . head) blocks ++ ", values " ++
--         (show . map blockVal) blocks

-- readStatSpecsFromFile :: FilePath -> IO [FStatSpec]
-- readStatSpecsFromFile statSpecsFile = do
--     let multiFstatSpecParser = fStatSpecParser `P.sepBy1` (P.newline *> P.spaces)
--     eitherParseResult <- P.parseFromFile (P.spaces *> multiFstatSpecParser <* P.spaces) statSpecsFile
--     case eitherParseResult of
--         Left err -> throwIO (PoseidonFStatsFormatException (show err))
--         Right r  -> return r

-- computeJackknife :: [Int] -> [Double] -> (Double, Double)
-- computeJackknife weights values =
--     let weights'    = map fromIntegral weights
--         sumWeights  = sum weights'
--         g           = fromIntegral (length weights)
--         theta       = sum [mj * val | (mj, val) <- zip weights' values] / sumWeights
--         sigmaSquare = sum [mj * (val - theta) ^ (2 :: Int) / (sumWeights - mj) | (mj, val) <- zip weights' values] / g
--     in  (theta, sqrt sigmaSquare)

-- collectStatSpecGroups :: [FStatSpec] -> [PopSpec]
-- collectStatSpecGroups statSpecs = nub . concat $ do
--     stat <- statSpecs
--     case stat of
--         F4Spec  a b c d -> return [a, b, c, d]
--         F3Spec  a b c   -> return [a, b, c]
--         F2Spec  a b     -> return [a, b]
--         PWMspec a b     -> return [a, b]

-- findRelevantPackages :: [PopSpec] -> [PoseidonPackage] -> IO [PoseidonPackage]
-- findRelevantPackages popSpecs packages = do
--     let indNamesStats   = [ind   | PopSpecInd   ind   <- popSpecs]
--         groupNamesStats = [group | PopSpecGroup group <- popSpecs]
--     fmap catMaybes . forM packages $ \pac -> do
--         inds <- getIndividuals pac
--         let indNamesPac   = [ind   | EigenstratIndEntry ind _ _     <- inds]
--             groupNamesPac = [group | EigenstratIndEntry _   _ group <- inds]
--         if   length (intersect indNamesPac indNamesStats) > 0 || length (intersect groupNamesPac groupNamesStats) > 0
--         then return (Just pac)
--         else return Nothing
