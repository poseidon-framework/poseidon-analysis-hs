{-# LANGUAGE OverloadedStrings #-}

import           FStats                     (FStatSpec (..), FstatsOptions (..),
                                             fStatSpecParser, runFstats)
import           RAS                        (RASOptions (..), runRAS)
import           Utils                      (JackknifeMode (..), runParser)

import           Control.Applicative        ((<|>))
import           Control.Exception          (catch)
import           Data.ByteString.Char8      (pack, splitWith)
import           Data.List                  (intercalate)
import           Data.Version               (showVersion)
import qualified Options.Applicative        as OP
import           Paths_poseidon_analysis_hs (version)
import           Poseidon.PoseidonVersion   (showPoseidonVersion,
                                             validPoseidonVersions)
import           Poseidon.Utils             (PoseidonException (..),
                                             renderPoseidonException)
import           SequenceFormats.Utils      (Chrom (..))
import           System.Exit                (exitFailure)
import           System.IO                  (hPutStrLn, stderr)
import           Text.Read                  (readEither)

data Options = CmdFstats FstatsOptions
    | CmdRAS RASOptions

main :: IO ()
main = do
    hPutStrLn stderr renderVersion
    hPutStrLn stderr ""
    cmdOpts <- OP.customExecParser p optParserInfo
    catch (runCmd cmdOpts) handler
    where
        p = OP.prefs OP.showHelpOnEmpty
        handler :: PoseidonException -> IO ()
        handler e = do
            hPutStrLn stderr $ renderPoseidonException e
            exitFailure

runCmd :: Options -> IO ()
runCmd o = case o of
    CmdFstats opts -> runFstats opts
    CmdRAS opts    -> runRAS opts

optParserInfo :: OP.ParserInfo Options
optParserInfo = OP.info (OP.helper <*> versionOption <*> optParser) (
    OP.briefDesc <>
    OP.progDesc "xerxes is an analysis tool for Poseidon packages. \
                \Report issues here: \
                \https://github.com/poseidon-framework/poseidon-analysis-hs/issues"
    )

versionOption :: OP.Parser (a -> a)
versionOption = OP.infoOption (showVersion version) (OP.long "version" <> OP.help "Show version number")

renderVersion :: String
renderVersion =
    "xerxes v" ++ showVersion version ++ " for poseidon v" ++
    intercalate ", v" (map showPoseidonVersion validPoseidonVersions) ++ "\n" ++
    "https://poseidon-framework.github.io" -- ++ "\n" ++
    --")<(({°> ~ ────E ~ <°}))>("

optParser :: OP.Parser Options
optParser = OP.subparser $
    OP.command "fstats" fstatsOptInfo <>
    OP.command "ras" rasOptInfo
  where
    fstatsOptInfo = OP.info (OP.helper <*> (CmdFstats <$> fstatsOptParser))
        (OP.progDesc "Compute f-statistics on groups and invidiuals within and across Poseidon packages")
    rasOptInfo = OP.info (OP.helper <*> (CmdRAS <$> rasOptParser))
        (OP.progDesc "Compute RAS statistics on groups and individuals within and across Poseidon packages")


fstatsOptParser :: OP.Parser FstatsOptions
fstatsOptParser = FstatsOptions <$> parseBasePaths
                                <*> parseJackknife
                                <*> parseExcludeChroms
                                <*> OP.many parseStatSpecsDirect
                                <*> parseStatSpecsFromFile
                                <*> parseRawOutput
                                <*> parseMaxSnps

parseBasePaths :: OP.Parser [FilePath]
parseBasePaths = OP.some (OP.strOption (OP.long "baseDir" <>
    OP.short 'd' <>
    OP.metavar "DIR" <>
    OP.help "a base directory to search for Poseidon Packages (could be a Poseidon repository)"))

parseJackknife :: OP.Parser JackknifeMode
parseJackknife = OP.option (OP.eitherReader readJackknifeString) (OP.long "jackknife" <> OP.short 'j' <>
    OP.help "Jackknife setting. If given an integer number, this defines the block size in SNPs. \
        \Set to \"CHR\" if you want jackknife blocks defined as entire chromosomes. The default is at 5000 SNPs" <> OP.value (JackknifePerN 5000))
  where
    readJackknifeString :: String -> Either String JackknifeMode
    readJackknifeString s = case s of
        "CHR"  -> Right JackknifePerChromosome
        numStr -> let num = readEither numStr
                  in  case num of
                        Left e  -> Left e
                        Right n -> Right (JackknifePerN n)

parseExcludeChroms :: OP.Parser [Chrom]
parseExcludeChroms = OP.option (map Chrom . splitWith (==',') . pack <$> OP.str) (OP.long "excludeChroms" <> OP.short 'e' <>
    OP.help "List of chromosome names to exclude chromosomes, given as comma-separated \
        \list. Defaults to X, Y, MT, chrX, chrY, chrMT, 23,24,90" <> OP.value [Chrom "X", Chrom "Y", Chrom "MT",
        Chrom "chrX", Chrom "chrY", Chrom "chrMT", Chrom "23", Chrom "24", Chrom "90"])

parseStatSpecsDirect :: OP.Parser FStatSpec
parseStatSpecsDirect = OP.option (OP.eitherReader readStatSpecString) (OP.long "stat" <>
    OP.help "Specify a summary statistic to be computed. Can be given multiple times. \
        \Possible options are: F4(name1, name2, name3, name4), and similarly F3 and F2 stats, \
        \as well as PWM(name1,name2) for pairwise mismatch rates. Group names are by default \
        \matched with group names as indicated in the PLINK or Eigenstrat files in the Poseidon dataset. \
        \You can also specify individual names using the syntax \"<Ind_name>\", so enclosing them \
        \in angular brackets. You can also mix groups and individuals, like in \
        \\"F4(<Ind1>,Group2,Group3,<Ind4>)\". Group or individual names are separated by commas, and a comma \
        \can be followed by any number of spaces, as in some of the examples in this help text.")

parseStatSpecsFromFile :: OP.Parser (Maybe FilePath)
parseStatSpecsFromFile = OP.option (Just <$> OP.str) (OP.long "statFile" <> OP.help "Specify a file with F-Statistics specified \
    \similarly as specified for option --stat. One line per statistics, and no new-line at the end" <> OP.value Nothing)

readStatSpecString :: String -> Either String FStatSpec
readStatSpecString s = case runParser fStatSpecParser () "" s of
    Left p  -> Left (show p)
    Right x -> Right x

parseRawOutput :: OP.Parser Bool
parseRawOutput = OP.switch (
    OP.long "raw" <>
    OP.help "output table as tsv without header. Useful for piping into grep or awk"
    )

rasOptParser :: OP.Parser RASOptions
rasOptParser = RASOptions <$>
    parseBasePaths <*>
    parseJackknife <*>
    parseExcludeChroms <*>
    parsePopConfigFile <*>
    parseMaxCutoff <*>
    parseMaxMissingness <*>
    parseFullTable <*>
    parseTableOutFile <*>
    parseMaxSnps

parsePopConfigFile :: OP.Parser FilePath
parsePopConfigFile = OP.option OP.str (OP.long "popConfigFile" <> OP.help "a file containing the population configuration")

parseMaxCutoff :: OP.Parser (Maybe Int)
parseMaxCutoff = parseK <|> parseNoCutoff
  where
    parseK = OP.option (Just <$> OP.auto) (OP.long "maxAlleleCount" <> OP.short 'k' <>
        OP.help "define a maximal allele-count cutoff for the RAS statistics. " <>
        OP.value (Just 10) <> OP.showDefault)
    parseNoCutoff = OP.flag' Nothing (OP.long "noMaxAlleleCount" <>
        OP.help "switch off the allele count filter. This effectively then mimics Outgroup-F3")

parseMaxMissingness :: OP.Parser Double
parseMaxMissingness = OP.option OP.auto (OP.long "maxMissingness" <> OP.short 'm' <>
    OP.help "define a maximal missingness for the right populations in the RAS statistics." <>
    OP.value 0.1 <> OP.showDefault)

parseFullTable :: OP.Parser Bool
parseFullTable = OP.switch (OP.long "fulltable" <> OP.help "If specified, output the full RAS table over all allele counts")

parseTableOutFile :: OP.Parser FilePath
parseTableOutFile = OP.option OP.str (OP.long "tableOutFile" <> OP.short 'f' <>
    OP.help "the file to which results are written as tab-separated file")

parseMaxSnps :: OP.Parser (Maybe Int)
parseMaxSnps = OP.option (Just <$> OP.auto) (OP.long "maxSnps" <>
    OP.help "Stop after a maximum nr of snps has been processed. Useful for short test runs" <>
    OP.value Nothing <> OP.hidden)
