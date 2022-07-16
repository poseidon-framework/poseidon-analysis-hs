{-# LANGUAGE OverloadedStrings #-}

import           Poseidon.Analysis.CLI.FStats     (FstatsOptions (..),
                                                   runFstats, runParser)
import           Poseidon.Analysis.CLI.RAS        (FreqSpec (..),
                                                   RASOptions (..), runRAS)
import           Poseidon.Analysis.FStatsConfig   (FStatInput (..),
                                                   fStatSpecParser)
import           Poseidon.Analysis.Utils          (JackknifeMode (..))
import           Poseidon.Generator.CLI.AdmixPops (AdmixPopsOptions(..), runAdmixPops)

import           Colog                            (logError)
import           Control.Applicative              ((<|>))
import           Control.Exception                (catch)
import           Data.ByteString.Char8            (pack, splitWith)
import           Data.List                        (intercalate)
import qualified Data.Text                        as T
import           Data.Version                     (showVersion)
import qualified Options.Applicative              as OP
import           Paths_poseidon_analysis_hs       (version)
import           Poseidon.PoseidonVersion         (showPoseidonVersion,
                                                   validPoseidonVersions)
import           Poseidon.Utils                   (LogMode (..),
                                                   PoseidonException (..),
                                                   PoseidonLogIO,
                                                   renderPoseidonException,
                                                   usePoseidonLogger)
import           SequenceFormats.Utils            (Chrom (..))
import           System.Exit                      (exitFailure)
import           System.IO                        (hPutStrLn, stderr)
import           Text.Read                        (readEither)
import Poseidon.Generator.Types (IndWithAdmixtureSet)
import Poseidon.Generator.Parsers (readIndWithAdmixtureSetString)
import Poseidon.GenotypeData (GenotypeFormatSpec (..))

data Options =
      CmdFstats FstatsOptions
    | CmdRAS RASOptions
    | CmdAdmixPops AdmixPopsOptions

main :: IO ()
main = do
    hPutStrLn stderr renderVersion
    hPutStrLn stderr ""
    cmdOpts <- OP.customExecParser p optParserInfo
    catch (usePoseidonLogger DefaultLog $ runCmd cmdOpts) handler
    where
        p = OP.prefs OP.showHelpOnEmpty
        handler :: PoseidonException -> IO ()
        handler e = do
            usePoseidonLogger DefaultLog . logError . T.pack $ renderPoseidonException e
            exitFailure

runCmd :: Options -> PoseidonLogIO ()
runCmd o = case o of
    CmdFstats opts -> runFstats opts
    CmdRAS opts    -> runRAS opts
    CmdAdmixPops opts -> runAdmixPops opts

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
    OP.command "ras" rasOptInfo <>
    OP.command "admixpops" admixPopsOptInfo
  where
    fstatsOptInfo = OP.info (OP.helper <*> (CmdFstats <$> fstatsOptParser))
        (OP.progDesc "Compute f-statistics on groups and invidiuals within and across Poseidon packages")
    rasOptInfo = OP.info (OP.helper <*> (CmdRAS <$> rasOptParser))
        (OP.progDesc "Compute RAS statistics on groups and individuals within and across Poseidon packages")
    admixPopsOptInfo = OP.info (OP.helper <*> (CmdAdmixPops <$> admixPopsOptParser))
        (OP.progDesc "Generate individuals with randomized genotype profiles based on admixture proportions")


fstatsOptParser :: OP.Parser FstatsOptions
fstatsOptParser = FstatsOptions <$> parseBasePaths
                                <*> parseJackknife
                                <*> parseExcludeChroms
                                <*> parseFstatInput
                                <*> parseMaxSnps
                                <*> parseTableOutFile

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

parseFstatInput :: OP.Parser [FStatInput]
parseFstatInput = OP.some (parseStatSpecsDirect <|> parseYamlInput <|> parseSimpleText)
  where
    parseStatSpecsDirect = OP.option (FStatInputDirect <$> OP.eitherReader readStatSpecString) (OP.long "stat" <>
        OP.help "Specify a summary statistic to be computed. Can be given multiple times. Possible options are: F4(a, \
            \b, c, d), F3(a, b, c), F2(a, b), PWM(a, b), FST(a, b), Het(a) and some more special options \
            \described at https://poseidon-framework.github.io/#/xerxes?id=fstats-command. Valid entities used in the \
            \statistics are group names as specified in the *.fam, *.ind or *.janno failes, individual names using the \
            \syntax \"<Ind_name>\", so enclosing them in angular brackets, and entire packages like \"*Package1*\" using the \
            \Poseidon package title. You can mix entity types, like in \
            \\"F4(<Ind1>,Group2,*Pac*,<Ind4>)\". Group or individual names are separated by commas, and a comma \
            \can be followed by any number of spaces.")
    parseYamlInput = OP.option (FStatInputYaml <$> OP.str) (OP.long "statConfig" <> OP.help "Specify a yaml file for the Fstatistics and group configurations")
    parseSimpleText = OP.option (FStatInputSimpleText <$> OP.str) (OP.long "statFile" <> OP.help "Specify a file with F-Statistics specified \
    \similarly as specified for option --stat. One line per statistics, and no new-line at the end")
    readStatSpecString s = case runParser fStatSpecParser () "" s of
        Left p  -> Left (show p)
        Right x -> Right x

rasOptParser :: OP.Parser RASOptions
rasOptParser = RASOptions <$>
    parseBasePaths <*>
    parseJackknife <*>
    parseExcludeChroms <*>
    parsePopConfigFile <*>
    parseMinFreq <*>
    parseMaxFreq <*>
    parseMaxMissingness <*>
    parseBlockTableFile <*>
    parseF4tableOutFile <*>
    parseMaxSnps <*>
    parseNoTransitions <*>
    parseBedFile

parsePopConfigFile :: OP.Parser FilePath
parsePopConfigFile = OP.option OP.str (OP.long "popConfigFile" <> OP.help "a file containing the population configuration")

parseMaxFreq :: OP.Parser FreqSpec
parseMaxFreq = parseK <|> parseX <|> parseNone
  where
    parseK = OP.option (FreqK <$> OP.auto) (OP.long "maxAC" <>
        OP.help "define a maximal allele-count cutoff for the RAS statistics. ")
    parseX = OP.option (FreqX <$> OP.auto) (OP.long "maxFreq" <>
        OP.help "define a maximal allele-frequency cutoff for the RAS statistics. ")
    parseNone = OP.flag' FreqNone (OP.long "noMaxFreq" <>
        OP.help "switch off the maximum allele frequency filter. This cam help mimic Outgroup-F3")

parseMinFreq :: OP.Parser FreqSpec
parseMinFreq = parseK <|> parseX <|> parseNone
  where
    parseK = OP.option (FreqK <$> OP.auto) (OP.long "minAC" <>
        OP.help "define a minimal allele-count cutoff for the RAS statistics. ")
    parseX = OP.option (FreqX <$> OP.auto) (OP.long "minFreq" <>
        OP.help "define a minimal allele-frequency cutoff for the RAS statistics. ")
    parseNone = OP.flag' FreqNone (OP.long "noMinFreq" <>
        OP.help "switch off the minimum allele frequency filter")

parseMaxMissingness :: OP.Parser Double
parseMaxMissingness = OP.option OP.auto (OP.long "maxMissingness" <> OP.short 'm' <>
    OP.help "define a maximal missingness for the right populations in the RAS statistics." <>
    OP.value 0.1 <> OP.showDefault)

parseTableOutFile :: OP.Parser (Maybe FilePath)
parseTableOutFile = OP.option (Just <$> OP.str) (OP.long "tableOutFile" <> OP.short 'f' <>
    OP.help "a file to which results are written as tab-separated file" <> OP.value Nothing)

parseF4tableOutFile :: OP.Parser (Maybe FilePath)
parseF4tableOutFile = OP.option (Just <$> OP.str) (OP.long "f4TableOutFile" <>
    OP.help "a file to which F4 computations are written as tab-separated file" <> OP.value Nothing)

parseBlockTableFile :: OP.Parser (Maybe FilePath)
parseBlockTableFile = OP.option (Just <$> OP.str) (OP.long "blockTableFile" <>
    OP.help "a file to which the per-Block results are written as tab-separated file" <> OP.value Nothing)

-- parseFstatConfig :: OP.Parser (Maybe FilePath)
-- parseFstatConfig = OP.option (Just <$> OP.str) (OP.long "config" <> OP.help "Config file in Dhall" <> OP.value Nothing)

parseMaxSnps :: OP.Parser (Maybe Int)
parseMaxSnps = OP.option (Just <$> OP.auto) (OP.long "maxSnps" <>
    OP.help "Stop after a maximum nr of snps has been processed. Useful for short test runs" <>
    OP.value Nothing <> OP.hidden)

parseNoTransitions :: OP.Parser Bool
parseNoTransitions = OP.switch (OP.long "noTransitions" <> OP.help "Skip transition SNPs and use only transversions")

parseBedFile :: OP.Parser (Maybe FilePath)
parseBedFile = OP.option (Just <$> OP.str) (OP.long "bedFile" <> OP.help "An optional bed file that gives sites to be \
    \included in the analysis." <> OP.value Nothing)

admixPopsOptParser :: OP.Parser AdmixPopsOptions
admixPopsOptParser = AdmixPopsOptions <$> parseBasePaths
                                      <*> parseIndWithAdmixtureSetDirect
                                      <*> parseIndWithAdmixtureSetFromFile
                                      <*> parseMarginalizeMissing
                                      <*> parseOutGenotypeFormat
                                      <*> parseOutPath

parseIndWithAdmixtureSetDirect :: OP.Parser [IndWithAdmixtureSet]
parseIndWithAdmixtureSetDirect = OP.option (OP.eitherReader readIndWithAdmixtureSetString) (
    OP.long "admixString" <>
    OP.short 'a' <>
    OP.value [] <>
    OP.help "Population setup of interest: Each setup is a string of the form \
            \\"[id:group](population1=10+population2=30+...)\". Multiple setups can be listed separated by ;. \
            \id and group are simple strings. \
            \The population fractions must be simple integers and sum to 100."
    )

parseIndWithAdmixtureSetFromFile :: OP.Parser (Maybe FilePath)
parseIndWithAdmixtureSetFromFile = OP.option (Just <$> OP.str) (OP.long "admixFile" <>
    OP.value Nothing <>
    OP.help "A file with a list of spatiotemporal positions. \
            \Works just as -p, but multiple values can be given separated by newline. \
            \-a and --admixFile can be combined."
    )

parseMarginalizeMissing :: OP.Parser Bool
parseMarginalizeMissing = OP.switch (
    OP.long "marginalizeMissing" <> 
    OP.help "ignore missing SNPs in the per-population genotype frequency calculation \
            \(except all individuals have missing information for a given SNP)"
    )

parseOutGenotypeFormat :: OP.Parser GenotypeFormatSpec
parseOutGenotypeFormat = OP.option (OP.eitherReader readGenotypeFormat) (
    OP.long "outFormat" <>
    OP.help "The format of the output genotype data: EIGENSTRAT or PLINK" <>
    OP.value GenotypeFormatEigenstrat
    )
    where
    readGenotypeFormat :: String -> Either String GenotypeFormatSpec
    readGenotypeFormat s = case s of
        "EIGENSTRAT" -> Right GenotypeFormatEigenstrat
        "PLINK"      -> Right GenotypeFormatPlink
        _            -> Left "must be EIGENSTRAT or PLINK"

parseOutPath :: OP.Parser FilePath
parseOutPath = OP.strOption (
    OP.long "outPath" <>
    OP.short 'o' <>
    OP.help "The output directory path"
    )
