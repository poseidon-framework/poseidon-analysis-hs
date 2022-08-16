module Poseidon.Generator.SampleGeno where

import           Poseidon.Generator.Types

import           Control.Monad.Random       (RandomGen, evalRand, fromList,
                                             getStdGen, newStdGen)
import           Data.List
import           Data.Ratio                 ((%))
import qualified Data.Vector                as V
import           Pipes
import           Pipes.Safe
import           SequenceFormats.Eigenstrat

-- admixpops

samplePerChunk :: 
       [IndAdmixpops]
    -> Producer (EigenstratSnpEntry, GenoLine) (SafeT IO) r
    -> Producer (EigenstratSnpEntry, GenoLine) (SafeT IO) r
samplePerChunk inds prod = do
    gen <- liftIO getStdGen
    let huhu = head $ map _popSet inds
    let bobo = zip (map _popName huhu) (map _popFrac huhu)
    let schuhu = sampleWeightedList gen bobo
    liftIO $ putStrLn $ show schuhu
    _ <- liftIO newStdGen
    for prod (handleEntry inds)
  where
    handleEntry :: [IndAdmixpops] -> (EigenstratSnpEntry, GenoLine) -> Producer (EigenstratSnpEntry, GenoLine) (SafeT IO) ()
    handleEntry inds_ x = do
        yield x

getFracPerInd :: IndAdmixpops -> [(Int, Rational)]
getFracPerInd x = 
    let popsSets = _popSet x
        popsFracs = map _popFrac popsSets
        popsInds = map _popInds popsSets
        lengthPopsInds = map length popsInds
        fracPerIndPerPop = map (\(x,y) -> x / (fromIntegral y)) $ zip popsFracs lengthPopsInds
        popFracsPerInd = map (\(x, y) -> replicate y x) $ zip fracPerIndPerPop lengthPopsInds
    in zip (concat popsInds) (concat popFracsPerInd)

samplePerSNP ::
       Bool
    -> [IndAdmixpops]
    -> (EigenstratSnpEntry, GenoLine)
    -> SafeT IO (EigenstratSnpEntry, GenoLine)
samplePerSNP marginalizeMissing inds (snpEntry, genoLine) = do
    entries <- mapM (\x -> samplePerSNPForOneOutInd marginalizeMissing x genoLine) inds
    return (snpEntry, V.fromList entries)

samplePerSNPForOneOutInd :: Bool -> IndAdmixpops -> GenoLine -> SafeT IO GenoEntry
samplePerSNPForOneOutInd marginalizeMissing ind genoLine = do
    gen <- liftIO getStdGen
    let indIDsAndFracs = map (\(PopAdmixpops _ frac_ inds_) -> (inds_, frac_)) $ _popSet ind
    let sampledGenotypesPerPop = map (\(x,y) -> (sampleWeightedList gen $ getGenotypeFrequency marginalizeMissing x genoLine, y)) indIDsAndFracs
    let sampledGenotypeAcrossPops = sampleWeightedList gen sampledGenotypesPerPop
    _ <- liftIO newStdGen
    return sampledGenotypeAcrossPops

getGenotypeFrequency :: Bool -> [Int] -> GenoLine -> [(GenoEntry, Rational)]
getGenotypeFrequency marginalizeMissing individualIndices genoLine =
    let relevantGenoEntries = [genoLine V.! i | i <- individualIndices]
    in  if all (Missing ==) relevantGenoEntries
        then [(Missing, 1)]
        else if marginalizeMissing
             then calcFractions $ filter (Missing /=) relevantGenoEntries
             else calcFractions relevantGenoEntries

calcFractions :: [GenoEntry] -> [(GenoEntry, Rational)]
calcFractions xs =
    let ls = toInteger $ length xs
    in map (\x -> (head x, toInteger (length x) % ls)) $ group $ sortBy compareGenoEntry xs

compareGenoEntry :: GenoEntry -> GenoEntry -> Ordering
compareGenoEntry Het x = case x of
    Het -> EQ
    _   -> LT
compareGenoEntry HomRef x = case x of
    Het    -> GT
    HomRef -> EQ
    _      -> LT
compareGenoEntry HomAlt x = case x of
    Het    -> GT
    HomRef -> GT
    HomAlt -> EQ
    _      -> LT
compareGenoEntry Missing x = case x of
    Missing -> EQ
    _       -> GT

sampleWeightedList :: RandomGen g => g -> [(a, Rational)] -> a
sampleWeightedList gen weights = head $ evalRand m gen
    where m = sequence . repeat . fromList $ weights

-- spacetime

sampleGenoForMultiplePOIs :: [([Int], [Rational])] -> (EigenstratSnpEntry, GenoLine) -> SafeT IO (EigenstratSnpEntry, GenoLine)
sampleGenoForMultiplePOIs infoForIndividualPOIs (snpEntry, genoLine) = do
    entries <- mapM (\(x,y) -> sampleGenoForOnePOI x y genoLine) infoForIndividualPOIs
    return (snpEntry, V.fromList entries)

sampleGenoForOnePOI :: [Int] -> [Rational] -> GenoLine -> SafeT IO GenoEntry
sampleGenoForOnePOI individualIndices weights genoLine = do
    let relevantGenoEntries = [genoLine V.! i | i <- individualIndices]
        -- count occurrence of GenoEntries
        genoEntryIndices = getGenoIndices relevantGenoEntries
        -- sum distance-based weight for each GenoEntry
        weightsPerGenoEntry = sumWeights genoEntryIndices weights
    -- sample GenoEntry based on weight
    gen <- liftIO getStdGen
    -- liftIO $ hPutStrLn stderr (show gen)
    let selectedGenoEntry = sampleWeightedList gen weightsPerGenoEntry
    _ <- liftIO newStdGen
    -- return
    return selectedGenoEntry

getGenoIndices :: Eq a => [a] -> [(a, [Int])]
getGenoIndices xs =
    let unique = nub xs
        indices = map (`elemIndices` xs) unique
    in  zip unique indices

sumWeights :: Num b => [(a, [Int])] -> [b] -> [(a, b)]
sumWeights xs weights = map (\(x, ys) -> (x, sum $ subset ys weights)) xs
    where
        subset :: [Int] -> [a] -> [a]
        subset indices ys = [ys !! i | i <- indices]
