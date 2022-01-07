module Utils (JackknifeMode(..)) where

-- | A datatype representing the two options for how to run the Block-Jackknife
data JackknifeMode = JackknifePerN Int -- ^ Run the Jackknife in blocks of N consecutive SNPs
    | JackknifePerChromosome -- ^ Run the Jackknife in blocks defined as entire chromosomes.
  deriving (Show)
