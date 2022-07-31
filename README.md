[![GitHub Workflow Status](https://img.shields.io/github/workflow/status/poseidon-framework/poseidon-analysis-hs/CI)](https://github.com/poseidon-framework/poseidon-analysis-hs/actions?query=workflow%3ACI)
[![Coverage Status](https://img.shields.io/codecov/c/github/poseidon-framework/poseidon-analysis-hs/main.svg)](https://codecov.io/github/poseidon-framework/poseidon-analysis-hs?branch=main)
[![GitHub release (latest by date including pre-releases)](https://img.shields.io/github/v/release/poseidon-framework/poseidon-analysis-hs?include_prereleases) ![GitHub all releases](https://img.shields.io/github/downloads/poseidon-framework/poseidon-analysis-hs/total)](https://github.com/poseidon-framework/poseidon-analysis-hs/releases)
[![Install with Bioconda](https://anaconda.org/bioconda/poseidon-xerxes/badges/installer/conda.svg)](https://anaconda.org/bioconda/poseidon-xerxes)

# poseidon-analysis-hs
Tools to analyse genotype data in the Poseidon package format. The main executable within this package is called `xerxes`.

**Detailed user documentation can be found on our [website](https://poseidon-framework.github.io/#/xerxes).**

***

## For (Haskell) developers

To install the development version of poseidon-analysis-hs/xerxes you can follow these steps:

1. Install the Haskell build tool [Stack](https://docs.haskellstack.org/en/stable/README/)
2. Clone this repository
3. Execute `stack install` inside the repository to build the tool and copy the executables to `~/.local/bin` (which you may want to add to your path). This will install the compiler and all dependencies into folders that won't interfere with any installation you might already have.
4. To run the tests, execute `stack test` inside the repository to build and run tests.
