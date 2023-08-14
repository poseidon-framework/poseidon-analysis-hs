[![CI](https://github.com/poseidon-framework/poseidon-analysis-hs/actions/workflows/main.yml/badge.svg)](https://github.com/poseidon-framework/poseidon-analysis-hs/actions/workflows/main.yml)
[![GitHub release (latest by date including pre-releases)](https://img.shields.io/github/v/release/poseidon-framework/poseidon-analysis-hs?include_prereleases) ![GitHub all releases](https://img.shields.io/github/downloads/poseidon-framework/poseidon-analysis-hs/total)](https://github.com/poseidon-framework/poseidon-analysis-hs/releases)
[![Install with Bioconda](https://anaconda.org/bioconda/poseidon-xerxes/badges/version.svg)](https://anaconda.org/bioconda/poseidon-xerxes)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/poseidon-xerxes/badges/downloads.svg)](https://anaconda.org/bioconda/poseidon-xerxes)

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

### Preparing a new stable release

The Github Actions script in `.github/workflows/release.yml` registers a new draft release and automatically builds and uploads xerxes binaries when a new Git tag with the prefix `v*` is pushed. 

```bash
# locally register a new tag (e.g. 0.3.1)
git tag -a v0.3.1 -m "see CHANGELOG.md"
# push tag
git push origin v0.3.1
```

In case of a failing build delete the tag and the release draft on Github and then delete the tag locally with

```bash
git tag -d v0.3.1
```

before rerunning the procedure above.