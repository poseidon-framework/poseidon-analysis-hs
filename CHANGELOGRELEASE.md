# V0.3.0.0: Added a first implementation of admixpops

For this release we integrated the very experimental artificial ancestry generation subcommand `admixpops` into xerxes (formerly a component of a tool called [paagen](https://github.com/nevrome/paagen)). We also added some innovations from trident for logging, error handling and code maintenance.

`admixpops` allows to generate individuals with randomized genotype profiles based on arbitrary input populations. Running it for example with the following `--admixString` would yield an artificial individual `Ind1` with 10% ancestry from the reference population Mbuti, 20% from Han, and 70% from French: `"[Ind1:Pop1](Mbuti=10+Han=20+French=70)"`. The implemented algorithm is (trivial) SNP-wise weighted random-sampling. `admixpops` is therefore only an experimental feature and the current implementation will likely change in the future.

In the process of adding it to xerxes, we also upgraded the latter with some of the (recent) improvements in trident. From a user-perspective the only relevant changes are the new general options `--logMode` and `--errLength`, which give you more control over the logging output and error messages.

# V0.2.1.0: Optimized fstats with allele-frequency lookups

This release internally optimizes the core loop by creating lookup tables for allele frequencies, so that they can be reused for the various statistics. 

Code has been refactored and unit tests introduced.

Logging is partly implemented with Co-Log now, but not fully yet.

# V0.2.0.0: New config file input, ascertainment and adhoc-groups Latest

This is a big release with three new features in Xerxes fstats:

1.) A new input option for F-Statistics using a YAML configuration file, supporting automatic looping through multiple statistics.
2.) The option of entering SNP ascertainment, using allele frequency bounds in a reference population
3.) The option of defining adhoc-groups, for example to exclude outlier individuals from populations without changing the genotype data files.

See documentation on https://poseidon-framework.github.io/#/xerxes
