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
