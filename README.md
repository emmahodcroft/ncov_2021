[![DOI](https://zenodo.org/badge/339133771.svg)](https://zenodo.org/badge/latestdoi/339133771)
<- Release version for manuscript

# About

This branch is a frozen version of the `ncov_2021` master branch used for the Alpha & Delta cluster introductions to Switzerland analysis.
The `ncov_2021` master branch is the copy of the Nextstrain `ncov` pipeline used by Emma Hodcroft to create cluster-focused builds, primarily for CoVariants.

Almost all code used for the Alpha & Delta introduction analysis (to generate the phylogenies) is the same, with a couple of exceptions specifically for the reruns with reduced background sequences:
- `scripts/reduceMetadata.py` was used to reduce the `metadata.tsv` which is one of the starting files of analysis to only have sequences with dates prior to 31 Aug 2021, since our analysis focuses only prior to this period. This reduces the size of the file and simply makes the rest of the pipeline more efficient
- `scripts/checkDatesOfMetadata.py` to check the above script works properly, this script outputs the min & max dates of the before & after metadata files

- `scripts/priorities.py` Main change. This adds an `exclusion_penalty` (of -10000) to the priority scores of half of the possible context sequences, essentially downsampling the global context sequences by half. This simulates a scenario where the rest of the world (exlcuding Switerland) had done half as much sequencing.

To run the analysis without the reduced background sequenes (as in the main Alpha/Delta analyses of the paper), simply revert the above change in `scripts/priorities.py`, which adds the -10000 penalty, and start with a full metadata file (not a reduced one as above) and the rest of the code remains the same. 
