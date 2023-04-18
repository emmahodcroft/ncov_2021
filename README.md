[![GitHub release (latest by date)](https://img.shields.io/github/v/release/nextstrain/ncov)](https://github.com/nextstrain/ncov/releases)
[![See recent changes](https://img.shields.io/badge/changelog-See%20recent%20changes-blue)](https://docs.nextstrain.org/projects/ncov/en/latest/reference/change_log.html)

# About

This branch is a frozen version of the `ncov_2021` master branch used for the Alpha & Delta cluster introductions to Switzerland analysis.
The `ncov_2021` master branch is the copy of the Nextstrain `ncov` pipeline used by Emma Hodcroft to create cluster-focused builds, primarily for CoVariants.

Almost all code used for the Alpha & Delta introduction analysis (to generate the phylogenies) is the same, with a couple of exceptions:
- `scripts/reduceMetadata.py` was used to reduce the `metadata.tsv` which is one of the starting files of analysis to only have sequences with dates prior to 31 Aug 2021, since our analysis focuses only prior to this period. This reduces the size of the file and simply makes the rest of the pipeline more efficient
- `scripts/checkDatesOfMetadata.py` to check the above script works properly, this script outputs the min & max dates of the before & after metadata files

- `scripts/priorities.py` Main change. This adds an `exclusion_penalty` (of -10000) to the priority scores of half of the possible context sequences, essentially downsampling the global context sequences by half. This simulates a scenario where the rest of the world (exlcuding Switerland) had done half as much sequencing.

