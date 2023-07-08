<!-- vim-markdown-toc GFM -->

* [CONSULT-II](#consult-ii)
    * [Overview](#overview)
    * [Getting started](#getting-started)
        * [System requirements](#system-requirements)
            * [Memory and disk space](#memory-and-disk-space)
            * [Dependencies](#dependencies)
        * [Installation](#installation)
        * [Testing](#testing)
    * [Guide](#guide)
        * [Constructing a reference library](#constructing-a-reference-library)
            * [Construction of the hash table](#construction-of-the-hash-table)
            * [Adding taxonomic LCA information](#adding-taxonomic-lca-information)
        * [Taxonomic identification](#taxonomic-identification)
            * [Classification of reads](#classification-of-reads)
            * [Abundance profiling](#abundance-profiling)
            * [Contamination removal](#contamination-removal)
        * [Input and output formats](#input-and-output-formats)
        * [Description of CONSULT-II arguments and usage](#description-of-consult-ii-arguments-and-usage)
        * [Public libraries](#public-libraries)
        * [Useful tips](#useful-tips)

<!-- vim-markdown-toc -->
# CONSULT-II
## Overview
CONSULT-II is a tool for taxonomic identification, and successor of CONSULT.
CONSULT-II extends CONSULT by enabling taxonomic classification and abundance profiling, in addition to contamination removal.

Relying on locality-sensitive hashing, CONSULT-II extracts *k*-mers from a query set and tests whether they fall within a user-specified hamming distance of *k*-mers in the reference dataset.
Using this invaluable information, it can derive taxon predictions for given query reads, or a profile estimation for a query sample.
Furthermore, for a given query read, CONSULT-II can report *k*-mer matches and corresponding Hamming distances, and optionally *probabilistic* taxonomic least common ancestor (LCA) of matched reference *k*-mers.
It supports the inclusion of approximately billions *k*-mers in its reference library, accommodating datasets with tens of thousands of microbial species.
CONSULT-II is efficiently parallelized, and can handle very large datasets.

Despite being memory hungry, our careful benchmarking shows that CONSULT-II outperforms popular *k*-mer based tools such as Kraken-2 and CLARK.
We provide reference libraries, so you can download them and jump into taxonomic identification of your queries.

## Getting started
### System requirements
#### Memory and disk space
Exact memory footprint depends on a number of *k*-mers in a reference set, and the parameters used.
CONSULT-II requires enough free memory to hold the entire reference library in RAM, and this is needed during library construction and searching query *k*mers.
During, classification and abundance profiling, CONSULT-II reads matching information from the disk.
Using the default values (or heuristic), some examples of approximately required memory (conservative upper bounds) with respect to the number of *k*-mers in a reference set can be listed as below;

* $2^{28}$ *k*-mers $\rightarrow$ $<5$ GB,

* $2^{30}$ *k*-mers $\rightarrow$ $<20$ GB,

* $2^{32}$ *k*-mers $\rightarrow$ $<80$ GB.

We note that during library construction the user will need slightly more RAM than the given values to accommodate intermediary processes (about additional 10\%).
Once the database is built, these values should be sufficient.
For instance, the main reference library, with more than 10k species and about 8 billion *k*-mers leads to memory usage varying between 140GB and 150GB.
In addition to reference library and metadata, CONSULT-II also needs to store *k*-mer match information per read to disk, which then be read during classification and profiling.
This shouldn't exceed input data file sizes (usually less than 10\% of input size).

#### Dependencies
CONSULT-II is a command-line tool implemented in C++11 with some x86 assembly code.
Many steps of CONSULT-II are embracingly parallelized, such as database reading, query search, classification, profiling, using [OpenMP](https://www.openmp.org).
Compilation requires g++ that supports C++11 (required).
For our tests, we have compiled versions 4.8.5 and 7.2.0, both of which work.
External tools such as [Jellyfish](http://www.genome.umd.edu/jellyfish.html) might be useful for database construction.

### Installation
1. Download using one of two approaches:
    - You can obtain the [zip file](https://github.com/noraracht/CONSULT/archive/main.zip) and extract the contents to a folder of your choice.
    Then, proceed to compilation.
    - Alternatively, you can clone the [github repository](https://github.com/noraracht/CONSULT.git) and continue with compilation.

2. To compile, go to the directory where core programs for map construction and query search are located and run the below commands.
    * You can use ``make`` to compile CONSULT.
    ```bash
    make all # for all components of CONSULT
    # OR
    make minimize # for the minimization script
    make map # for consult_map to construct a library
    make search # for consult_search to make queries
	make classify # for consult_classify to perform read-level classification from match info
	make profile # for consult_profile to perform abundance profiling from match info
    ```
    * Alternatively, you can run ``g++`` directly.
    ```bash
    g++ minimize.cpp -std=c++11 -o minimize # for the minimization script
    g++ consult_map.cpp -std=c++11 -O3 -o consult_map  # for consult_map to construct a library
    g++ consult_search.cpp -std=c++11 -fopenmp -O3 -o consult_search # for consult_search to make queries
	g++ consult_classify.cpp -std=c++11 -fopenmp -O3 -o consult_classify # for consult_classify to perform read-level classification from match info
	g++ consult_profile.cpp -std=c++11 -fopenmp -O3 -o consult_profile # for consult_profile to perform abundance profiling from match info
    ```

### Testing
To test your installation, and gain a better insight about input/output files and CONSULT-II's workflow, see the directory `example`.
To test all functionality, run `make all` in `example`.
As a results, following files should be generated:

- `G000307305-k32C_minimized.fa`: minimized *k*-mers, in FASTA format,

- `G000307305_nbr_mapping`: a directory to store reference library files, such as *k*-mer encodings, taxonomic LCAs, metadata,

- `classified-seq_query000`: list of reads with at least one (default `c` value) *k*-mer match,

- `unclassified-seq_query000`: list of reads without any *k*-mer match,

- `kmer-distances_query000`: number of *k*-mer matches across different Hamming distances for each read,

- `match-info_query000`: list of matched *k*-mers with taxonomic LCA, this file is used by `consult_classify` and `consult_profile`,

- `classification_query000`: taxon predictions for each read,

- `profile_query000-*`: an abundance profile for each rank.


## Guide

### Constructing a reference library
#### Construction of the hash table

#### Adding taxonomic LCA information

### Taxonomic identification
#### Classification of reads

#### Abundance profiling

#### Contamination removal

### Input and output formats

### Description of CONSULT-II arguments and usage

### Public libraries

### Useful tips
