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
        * [Preprocessing](#preprocessing)
        * [Constructing a reference library](#constructing-a-reference-library)
            * [Construction of the hash table](#construction-of-the-hash-table)
            * [Adding taxonomic LCA information to reference library](#adding-taxonomic-lca-information-to-reference-library)
        * [Taxonomic identification](#taxonomic-identification)
            * [Searching a query against a reference library](#searching-a-query-against-a-reference-library)
            * [Classification of reads](#classification-of-reads)
            * [Abundance profiling](#abundance-profiling)
            * [Contamination removal](#contamination-removal)
        * [Input and output formats](#input-and-output-formats)
            * [FASTA file for library construction](#fasta-file-for-library-construction)
            * [Query files for taxonomic identification](#query-files-for-taxonomic-identification)
            * [Taxonomy lookup and filename map](#taxonomy-lookup-and-filename-map)
            * [Report for number of matched *k*-mers and their Hamming distances](#report-for-number-of-matched-k-mers-and-their-hamming-distances)
            * [Information for match distances and corresponding taxonomic LCAs](#information-for-match-distances-and-corresponding-taxonomic-lcas)
        * [Description of CONSULT-II arguments and usage](#description-of-consult-ii-arguments-and-usage)
            * [``minimize``](#minimize)
            * [``consult_map``](#consult_map)
            * [``consult_search``](#consult_search)
        * [Public libraries](#public-libraries)

<!-- vim-markdown-toc -->
# CONSULT-II
## Overview
CONSULT-II is a tool for taxonomic identification, and successor of CONSULT.
CONSULT-II extends CONSULT-II by enabling taxonomic classification and abundance profiling, in addition to contamination removal.

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
### Preprocessing

### Constructing a reference library
#### Construction of the hash table
To construct a standard reference database, go to the place where scripts were compiled and use the following command to run ``consult_map``.
```bash
./consult_map -i $INPUT_FASTA_FILE -o $DB_NAME
# For example:
./consult_map -p 3 -t 2 -i k32C_af_mininimization.fa -o G000307305_nbr_map
```
``-t`` is tag size in bits, and determines the number of partitions ( $2^t$ ).
``-p`` is the Hamming distance threshold for a match.
`$INPUT_FASTA_FILE` is supposed to be the path of input FASTA file, formatted as shown in the corresponding [section](#fasta-file-for-library-construction).
Replace `$DB_NAME` above with your preferred database name, which also will be the name of the database directory.
You can also give a path to some other valid directory other than current working directory: `$/path/to/$DB_NAME`.
If this working directory already contains a database with the same name, the software will throw an exception.
This feature is included to prevent existing databases from being overwritten.
In the database directory, many non human-readable binary files will be stored, such as metadata, encoding arrays, tag array (in chunks).

#### Adding taxonomic LCA information to reference library
After constructing the hash table, and populating it with reference *k*-mers, we need to add taxonomic LCA information to each *k*-mer for abundance profiling and taxonomic classification.
For contamination removal, this step is needed since matching to *any* reference sequence is sufficient in that case.
Two consecutive commands must be run, and filename and taxonomy look-up files must be given as argunments:
```bash
./consult_search -q $REFERENCE_GENOMES -i /path/to/$DB_NAME -o . --taxonomy-lookup-path $TAXONOMY_LOOKUP --filename-map-path $FILENAME_MAP --init-ID
./consult_search -q $REFERENCE_GENOMES -i /path/to/$DB_NAME -o . --taxonomy-lookup-path $TAXONOMY_LOOKUP --filename-map-path $FILENAME_MAP --update-ID
```
As the name suggests, `$DB_NAME` is the name of the library (also the name of the directory), and `$REFERENCE_GENOMES` is the path to the directory in which all the reference genomes are stored in FASTQ format separately with appropriate file names.
Note that, filenames are going to be used to map each genome file to a species.
See the [taxonomy lookup and filename map section](#taxonomy-lookup-and-filename-map) for a detailed description and how to generate them.
The command with flag `--initID`, initializes each label associated with *k*-mers, and counts how many distinct genomes have each *k*-mer.
Then, `update-ID` computes and stores the probabilistic LCA labels for each *k*-mer.
Also, make sure to set `--thread-count` to a value as high as possible as this step might be pretty slow but efficiently parallelized.

### Taxonomic identification
#### Searching a query against a reference library
To query a set of sequences against a reference, go to the directory where binaries are and execute the CONSULT-II command:
```bash
 ./consult_search -i $DB_NAME -q $QUERY_PATH -o $OUTPUT_DIRECTORY
# For example:
./consult_search -q query000.fq -i G000307305_nbr_mapping -o . -c 1
```
The files containing query sequences to be searched should be located in ``$QUERY_PATH`` and be in a FASTQ format (one uncompressed ``.fq``/``.fastq`` file per each sample).
See [section](#query-files-for-taxonomic-identification) for a detailed description of the query path.
This step is required for all taxonomic identification tasks: classification, profiling and contamination removal.
But, as described below, different flags have to be used for different tasks, above `consult_search` command is not useful by itself.

#### Classification of reads

#### Abundance profiling

#### Contamination removal
For contamination removal, there is no need to run an additional command: `consult_search` also is able to generate a file that contains the **classified reads** and **unclassified read**
To make CONSULT-II behave this way, give the ``--classified-out`` and ``--unclassified-out`` flags, correspondingly. The output file name will be prefixed with *"classified-seq_"* and *"classified-seq_"*  .
Files are stored in the directory given in ``$OUTPUT_DIRECTORY``, and the default is where software is run.
Every sample retains its original file name prefixed with *"classified-seq_"* or *"unclassified-seq_"*.
```bash
 ./consult_search -i $DB_NAME -q $QUERY_PATH -o $OUTPUT_DIRECTORY --classified-out --unclassified-out --save-distances
# For example:
./consult_search -q query000.fq -i G000307305_nbr_mapping -o . -c 1 --classified-out --unclassified-out --save-distances
```
Another output that CONSULT-II can generate is a tab-separated file, in which, each row is a read and the column values are the total number of *k*-mers in that read which matches with some *k*-mer in the reference library with Hamming distance $d$.
Each column corresponds to a distance value $d$.
CONSULT-II outputs this file, named file name prefixed with *"kmer-distances_*" if the ``--save-distances.`` flag is given.
The maximum distance value included in this file is determined by the $p$ value (Hamming distance threshold for a match) as a default ( $\lceil 1.5p \rceil$ ), but can be set to some other value with argument ``--maximum-distance``.

Considering the example, unclassified reads would be stored in ``unclassified-seq_G000307305.fq``.
If the flag ``--classified-out`` is given, classified reads would be stored in ``classified-seq_G000307305.fq``.
If the flag ``--save-distances`` is given, distance values would be stored in ``kmer-distances_G000307305.fq``.

### Input and output formats
#### FASTA file for library construction
Specifically, CONSULT-II is designed to accept [Jellyfish](http://www.genome.umd.edu/jellyfish.html) output files that represent a list of 32 bp *k*-mers associated with their counts.
We tested with [Jellyfish](http://www.genome.umd.edu/jellyfish.html) 2.3.0.
See the [preprocessing section](#preprocessing) for details on how to generate the input file.
Note that CONSULT-II does not use the count values and the only relevant information is the sequence itself.
Jellyfish output is pseudo-randomly ordered, and thus, further randomization is not needed.
The sequences should not be repeated.

Example FASTA:
```
> FASTA sequence 1 label
AGACGAGCTTCTTCATCTAAAATGAATTCTCC
> FASTA sequence 2 label
CCAGCTGCATTGATGGTGGGAGTTGGTAAAGG
> FASTA sequence 3 label
GGACCTTGATTTTGACAAGATAGTCGGTAGAC
> FASTA sequence 4 label
ACCACATTTTATACATCGTAAGACAAGCGGCT
```

#### Query files for taxonomic identification
The path ``$QUERY_PATH`` can be a directory containing ``.fastq`` files or a ``.fastq``.
If it is a directory, each query file in the directory will be queried against the library, and separate outputs will be generated for each.
FASTA format is not supported at the moment.
Note, if you need to query FASTA files you can convert ``.fasta``/``.fa`` to ``.fastq``/``.fq`` using [seqtk](https://github.com/lh3/seqtk) ``seqtk seq -F CHAR`` command which attaches fake quality scores to the sequences.
Quality factors are not being utilized by CONSULT-II but FASTQ labels will be used to identify the sequences in the output file.

Example FASTQ:
```
@ FASTQ sequence 1 label
CATCGAGCAGCTATGCAGCTACGAGT
+
-$#'%-#.&)%#)"".)--'*()$)%
@ FASTQ sequence 2 label
TACTGCTGATATTCAGCTCACACC
+
,*#%+#&*$-#,''+*)'&.,).,
```

#### Taxonomy lookup and filename map

#### Report for number of matched *k*-mers and their Hamming distances

#### Information for match distances and corresponding taxonomic LCAs

### Description of CONSULT-II arguments and usage
#### ``minimize``

- ``-i`` or ``--input-fasta-file``: input ``.fasta`` file containing canonical *k*-mers, default length is 35.
- ``-o`` or ``--output-fasta-file``: ``.fasta`` file to output minimizers of the given canonical *k*-mers, default length is 32.

#### ``consult_map``

- ``-i`` or ``--input-fasta-file``: input ``.fasta`` file to construct library.
- ``-o`` or ``--output-library-dir``: output path to the directory that will constitute the CONSULT-II library.
- ``-h`` or ``--number-of-positions``: number of randomly positioned bits to compute LSH.
- ``-t`` or ``--tag-size``: number of bits to be used as tag.
- ``-l`` or ``--number-of-tables``: number of tables, i.e., number of hash functions.
- ``-b`` or ``--column-per-tag``: number of columns per each tag partition, i.e., number of *k*-mers each encoding can map to.

#### ``consult_search``
- ``-i`` or ``input-library-dir``: directory of the CONSULT-II library that will be used as the reference database.
- ``-o`` or ``output-result-dir``: directory in which all output files (classified and unclassified reads, matching information) will be saved.
- ``-q`` or ``--query-path``: the path to the query file, or to the directory containing query files (or reference genomes to add taxonomic information to input library).
- ``-c`` or ``--number-of-matches``: the minimum number of matched *k*-mers that is required to call sequencing read classified.
For instance, if at least one *k*-mer match is enough to classify a read (default setting mentioned in a paper), ``-c`` should be set to 1 in the software.
If at least two *k*-mer matches are required to call the entire read a match, ``-c`` should be set to 2.
Default value is 1.
- ``--thread-count``: number of threads to be used.
- ``--unclassified-out``: to output reads that are unclassified in a file with a name query file name prefixed with *"unclassified-seq_"*.
This is given by default.
- ``--classified-out``: to output reads that are classified in a file with a name query file name prefixed with *"classified-seq_"*.
- ``--save-distances``: to save number of matched *k*-mers in a tab-separated file where columns are the distances of corresponding counts.
File name is query file name prefixed with *"kmer-distances_"*. See the corresponding [section](#report-for-number-of-matched-k-mers-and-their-hamming-distances).
- ``--maximum-distance``: maximum distance to be included as a column in the file containing *k*-mer match counts with respect varying Hamming distance values.
Note that, when the ``--maximum-distance`` value is too large compared to $p$, CONSULT-II does not necessarily aim to find such *k*-mers to compute column values.
This is because the library size and the corresponding theoretical guarantees only consider $p$ value.
So, if one would like to use a large ``--maximum-distance`` value, $p$ value should be increased proportionally.
- ``init-ID``: to initialize taxonomic LCA labels (set them to the root of the tree), and to count distinct genomes in which each *k*-mer appears.
- ``update-ID``: to compute and store probabilistic LCA labels for each *k*-mer in the reference library.
It is required to run `consult_search` with this flag (and `--init-ID` before classification and profiling.
- ``--taxonomy-lookup-path``: path of the taxonomic LCA lookup table in a human-readable format.
See the corresponding [section](#taxonomy-lookup-and-filename-map) for details.
- ``--filename-map-path``: path of the filename/genome to species ID map in a human-readable format.
See the corresponding [section](#taxonomy-lookup-and-filename-map) for details.
- ``--save-matches``: to output detailed *k*-mer match information of query sequences consisting of Hamming distances and taxonomic LCAs.
See the corresponding [section](#information-for-match-distances-and-corresponding-taxonomic-lcas) for details.

### Public libraries
- Soon.
