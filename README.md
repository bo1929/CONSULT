# Announcement
We recommend using [KRANK](https://github.com/bo1929/KRANK).
KRANK uses the **same algorithm** as CONSULT-II for sequence classification and abundance profiling.
It is faster, much more user-friendly, and is simply better software.
It wraps many redundancies of CONSULT-II's procedures into single commands.
KRANK and CONSULT-II only differ in terms of *k*-mers selected to keep in the library.
KRANK develops some heuristics to select a **better** *k*-mer subset and hence requires less memory to maintain the same performance as CONSULT-II.
It also is very easy to emulate CONSULT-II using the `--fast-mode` option of KRANK, which is much faster!
You may want to cite both if you use KRANK as it borrows its classification algorithm from CONSULT-II.

# CONSULT-II
*Ali Osman Berk Şapcı, Eleonora Rachtman, Siavash Mirarab, CONSULT-II: Accurate taxonomic identification and profiling using locality-sensitive hashing, Bioinformatics, 2024;, btae150, https://doi.org/10.1093/bioinformatics/btae150*

- [CONSULT-II](#consult-ii)
  * [Overview](#overview)
  * [Getting started](#getting-started)
    + [System requirements](#system-requirements)
      - [Memory and disk space](#memory-and-disk-space)
      - [Dependencies](#dependencies)
    + [Installation](#installation)
    + [Testing](#testing)
  * [Guide](#guide)
    + [Constructing a reference library](#constructing-a-reference-library)
      - [Preprocessing](#preprocessing)
      - [Construction of the hash table](#construction-of-the-hash-table)
      - [Adding taxonomic LCA information to the reference library](#adding-taxonomic-lca-information-to-the-reference-library)
    + [Taxonomic identification](#taxonomic-identification)
      - [Searching a query against a reference library](#searching-a-query-against-a-reference-library)
      - [Classification of reads](#classification-of-reads)
      - [Abundance profiling](#abundance-profiling)
      - [Contamination removal](#contamination-removal)
    + [Input and output formats](#input-and-output-formats)
      - [FASTA file for library construction](#fasta-file-for-library-construction)
      - [Query files for taxonomic identification](#query-files-for-taxonomic-identification)
      - [Taxonomy lookup and filename map](#taxonomy-lookup-and-filename-map)
      - [Report for number of matched *k*-mers and their Hamming distances](#report-for-number-of-matches-and-their-hamming-distances)
      - [Information for match distances and corresponding taxonomic LCAs](#information-for-match-distances-and-corresponding-taxonomic-lcas)
      - [Result files for classification and profiling](#result-files-for-classification-and-profiling)
    + [Description of CONSULT-II arguments and usage](#description-of-consult-ii-arguments-and-usage)
      - [`minimize`](#-minimize-)
      - [`consult_map`](#-consult-map-)
      - [`consult_search`](#-consult-search-)
      - [`consult_profile` & `consult_classify`](#-consult-profile-----consult-classify-)
    + [Public libraries](#public-libraries)
## Overview
CONSULT-II is a tool for taxonomic identification and successor of CONSULT.
CONSULT-II extends CONSULT by enabling taxonomic classification and abundance profiling, in addition to contamination removal.
The paper[^1] where we have described the design of the algorithm and software architecture is available [online](https://link.springer.com/chapter/10.1007/978-3-031-36911-7_13).

Relying on locality-sensitive hashing, CONSULT-II extracts *k*-mers from a query set and tests whether they fall within a user-specified hamming distance of *k*-mers in the reference dataset.
Using this invaluable information, it can derive taxonomic predictions for given query reads, or abundance estimation for samples.
Furthermore, for a given query read, CONSULT-II can report *k*-mer matches, Hamming distances to corresponding references, and *probabilistic* taxonomic least common ancestor (LCA) of matched reference *k*-mers.
It supports the inclusion of approximately billions *k*-mers in its reference library, accommodating datasets with tens of thousands of microbial species.
CONSULT-II is efficiently parallelized and can handle very large datasets.

Despite being memory-hungry, our careful benchmarking shows that CONSULT-II outperforms popular *k*-mer-based tools such as Kraken-2 and CLARK in accuracy.
We provide reference libraries, so you can download them and jump into taxonomic identification of your queries.

## Getting started
### System requirements
#### Memory and disk space
The exact memory footprint depends on the number of *k*-mers in a reference set, and the parameters used.
CONSULT-II requires enough free memory to hold the entire reference library in RAM, and this is needed during library construction and searching query *k*mers.
During, classification and abundance profiling, CONSULT-II reads matching information from the disk.
Using the default values (or heuristic), some examples of approximately required memory (conservative upper bounds) with respect to the number of *k*-mers in a reference set can be listed below;

* $2^{28}$ *k*-mers $\rightarrow$ $<5$ GB,

* $2^{30}$ *k*-mers $\rightarrow$ $<20$ GB,

* $2^{32}$ *k*-mers $\rightarrow$ $<80$ GB.

We note that during library construction the user will need slightly more RAM than the given values to accommodate intermediary processes (about an additional 10\%).
Once the reference library is built, these values should be sufficient.
For instance, the main reference library, with more than 10k species and about 8 billion *k*-mers leads to memory usage varying between 140GB and 150GB.
In addition to the reference library and metadata, CONSULT-II also needs to store *k*-mer match information per read to disk, which then be read during classification and profiling.
This shouldn't exceed input data file sizes (usually less than 10\% of input size).

#### Dependencies
CONSULT-II is a command-line tool implemented in C++11 with some x86 assembly code.
Many steps of CONSULT-II are embarrassingly parallelized, such as reference library reading, query search, classification, and profiling, using [OpenMP](https://www.openmp.org).
Compilation requires g++ that supports C++11.
For our tests, we have compiled versions 4.8.5 and 7.2.0, both of which work.
External tools such as [Jellyfish](http://www.genome.umd.edu/jellyfish.html) might be useful for reference library construction.

### Installation
1. Download using one of two approaches:
    - You can obtain the [zip file](https://github.com/noraracht/CONSULT/archive/main.zip) and extract the contents to a folder of your choice.
    Then, proceed to compilation.
    - Alternatively, you can clone the GitHub repository and continue with compilation.

2. To compile, go to the directory where core programs (*consult_\*.cpp*) are located and run the below commands.
    * You can use `make` to compile CONSULT.
    ```bash
    make all # for all components of CONSULT
    # OR
    make minimize # for the minimization script
    make map # for consult_map to construct a library
    make search # for consult_search to make queries
    make classify # for consult_classify to perform read-level classification from match info
    make profile # for consult_profile to perform abundance profiling from match info
    ```
    * Alternatively, you can run `g++` directly.
    ```bash
    g++ minimize.cpp -std=c++11 -o minimize # for the minimization script
    g++ consult_map.cpp -std=c++11 -O3 -o consult_map  # for consult_map to construct a library
    g++ consult_search.cpp -std=c++11 -fopenmp -O3 -o consult_search # for consult_search to make queries
    g++ consult_classify.cpp -std=c++11 -fopenmp -O3 -o consult_classify # for consult_classify to perform read-level classification from match info
    g++ consult_profile.cpp -std=c++11 -fopenmp -O3 -o consult_profile # for consult_profile to perform abundance profiling from match info
    ```

### Testing
To test your installation, and gain a better insight into input/output files and CONSULT-II's workflow, see the directory `example`.
To test all functionality, run `make all` in `example`.
As a result, the following files should be generated:

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
We made some libraries available online, and many more will be in the future. See [Public libraries](#public-libraries) for a list, just click on them to download. You may want to use `wget` or a similar tool. If you would like to use a pre-built library, skip this subsection and jump to [Taxonomic identification](#taxonomic-identification) subsection.
#### Preprocessing
We suggest the following workflow to obtain the *k*-mer lists from FASTA files.
To construct the CONSULT-II library from multiple assembly references, follow the steps below.

1. **To combine assembly references** into single file use ``cat`` as follows.
```bash
cat /path/to/folder/*.fna > combined.fna
```

2. **Extraction of canonical 35 bp *k*-mers** from FASTA genomic reference can be done with [Jellyfish](http://www.genome.umd.edu/jellyfish.html).
 - To compute *k*-mer profile run ``jellyfish count``.
```bash
jellyfish count -m 35 -s 100M -t 24 -C combined.fna -o counts.jf
```
 - To output a list of *k*-mers associated with their counts use ``jellyfish dump``.
 ```bash
jellyfish dump counts.jf > 35bp_kmer_lst.fa
```

3. **Minimization could be performed using custom C++11 script we provide.**
The script accepts a FASTA file (output of [Jellyfish](http://www.genome.umd.edu/jellyfish.html) as an input containing 35 bp canonical *k*-mers extracted from reference and outputs their 32 bp minimizers in FASTA format.
Run ``minimize`` as below.
```bash
./minimize -i 35bp_kmer_lst.fa -o 32bp_minzer_lst.fa
```

4. **Extraction of unique 32 bp *k*-mers** from the minimizers could be useful to further reduce *k*-mer count and remove duplicate sequences.
 - To compute *k*-mer profile run ``jellyfish count`` again.
```bash
jellyfish count -m 32 -s 100M -t 24 -C 32bp_minzer_lst.fa -o counts.jf
```
 - To export the list of *k*-mers and their counts use ``jellyfish dump``.
```bash
jellyfish dump counts.jf > 32bp_kmer_lst.fa
```

We note that, in our experiments, we used the minimization technique to decrease the total *k*-mer count in the original datasets and to reduce memory requirements.
Alternatively, if a dataset is small enough, and minimization is not needed, the user can use the last command directly to obtain a list of all 32 bp *k*-mers and utilize it as input into CONSULT-II software.

#### Construction of the hash table
To construct a standard reference library, use the following command to run `consult_map` as follows:
```bash
./consult_map -i /path/to/fasta_file -o $LIBRARY_NAME
# for example:
./consult_map -p 3 -t 2 -i k32C_af_mininimization.fa -o G000307305_nbr_map
```
where `-t` is tag size in bits, and determines the number of partitions ( $2^t$ ), `-p` is the Hamming distance threshold for a match, and `/path/to/fasta_file` is supposed to be the path of the input FASTA file, formatted as shown in the corresponding [section](#fasta-file-for-library-construction).
Replace `$LIBRARY_NAME` above with your preferred library name, which also will be the path of the reference library directory relative to the current working directory.
You can also give a path to another valid directory other than the current working directory: `$/path/to/$LIBRARY_NAME`.
Note that, CONSULT-II will not create parent directories of the given library path for you.
If the given directory path already contains a reference library with the same name, the software will throw an exception.
This feature is to prevent existing reference libraries from being overwritten.
In the reference library directory, several non-human-readable binary files will be stored, such as metadata, encoding arrays, and tag array (in chunks).

#### Adding taxonomic LCA information to the reference library
After constructing the hash table, and populating it with reference *k*-mers, we need to add taxonomic LCA label to each *k*-mer for abundance profiling and taxonomic classification.
For contamination removal, this step is not needed since matching to *any* reference sequence is sufficient in that case.
Two consecutive commands must be run, and filename and taxonomy look-up files must be given as arguments:
```bash
./consult_search -q /path/to/reference_genomes -i /path/to/$LIBRARY_NAME -o . --taxonomy-lookup-path /path/to/taxonomy-lookup --filename-map-path /path/to/filename-map --init-ID
./consult_search -q /path/to/reference_genomes -i /path/to/$LIBRARY_NAME -o . --taxonomy-lookup-path /path/to/taxonomy-lookup --filename-map-path /path/to/filename-map --update-ID
```
As the name suggests, `$LIBRARY_NAME` is the name of the library (also the name of the directory), and `/path/to/reference_genomes` is the path to the directory in which all the reference genomes are stored in FASTQ format separately with appropriate file names.
Note that, filenames are going to be used to map each genome file to a taxon (preferablly species), and parents of each taxon have to be listed in the lookup table.
See the [taxonomy lookup and filename map section](#taxonomy-lookup-and-filename-map) for a detailed description and instructions for generating them.
The fist command (with the flag `--initID`) initializes each label associated with *k*-mers, and counts how many distinct genomes have each *k*-mer.
Then, the second command (with the flag `update-ID`) computes and stores the probabilistic soft LCA labels for each *k*-mer.
Also, make sure to set `--thread-count` to a value as high as possible because this step might be pretty slow but efficiently parallelized.

### Taxonomic identification
#### Searching a query against a reference library
To query a set of sequences against a reference library, go to the directory where binaries are and execute the following CONSULT-II command:
```bash
 ./consult_search -i /path/to/$LIBRARY_NAME -q /path/to/query -o /path/to/output_directory
# for example:
./consult_search -q query000.fq -i g000307305_nbr_mapping -o . -c 1
```
The files containing query sequences should be in `/path/to/query` and in FASTQ format (one uncompressed `.fq`/`.fastq` file per each sample).
See [section](#query-files-for-taxonomic-identification) for a detailed description of the query path.
This step is required for all taxonomic identification tasks: classification, profiling, and contamination removal.
But, as described below, different flags have to be used for different tasks, and the above `consult_search` command is not useful per se.

#### Classification of reads
Classification of reads consists of two stages: i) finding *k*-mer matches of a given query (with Hamming distances and soft LCA labels) and ii) predicting a taxon based on the matching information.
For the first part run `consult_search` with the `--save-matches` flag, which will result in outputting *k*-mer matches, Hamming distances, and soft LCA labels in a file formatted as described in [this section](#information-for-match-distances-and-corresponding-taxonomic-lcas).
See an example command with needed flags below.
```bash
 ./consult_search -i /path/to/$LIBRARY_NAME -q /path/to/query -o /path/to/output_directory taxonomy-lookup-path /path/to/taxonomy-lookup --filename-map-path /path/to/filename-map --save-matches
# for example:
 ./consult_search -q query000.fq -i g000307305_nbr_mapping --taxonomy-lookup-path taxonomy-lookup --filename-map-path filename-map --save-matches
```

Then, after obtaining a list of matches with reference *k*-mers for each query, run `consult_classify` to summarize matching information with a final classification:
```bash
 ./consult_classify -i /path/to/match-info -o /path/to/output_directory taxonomy-lookup-path /path/to/taxonomy-lookup
# for example:
 ./consult_classify -i match-info_query000 -o . --taxonomy-lookup-path taxonomy-lookup
```
If successful, this command will generate a file populated with taxon predictions with corresponding total votes and the classified form of the reads (reverse-complement or original).
See [the corresponding section](result-files-for-classification-and-profiling) for the output format.
Filenames are input query filenames prefixed with *"classification_"*, for example, *"classification_query000"*.

#### Abundance profiling
Similar to classification, after finding *k*-mer matches of queries using `consult_search` with the `--save-matches` flag, we need to run `consult_profile` by passing the output of `consult_search` (in this case matching information) as input.
```bash
 ./consult_profile -i /path/to/match-info -o /path/to/output_directory taxonomy-lookup-path /path/to/taxonomy-lookup
# for example:
 ./consult_profile -i match-info_query000 -o . --taxonomy-lookup-path taxonomy-lookup
```
This command will output profile reports to `/path/to/output_directory`.
There will be independent profile vectors for species, genus, family, class, order, phylum, and kingdom.
Each rank will have its separate profile vector.
Note that, profiles of different ranks do not necessarily agree.
See [the corresponding section](result-files-for-classification-and-profiling) for the output format.
Filenames are input query filenames prefixed with *"profile_"* and postfixed with *"-"* & rank, for example, *"profile_query000-genus"*.

#### Contamination removal
For contamination removal, there is no need to run an additional command: `consult_search` also is able to generate files that contain the **classified reads** and **unclassified read**
To make CONSULT-II behave this way, give the `--classified-out` and `--unclassified-out` flags, correspondingly. The output file name will be prefixed with *"classified-seq_"* and *"classified-seq_"*.
Files are stored in the directory given in `/path/to/output_directory`, and the current working directory.
Every sample retains its original file name prefixed with *"classified-seq_"* or *"unclassified-seq_"*.
```bash
 ./consult_search -i $LIBRARY_NAME -q /path/to/query -o /path/to/output_directory --classified-out --unclassified-out --save-distances
# for example:
./consult_search -q query000.fq -i G000307305_nbr_mapping -o . -c 1 --classified-out --unclassified-out --save-distances
```
Another output that CONSULT-II can output (with `--save-distances` flag) is a tab-separated file, in which, each row is a read and the column values are the total numbers of *k*-mers of that read which matches with some *k*-mer in the reference library with Hamming distance $d$.
Each column corresponds to a distance value $d$.
The maximum distance value included in this file is determined by the $p$ value (Hamming distance threshold for a match) as a default ( $\lceil 1.5p \rceil$ ), but can be set to some other value with argument `--maximum-distance`. Filenames of the output will be query names prefixed with *"kmer-distances_"*.

Considering the example, unclassified reads would be stored in `unclassified-seq_G000307305.fq`.
If the flag `--classified-out` is given, classified reads would be stored in `classified-seq_G000307305.fq`.
If the flag `--save-distances` is given, distance values would be stored in `kmer-distances_G000307305.fq`.

### Input and output formats
#### FASTA file for library construction
For library construction, CONSULT was designed with [Jellyfish](http://www.genome.umd.edu/jellyfish.html)'s output in mind, which is a FASTA file representing a list of *k*-mers associated with their counts.
We tested with [Jellyfish](http://www.genome.umd.edu/jellyfish.html) 2.3.0.
See the [preprocessing section](#preprocessing) for details on how to generate the input file.
Note that CONSULT-II does not use the count values and the only relevant information is the sequence itself.
Jellyfish output is pseudo-randomly ordered, and thus, further randomization is not needed.
Note that the sequences should not be repeated.

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
The path `/path/to/query` can be a directory containing `.fastq` files or a `.fastq` file.
If it is a directory, each query file in that directory will be queried against the library, and separate outputs will be generated.
FASTA format is not supported at the moment.
Note, if you need to query FASTA files you can convert `.fasta`/`.fa` to `.fastq`/`.fq` using [seqtk](https://github.com/lh3/seqtk) `seqtk seq -F CHAR` command which attaches fake quality scores to the sequences.
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
A filename map consists of two space-separated columns and a row per genome/file.
The first column is for filenames without their extension, and the second column is for their taxa IDs (species if known).
This file will be used to add taxonomic information (soft LCA labels) with `consult_search` using `--update-ID` and `--init-ID` flags.
Essentially, this is to map reference genomes in `/path/to/reference_genomes` to taxa, an example row would be as below.
```
G000307305-k32C_minimized 1214065
```

Taxonomy lookup is an auxiliary file and it is a lookup table for parents of each taxon.
It is essentially the entire taxonomy, but formatted pragmatically.
Each row is for a taxon and consists of taxon ID, white space, and a comma-separated list of parents (from bottom to up).
Two example rows would look as below.
```
79885 79885,1386,186817,1385,91061,1239,2
1386 1386,186817,1385,91061,1239,2
```

Fortunately, we provide a script to create this lookup table from a taxonomy (`.dmp` file).
Simply run `scripts/construct_taxonomy_lookup.py --help` to start using the script.

#### Report for number of matches and their Hamming distances
CONSULT-II can report the number of *k*-mer matches for given query reads in a tab-separated file.
There will be `--maximum-distance`+1 many columns in addition to read ID and read form (original or reverse complement) columns.
Under each distance column, the number of *k*-mers matched with some reference *k*-mers is given.
An example (with `--maximum-distance` 5) is given below.
```
READ_ID    SEQ_TYPE    0    1    2    3    4    5
@NZ_DF158882.1-8256718    --    0    0    0    0    0    0
@NZ_DF158882.1-8256718    rc    1    0    0    0    0    0
```

#### Information for match distances and corresponding taxonomic LCAs
For classification and profiling, CONSULT-II relies on `consult_search`'s output of detailed information of *k*-mer matches, consisting of Hamming distance and taxonomic LCA for each *k*-mer match of a given query.
In this output of `consult_search`, each query read has three rows: i) read ID, ii) *k*-mer matches of the original direction and form iii) *k*-mer matches of the reverse complement.
You can find an example of *"match-info_\*"*.
```
@NZ_DF158883.1-1904898
-- 1214065:0 1214065:0 1214065:0
rc
```
Here, the read "@NZ_DF158883.1-190489", has three *k*-mer matches in its original form, and it has no match in its reverse complement form.
All three reference *k*-mers have the LCA taxon ID "1214065" (*Enterobacteriaceae bacterium B14*, species level LCA in this case), and the Hamming distances between query *k*-mers and matched reference *k*-mers are 0, i.e., we have exact matches.
As it can be easily seen, the format is `LCA_TAXON_ID:HAMMING_DISTANCE` per *k*-mer match.
In this case, the total vote for "1214065" would be $3.0$.

#### Result files for classification and profiling
Classification reports consist of five tab-separated columns and a row per given query read.
The first row is for the read ID, and the second column is to state the form (reverse complement or original) in which the read is classified.
Next to them, each row reports the classified taxon ID, the total-vote value of the classified taxon, and the total vote value of the root.
An example line from a classification report would look as below.
```
@NZ_DF158884.1-433795    rc    1214065    1.000000    1.000000
```

CONSULT-II reports separate profiles for each rank, and each file contains rows for all the taxa estimated to be present in the sample.
These reports contain three tab-separated columns: the first column is for taxa IDs, the second column is for taxonomic ranks, and the final column is for abundance value.
The abundance value can be interpreted as the total number of reads associated with the corresponding taxon.
The sum of abundance values gives the number of reads with at least a single *k*-mer match.
You can simply normalize this column to have a vector whose values sum up to 1.
An example is given below.
```
TAXONOMY_ID    TAXONOMY_LEVEL    FRACTION_TOTAL
1224    phylum    28472.0000000000
```

### Description of CONSULT-II arguments and usage
#### `minimize`

- `-i` or `--input-fasta-file`: input `.fasta` file containing canonical *k*-mers, default length is 35.
- `-o` or `--output-fasta-file`: `.fasta` file to output minimizers of the given canonical *k*-mers, default length is 32.

#### `consult_map`
- `-i` or `--input-fasta-file`: input `.fasta` file to construct library.
- `-o` or `--output-library-dir`: output path to the directory that will constitute the CONSULT-II library.
- `-h` or `--number-of-positions`: number of randomly positioned bits to compute LSH.
- `-t` or `--tag-size`: number of bits to be used as tag.
- `-l` or `--number-of-tables`: number of tables, i.e., number of hash functions.
- `-b` or `--column-per-tag`: number of columns per each tag partition, i.e., number of *k*-mers each encoding can map to.

#### `consult_search`
- `-i` or `input-library-dir`: directory of the CONSULT-II library that will be used as the reference library.
- `-o` or `output-result-dir`: directory in which all output files (classified and unclassified reads, matching information) will be saved.
- `-q` or `--query-path`: the path to the query file, or to the directory containing query files (or reference genomes to add taxonomic information to input library).
- `-c` or `--number-of-matches`: the minimum number of matched *k*-mers that is required to call sequencing read classified.
For instance, if at least one *k*-mer match is enough to classify a read (default setting mentioned in a paper), `-c` should be set to 1 in the software.
If at least two *k*-mer matches are required to call the entire read a match, `-c` should be set to 2.
The default value is 1.
- `--thread-count`: number of threads to be used, default is 1.
- `--unclassified-out`: to output reads that are unclassified in a file with a name query file name prefixed with *"unclassified-seq_"*.
This is given by default.
- `--classified-out`: to output reads that are classified in a file with a name query file name prefixed with *"classified-seq_"*.
- `--save-distances`: to save the number of matched *k*-mers in a tab-separated file where columns are the distances of corresponding counts.
The file name is a query file name prefixed with *"kmer-distances_"*. See the corresponding [section](#report-for-number-of-matched-k-mers-and-their-hamming-distances).
- `--maximum-distance`: maximum distance to be included as a column in the file containing *k*-mer match counts with respect to varying Hamming distance values.
Note that, when the `--maximum-distance` value is too large compared to $p$, CONSULT-II does not necessarily aim to find such *k*-mers to compute column values.
This is because the library size and the corresponding theoretical guarantees only consider $p$ value.
So, if one would like to use a large `--maximum-distance` value, the $p$ value should be increased proportionally.
- `init-ID`: to initialize taxonomic LCA labels (set them to the root of the tree), and to count distinct genomes in which each *k*-mer appears.
- `update-ID`: to compute and store probabilistic LCA labels for each *k*-mer in the reference library.
It is required to run `consult_search` with this flag (and `--init-ID` before classification and profiling.
- `--taxonomy-lookup-path`: the path of the taxonomic LCA lookup table in a human-readable format.
See the corresponding [section](#taxonomy-lookup-and-filename-map) for details.
- `--filename-map-path`: the path of the filename/genome to species ID map in a human-readable format.
See the corresponding [section](#taxonomy-lookup-and-filename-map) for details.
- `--save-matches`: to output detailed *k*-mer match information of query sequences consisting of Hamming distances and taxonomic LCAs.
Output filename convention is query file name prefixed with *"kmer-distances_"*.
See the corresponding [section](#information-for-match-distances-and-corresponding-taxonomic-lcas) for details.

#### `consult_profile` & `consult_classify`
- `-i` or `--input-matches-path`: path to input match information file, or path to a directory for multiple files.
Note that if this is a directory, it should only contain *" match-info_"* files as it tries to iterate over all files in the directory.
- `-o` or `--output-predictions-dir`: output directory to write classification results for each input query match information. Query names will be prefixed with *"classification_"*.
- `--taxonomy-lookup-path`: the path of the taxonomic LCA lookup table in a human-readable format.
- `--thread-count`: number of threads to be used, default is 1.

### Public libraries
- [WoL: Reference Phylogeny for Microbes (bacteria and archaea) (140 Gb - large but performant with defaults)](https://ter-trees.ucsd.edu/data/consult/CONSULT-II/library-v020-WoL140G.tar.gz)
- [WoL: Reference Phylogeny for Microbes (bacteria and archaea) (32 Gb - lighter-weight but still highly accurate)](https://ter-trees.ucsd.edu/data/consult/CONSULT-II/library-v030-WoL32G.tar.gz)
- [WoL: Reference Phylogeny for Microbes (bacteria and archaea) (18 Gb - lightweight and robust)](https://ter-trees.ucsd.edu/data/consult/CONSULT-II/library-v030-WoL18G.tar.gz)
- [Taxonomy lookup tables for all WoL-v1 libraries](https://github.com/bo1929/shared.CONSULT-II/raw/master/misc/taxonomy_lookup-RefSeq2019.tar.gz)

[^1]: Şapcı, A.O.B., Rachtman, E., Mirarab, S. (2023). CONSULT-II: Taxonomic Identification Using Locality Sensitive Hashing. In: Jahn, K., Vinař, T. (eds) Comparative Genomics. RECOMB-CG 2023. Lecture Notes in Computer Science, vol 13883. Springer, Cham. https://doi.org/10.1007/978-3-031-36911-7_13
