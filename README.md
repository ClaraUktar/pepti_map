# pepti_map

`pepti_map` is a tool for mapping peptide sequences to their possible genomic loci. It does so based only on the sequence information, so that no further peptide information besides the amino acid sequence and no annotation of the genome is required, enabling a mapping to personal genomes. `pepti_map` utilizes RNA-seq reads of the same sample to facilitate the mapping: Each peptide is first matched to RNA-seq reads, which are then assembled into longer contigs and aligned onto the genome.

## Setup

To setup `pepti_map`, first install all dependencies listed in the `environment_<os>.yml` (`environment_macos.yml` or `environment_linux.yml`). We recommend using [Mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html#mamba-install). For example:

```
$ mamba env create -f environment_linux.yml
```

`pepti_map` relies on the following tools to be installed:
- [`Trinity`](https://github.com/trinityrnaseq/trinityrnaseq/wiki)
- [`GMAP`](http://research-pub.gene.com/gmap/)
- [`PoGo`](https://www.sanger.ac.uk/tool/pogo/)

If you install the dependencies via the given `environment_<os>.yml`, both `Trinity` and `GMAP` should already be installed. To install `PoGo`, download the latest release from its [GitHub page](https://github.com/cschlaffner/PoGo/releases).

You will then need to add the path to the directory in which your `PoGo` installation is located in an `.env` file via the `POGO_PATH` environment variable, e.g.:

```
POGO_PATH=/Users/me/Tools/PoGo_v1.2.3/Linux
```

This is the only environment variable that needs to be set for `pepti_map` to work. You can, however, set additional environment variables. Below you will find a table listing all environment variables.

| Environment Variable Name | Usage |
| ------------------------- | ----- |
| `POGO_PATH`               | The path to the directory in which the `PoGo` installation is located. |
| `IO_N_PROCESSES`          | The number of processes to use when generating the input files for `PoGo`. If not set, defaults to `multiprocessing.cpu_count()`. |
| `TRINITY_USE_DOCKER`      | Whether to run a dockerized version of Trinity. Value must be `True` or `False`. If not set, defaults to `False`. If set to `True`, a dockerized version of Trinity must be installed on the system. |
| `TRINITY_PATH`            | The path to the `Trinity` installation. If not given, expects `Trinity` to be executable from the working directory (e.g. by using an installation via a `Mamba` environment). |
| `TRINITY_N_PROCESSES`     | The number of processes with which to run `Trinity` in parallel. If not set, defaults to `multiprocessing.cpu_count()`. |
| `GMAP_N_THREADS`          | The number of threads with which to run `GMAP` during the alignment. Corresponds to the `-t` option of `gmap`. If not set, defaults to `multiprocessing.cpu_count()`. |
| `GMAP_BATCH_MODE`         | The batch mode in which to run `GMAP` during the alignment. Corresponds to the `-B` option of `gmap`. If not set, defaults to 2.      |


## Usage

To run `pepti_map`, first activate your `Mamba` environment.

The general usage of `pepti_map` is as follows:
```
$ (peptimap-env) python -m pepti_map.main [OPTIONS]
```
A minimal command would look like this:
```
$ (peptimap-env) python -m pepti_map.main -p path/to/peptide/file -r path/to/rna/reads/file -g path/to/genome/fasta
```
In order for `pepti_map` to be able to run, at least the `-p` and `-r` options need to be set, as well as one of the `-g` or `-x` option.

For further specification of these and further options, see the table below or use
```
$ (peptimap-env) python -m pepti_map.main --help
```

Overview of the options for running `pepti_map`:

| Option | Usage |
| ------ | ----- |
| `-p` / `--peptide-file` | The path to the peptide file (for format, see below). |
| `-r` / `--rna-file` | The path to the RNA-seq file. In case of paired-end sequencing, this file is expected to be in forward orientation. |
| `-pa` / `--paired-end-file` | The path to the second RNA-seq file in case of paired-end sequencing. This file is expected to be in reverse orientation. If none is given, the RNA-seq file given with the `-r` option is assumed to result from single-end sequencing. |
| `-c` / `--cutoff` | The position of the last base in the reads after which a cutoff should be performed (starting at 1). The cutoff is applied to all reads. If the value is equal to or smaller than 0, no cutoff is performed. (Default: -1) |
| `-k` / `--kmer-length` | The k-mer size used during the mapping of peptides to RNA-seq reads. As the RNA-seq reads are 3-frame translated for the mapping, the k-mer size refers to amino acids. (Default: 7) |
| `-o` / `--output-dir` | The path to the output directory for all generated files. (Default: `./`)|
| `-pi` / `--precompute-intersections` | If used, the intersection sizes for the Jaccard Index calculation are precomputed during the matching phase. |
| `-j` / `--jaccard-index-threshold` | Sets of matched RNA-seq reads per peptide will only be merged together if their Jaccard Index has a value above the given threshold. (Default: 0.5) |
| `-m` / `--merging-method` | Which merging method to use for sets of matched RNA-seq reads. Must be one of `agglomerative-clustering`, `full-matrix`. (Default: `full-matrix`) |
| `-cl` / `--min-contig-length` | Sets the `--min_contig_length` option for Trinity during assembly. A value below 100 is not possible. (Default: 100) |
| `-g` / `--genome` | The path to the genome file(s) to align to. In case of multiple files, the paths must be separated by comma. |
| `-x` / `--gmap-index` | The path to an existing GMAP index that should be used instead of building a new one. If this option is set, the `-g` / `--genome` option is ignored.|

### Input file formats

#### Peptide file ( `-p` / `--peptide-file`):
This file should contain a list of peptides in form of amino acid sequences, with one peptide sequence per line. Optionally, the file may contain protein group information per peptide, with this information being on the same line as the peptide, separated by tab. A header is not needed. One line may thus look as follows:
```
<peptide sequence>
```
or
```
<peptide sequence>  <protein group information>
```

If the protein group information is given, peptides with the same protein group will be grouped together, with matches to the RNA-seq reads being allocated per group. If not given, each peptide is treated as a separate group.

#### RNA-seq file(s) (`-r` / `--rna-file` and `-pa` / `--paired-end-file`):
These files should contain the RNA-seq reads of the same sample/subject as the peptide data in FASTQ or gzipped FASTQ format.

#### Genome file(s) (`-g` / `--genome`):
These files should contain the genomic sequences to align to in FASTA format.

### `pepti_map` output
`pepti_map` will output a GTF (`pepti_map_output.gtf`) and a BED (`pepti_map_output.bed`) file containing the mappings of the peptides to genomic loci. Peptides that were not matched will not appear in these files. For further information on these formats, please see the `PoGo` [documentation](https://github.com/cschlaffner/PoGo), as these outputs correspond directly to the respective `PoGo` output formats.

Additionally, `pepti_map` outputs a quantification file (`peptide_read_quant.tsv`), containing the group id assigned to each peptide, the number of reads that were matched to the peptide, and the number of reads that were matched to the group the peptide belongs to. This means, for all peptides belonging to the same group, their individual read counts add up exactly to the read count for the group. In case a peptide is too short (smaller than the k-mer length), it will be assigned a group id of -1 and excluded from further calculations. The format per line is as follows, with the individual values being tab-separated:
```
<peptide sequence>  <numeric group id>  <number of matched reads for peptide>   <number of matched reads per group>
```
