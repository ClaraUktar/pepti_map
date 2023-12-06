# pepti_map

TODO: Short explanation of what it does

## Setup

To setup `pepti_map`, first install all dependencies listed in the `environment.yml`. We recommend using [Mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html#mamba-install).

```
$ mamba env create -f environment.yml
```

`pepti_map` relies on the following tools to be installed:
- [`Trinity`](https://github.com/trinityrnaseq/trinityrnaseq/wiki)
- [`GMAP`](http://research-pub.gene.com/gmap/)
- [`PoGo`](https://www.sanger.ac.uk/tool/pogo/)

If you install the dependencies via the given `environment.yml`, both `Trinity` and `GMAP` should already be installed. To install `PoGo`, download the latest release from its [`GitHub` page](https://github.com/cschlaffner/PoGo/releases).

You will then need to add the path to the directory in which your `PoGo` installation is located in an `.env` file via the `POGO_PATH` environment variable, e.g.:

```
POGO_PATH=/Users/me/Tools/PoGo_v1.2.3/Linux
```

This is the only environment variable that needs to be set for pepti_map to work. You can, however, set additional environment variables. Below you will find a table listing all environment variables.

| Environment Variable Name | Usage |
| ------------------------- | ----- |
| `POGO_PATH`               | The path to the directory in which the `PoGo` installation is located. |
| `IO_N_PROCESSES`          | The number of processes to use when generating the input files for `PoGo`. If not set, defaults to `multiprocessing.cpu_count()`. |
| `TRINITY_USE_DOCKER`      | Whether to run a dockerized version of Trinity. Value must be `True` or `False`. If not set, defaults to `False`. If set to `True`, a dockerized version of Trinity must be installed on the system. |
| `TRINITY_PATH`            | The path to the `Trinity` installation. If not given, expects `Trinity` to be executable from the working directory (e.g. by using an installation via a `Mamba` environment). |
| `TRINITY_N_PROCESSES`     | The number of processes with which to run `Trinity` in parallel. If not set, defaults to `multiprocessing.cpu_count()`. |
| `GMAP_N_THREADS`          | The number of threads with which to run `GMAP` during the alignment. Corresponds to the `-t` option of `gmap`. If not set, defaults to `multiprocessing.cpu_count()`. |
| `GMAP_BATCH_MODE`         | The batch mode in which to run `GMAP` during the alignment. Corresponds to the `-B` option of `gmap`. If not set, defaults to `2`.      |


## Usage

TODO
