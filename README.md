# mphil_diss_cp30

## Install

```sh
conda env create -f environment.yaml
conda activate calliope-CP30-rep
```

## Dry run (check what will run)

```sh
snakemake --dry-run data/intermediates/GSP_timeseries.csv
```

## Run

This will run the `default_target` in the Snakefile.

```sh
snakemake
```

## Create a rulegraph to show the order of events

```sh
snakemake rulegraph
```
