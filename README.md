# mphil_diss_cp30


## Requirements

```sh
This model uses weather data from the Climate Data Store, which requires a personal API key.
Instructions for setup can be found here:
https://cds.climate.copernicus.eu/how-to-api

Additionally, the terms of use for each dataset must be accepted for the key to be valid.
This model uses era5 datasets, for which terms of use can be found and accepted here:
https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels?tab=overview


```

## Install

```sh
conda env create -f environment.yaml
conda activate calliope-CP30-rep
```

## Dry run (check what will run)

```sh
snakemake --dry-run --cores 1 results/model_results_B1.nc results/model_results_S1_NC.nc
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
