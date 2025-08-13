import sys
from pathlib import Path

import atlite
import geopandas as gpd

cdsapirc = Path.home() / ".cdsapirc"
if not cdsapirc.exists():
    sys.exit(
        "ERROR: CDS API key not found.\n"
        "This model relies on CDS data for which an API key is necessary.\n"
        "You can find instructions on getting one at:\n"
        "https://cds.climate.copernicus.eu/how-to-api"
    )


print(f"Using atlite version: {atlite.__version__}")
print(
    "\n\n If the weather year has not been run before or the output files are missing, \n",
    "this process will take a while due to the large amount of data needed. \n",
)


zones = gpd.read_file(snakemake.input[0])

cutout = atlite.Cutout(
    path=snakemake.output[0],
    module="era5",
    bounds=zones.total_bounds.tolist(),
    time=f"{snakemake.wildcards.weather_year}-{snakemake.wildcards.month}",
    features=snakemake.params.features,
)
cutout.prepare()
