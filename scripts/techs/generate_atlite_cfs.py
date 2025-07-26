
import geopandas as gpd
import pandas as pd
from shapely.geometry import Point
from pathlib import Path 
import sys

cdsapirc = Path.home() / ".cdsapirc"
if not cdsapirc.exists():
    sys.exit(
        "ERROR: CDS API key not found.\n"
        "This model relies on CDS data for which an API key is necessary.\n"
        "You can find instructions on getting one at:\n"
        "https://cds.climate.copernicus.eu/how-to-api"
    )


import atlite
print(f"Using atlite version: {atlite.__version__}")



zones = gpd.read_file(snakemake.input[0])
year = int(snakemake.wildcards.year)  # weather year
cutout1_path, cutout2_path, solar_out, onshore_out, offshore_out = snakemake.output


time_range = slice(f"{year}-01-01", f"{year}-01-31")

cutout = atlite.Cutout(
                path=snakemake.output[0],
                module="era5",
                bounds=zones.total_bounds.tolist(),
                time=time_range)
cutout.prepare()

solar_cf = cutout.pv(
    panel='CSi',  # defaults
    orientation='latitude_optimal',
    shapes=zones
)
onshore_wind_cf = cutout.wind(
    turbine='Vestas_V112_3MW',
    shapes=zones,
    add_cutout_windspeed=True
)

solar_cf.to_series().to_csv(snakemake.output[2], index=None)
onshore_wind_cf.to_series().to_csv(snakemake.output[3], index=None)


#  offshore wind farm to enable calliope to decide which to build
offshore_df = pd.read_csv(snakemake.input[1])
offshore_gdf = gpd.GeoDataFrame(
    offshore_df,
    geometry=[Point(xy) for xy in zip(offshore_df.Longitude, offshore_df.Latitude)],
    crs="EPSG:4326"
    )
bounds2 = offshore_gdf.total_bounds.tolist()

utm_crs = "EPSG:3035" # expanding small radius around points
offshore_gdf_proj = offshore_gdf.to_crs(utm_crs)#
offshore_gdf_proj['geometry'] = offshore_gdf_proj.geometry.buffer(5000)
offshore_gdf_buffered = offshore_gdf_proj.to_crs("EPSG:4326")




cutout2 = atlite.Cutout(
    path=snakemake.output[1],
    module="era5",
    bounds=bounds2,
    time=time_range
)
cutout2.prepare()


cf_offshore = cutout2.wind(
turbine="NREL_ReferenceTurbine_2020ATB_15MW_offshore",
shapes=offshore_gdf_buffered,
add_cutout_windspeed=True
)

df = cf_offshore.to_pandas()
df.columns = offshore_gdf_buffered['Site Name'].values  
df.to_csv(snakemake.output[4])

