
import geopandas as gpd
import pandas as pd
from shapely.geometry import Point
from pathlib import Path 
import sys
import calendar
import os 

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
print(f'\n\n If the weather year has not been run before or the output files are missing, \n',
      "this process will take a while due to the large amount of data needed. \n")


zones = gpd.read_file(snakemake.input[0])
year = int(snakemake.wildcards.weather_year)

cutout_paths = snakemake.output["onshore_cutouts"]
cutout_dir = os.path.dirname(cutout_paths[0])
os.makedirs(cutout_dir, exist_ok=True)

time_ranges = []
for month in range(1, 13):
    last_day = calendar.monthrange(year, month)[1]  # returns correct last day of month
    time_ranges.append(slice(f"{year}-{month:02d}-01", f"{year}-{month:02d}-{last_day:02d}"))
months = [f"{m:02d}" for m in range(1, 13)]

# Checking if cutouts exist so they dont have to be recomputed
existing_cutouts = [os.path.exists(p) for p in cutout_paths]


solar_dfs = []
wind_dfs = []

for i, time_range in enumerate(time_ranges):
    cutout_path = cutout_paths[i]

    if existing_cutouts[i]:
        print(f"Using existing cutout: {cutout_path}")
        cutout = atlite.Cutout(path=cutout_path)
    else:
        print(f"Creating cutout: {cutout_path}")
        cutout = atlite.Cutout(
            path=cutout_path,
            module="era5",
            bounds=zones.total_bounds.tolist(),
            time=time_range,
            features=["pv", "onshore_wind"]
        )
        cutout.prepare()


    solar_cf = cutout.pv(
        panel="CSi",
        orientation="latitude_optimal",
        shapes=zones,
    ).to_pandas()
    solar_cf.columns = zones['z1'].values
    solar_cf = solar_cf[sorted(solar_cf.columns, key=lambda x: int(x[1:]))]
    solar_dfs.append(solar_cf)


    wind_cf = cutout.wind(
        turbine="NREL_ReferenceTurbine_2020ATB_7MW",
        shapes=zones,
        add_cutout_windspeed=True
    ).to_pandas()
    wind_cf.columns = zones['z1'].values
    wind_cf = wind_cf[sorted(wind_cf.columns, key=lambda x: int(x[1:]))]
    wind_dfs.append(wind_cf)

pd.concat(solar_dfs).div(100).to_csv(snakemake.output["solar_cf"])
pd.concat(wind_dfs).div(100).to_csv(snakemake.output["onshore_cf"]) 




## Offshore (iterate over projects instead of zones)

offshore_df = pd.read_csv(snakemake.input[1])
offshore_gdf = gpd.GeoDataFrame(
    offshore_df,
    geometry=[Point(xy) for xy in zip(offshore_df.Longitude, offshore_df.Latitude)],
    crs="EPSG:4326"
)

utm_crs = "EPSG:3035"  # projected CRS for buffering
offshore_gdf_proj = offshore_gdf.to_crs(utm_crs)
offshore_gdf_proj['geometry'] = offshore_gdf_proj.geometry.buffer(1000)  # 5km buffer
offshore_gdf_buffered = offshore_gdf_proj.to_crs("EPSG:4326")

bounds2 = offshore_gdf_buffered.total_bounds.tolist()

offshore_cutout_paths = snakemake.output["offshore_cutouts"]
offshore_cutout_dir = os.path.dirname(offshore_cutout_paths[0])
os.makedirs(offshore_cutout_dir, exist_ok=True)

months = [f"{m:02d}" for m in range(1, 13)]

existing_offshore_cutouts = [os.path.exists(p) for p in offshore_cutout_paths]


offshore_dfs = []

for i, month in enumerate(range(1, 13)):
    last_day = calendar.monthrange(year, month)[1]
    time_range_c = slice(f"{year}-{month:02d}-01", f"{year}-{month:02d}-{last_day:02d}")

    cutout_path = offshore_cutout_paths[i]

    if existing_offshore_cutouts[i]:
        print(f"Using existing offshore cutout: {cutout_path}")
        cutout2 = atlite.Cutout(path=cutout_path)
    else:
        print(f"Creating offshore cutout: {cutout_path}")
        cutout2 = atlite.Cutout(
            path=cutout_path,
            module="era5",
            bounds=bounds2,
            time=time_range_c,
            features=["offshore_wind"]
        )
        cutout2.prepare()

    cf_offshore = cutout2.wind(
        turbine="NREL_ReferenceTurbine_2020ATB_15MW_offshore",
        shapes=offshore_gdf_buffered,
        add_cutout_windspeed=True,
        capacity_factor=True
    )

    df_month = cf_offshore.to_pandas()
    df_month.columns = offshore_gdf_buffered['Site Name'].values
    offshore_dfs.append(df_month)


pd.concat(offshore_dfs).to_csv(snakemake.output["offshore_cf"])

