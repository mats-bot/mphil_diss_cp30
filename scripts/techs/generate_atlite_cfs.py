
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
      "This process will take a while due to the large amount of data needed. \n")


zones = gpd.read_file(snakemake.input[0])
year = int(snakemake.wildcards.year)

cutout_dir = snakemake.output[0]
os.makedirs(cutout_dir, exist_ok=True)

time_ranges = []
for month in range(1, 13):
    last_day = calendar.monthrange(year, month)[1]  # returns correct last day of month
    time_ranges.append(slice(f"{year}-{month:02d}-01", f"{year}-{month:02d}-{last_day:02d}"))
months = [f"{m:02d}" for m in range(1, 13)]


solar_dfs = []
wind_dfs = []

for i, time_range in enumerate(time_ranges):
    cutout_path = os.path.join(cutout_dir, f"cutout_{months[i]}.nc")

    if os.path.exists(cutout_path):
        print(f"Using existing cutout: {cutout_path}")
        cutout = atlite.Cutout(path=cutout_path)
    else:
        print(f"Creating cutout: {cutout_path}")
        cutout = atlite.Cutout(
            path=cutout_path,
            module="era5",
            bounds=zones.total_bounds.tolist(),
            time=time_range
        )
        cutout.prepare()

    # --- Solar ---
    solar_cf = cutout.pv(
        panel="CSi",
        orientation="latitude_optimal",
        shapes=zones
    ).to_pandas()
    solar_cf.columns = zones['z1'].values
    solar_cf = solar_cf[sorted(solar_cf.columns, key=lambda x: int(x[1:]))]
    solar_dfs.append(solar_cf)

    # --- Wind ---
    wind_cf = cutout.wind(
        turbine="Vestas_V112_3MW",
        shapes=zones,
        add_cutout_windspeed=True
    ).to_pandas()
    wind_cf.columns = zones['z1'].values
    wind_cf = wind_cf[sorted(wind_cf.columns, key=lambda x: int(x[1:]))]
    wind_dfs.append(wind_cf)

# --- Write final merged files ---
pd.concat(solar_dfs).to_csv(snakemake.output[1])  # solar CSV
pd.concat(wind_dfs).to_csv(snakemake.output[2])   # wind CSV




## Offshore (iterate over projects instead of zones)

offshore_df = pd.read_csv(snakemake.input[1])
offshore_gdf = gpd.GeoDataFrame(
    offshore_df,
    geometry=[Point(xy) for xy in zip(offshore_df.Longitude, offshore_df.Latitude)],
    crs="EPSG:4326"
)

utm_crs = "EPSG:3035"  # projected CRS for buffering
offshore_gdf_proj = offshore_gdf.to_crs(utm_crs)
offshore_gdf_proj['geometry'] = offshore_gdf_proj.geometry.buffer(5000)  # 5km buffer
offshore_gdf_buffered = offshore_gdf_proj.to_crs("EPSG:4326")

bounds2 = offshore_gdf_buffered.total_bounds.tolist()

offshore_cutout_dir = snakemake.output[3]
os.makedirs(offshore_cutout_dir, exist_ok=True)

months = [f"{m:02d}" for m in range(1, 13)]

offshore_dfs = []

for month in range(1, 13):
    last_day = calendar.monthrange(year, month)[1]
    time_range_c = slice(f"{year}-{month:02d}-01", f"{year}-{month:02d}-{last_day:02d}")

    cutout_path = os.path.join(offshore_cutout_dir, f"offshore_cutout_{month:02d}.nc")

    if os.path.exists(cutout_path):
        print(f"Using existing offshore cutout: {cutout_path}")
        cutout2 = atlite.Cutout(path=cutout_path)
    else:
        print(f"Creating offshore cutout: {cutout_path}")
        cutout2 = atlite.Cutout(
            path=cutout_path,
            module="era5",
            bounds=bounds2,
            time=time_range_c
        )
        cutout2.prepare()

    cf_offshore = cutout2.wind(
        turbine="NREL_ReferenceTurbine_2020ATB_15MW_offshore",
        shapes=offshore_gdf_buffered,
        add_cutout_windspeed=True
    )

    df_month = cf_offshore.to_pandas()
    df_month.columns = offshore_gdf_buffered['Site Name'].values
    offshore_dfs.append(df_month)

# Concatenate all months and save final offshore CF CSV
pd.concat(offshore_dfs).to_csv(snakemake.output[4])



# time_range_c = slice(f"{year}-01-01", f"{year}-12-31")

# #  offshore wind farm to enable calliope to decide which to build
# offshore_df = pd.read_csv(snakemake.input[1])
# offshore_gdf = gpd.GeoDataFrame(
#     offshore_df,
#     geometry=[Point(xy) for xy in zip(offshore_df.Longitude, offshore_df.Latitude)],
#     crs="EPSG:4326"
#     )

# utm_crs = "EPSG:3035" # expanding small radius (5km) around points
# offshore_gdf_proj = offshore_gdf.to_crs(utm_crs)#
# offshore_gdf_proj['geometry'] = offshore_gdf_proj.geometry.buffer(5000)
# offshore_gdf_buffered = offshore_gdf_proj.to_crs("EPSG:4326")

# bounds2 = offshore_gdf_buffered.total_bounds.tolist()

# cutout2 = atlite.Cutout(
#     path=snakemake.output[6],
#     module="era5",
#     bounds=bounds2,
#     time=time_range_c
# )
# cutout2.prepare()


# cf_offshore = cutout2.wind(
# turbine="NREL_ReferenceTurbine_2020ATB_15MW_offshore",
# shapes=offshore_gdf_buffered,
# add_cutout_windspeed=True
# )

# df = cf_offshore.to_pandas()
# df.columns = offshore_gdf_buffered['Site Name'].values  
# df.to_csv(snakemake.output[4])

