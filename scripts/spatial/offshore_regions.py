import geopandas as gpd
import pandas as pd
from shapely.geometry import Point

offshore_df = pd.read_csv(snakemake.input[0])
offshore_gdf = gpd.GeoDataFrame(
    offshore_df,
    geometry=[Point(xy) for xy in zip(offshore_df.Longitude, offshore_df.Latitude)],
    crs="EPSG:4326",
)

utm_crs = "EPSG:3035"  # projected CRS for buffering
offshore_gdf_proj = offshore_gdf.to_crs(utm_crs)
offshore_gdf_proj["geometry"] = offshore_gdf_proj.geometry.buffer(1000)  # 5km buffer
offshore_gdf_buffered = offshore_gdf_proj.to_crs("EPSG:4326")

offshore_gdf_buffered.to_file(snakemake.output[0], driver="GPKG")
