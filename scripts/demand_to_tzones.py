import pandas as pd
import geopandas as gpd
from shapely.geometry import Point


demand_df = pd.read_csv(snakemake.input[0])
gdf =gpd.read_file(snakemake.input[1])

demand_df['geometry'] = demand_df.apply(lambda row: Point(row['Longitude'], row['Latitude']), axis=1)
FES_gdf = gpd.GeoDataFrame(demand_df, geometry='geometry', crs="EPSG:4326")

print(gdf.columns)
print(gdf.head())

print(gdf.crs)           # Coordinate Reference System
print(gdf.geometry.name)


demand_df.to_csv(snakemake.output[0], index=False)