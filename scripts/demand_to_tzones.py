import pandas as pd
import geopandas as gpd
from shapely.geometry import Point



demand_df = pd.read_csv(snakemake.input[0])
zones_gdf =gpd.read_file(snakemake.input[1])

demand_gdf = gpd.GeoDataFrame(
    demand_df,
    geometry=gpd.points_from_xy(demand_df['Longitude'], demand_df['Latitude']),
    crs="EPSG:4326"
)

# Ensure same reference system
zones_gdf = zones_gdf.to_crs("EPSG:4326")

# Assign zone to each data point
joined_gdf = gpd.sjoin(
    demand_gdf, 
    zones_gdf[['z1', 'geometry']],  
    how='left',  
    predicate='within'  
)

# Handle points outside of the zones 
unassigned_mask = joined_gdf['z1'].isna()
unassigned_gdf = demand_gdf[unassigned_mask].copy()

if not unassigned_gdf.empty:  
    # Iterates over unassigned points to find closest zone   
    for idx, point in unassigned_gdf.geometry.items():
        nearest_zone = None
        min_distance = float('inf')
        

        for zone_idx, zone in zones_gdf.iterrows():
            distance = point.distance(zone.geometry)
            if distance < min_distance:
                min_distance = distance
                nearest_zone = zone['z1']
        joined_gdf.loc[idx, 'z1'] = nearest_zone



demand_df['zone'] = joined_gdf['z1']

demand_df.to_csv(snakemake.output[0], index=False)