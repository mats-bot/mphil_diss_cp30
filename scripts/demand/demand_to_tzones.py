import pandas as pd
import geopandas as gpd
from shapely.geometry import Point

def assign_zones(demand_csv, zones_file, output_csv):
    """
    Assigns each demand point to CP30 transmission zones.
    
    Parameters:
    - demand_csv: CSV with 'Longitude' and 'Latitude' columns.
    - zones_file: Path to a GeoPackage or shapefile with zone geometries and a 'z1' column.
    - output_csv: Path to write the resulting CSV with an additional 'zone' column.
    """
    demand_df = pd.read_csv(demand_csv)
    zones_gdf = gpd.read_file(zones_file)
    
    demand_gdf = gpd.GeoDataFrame(
        demand_df,
        geometry=gpd.points_from_xy(demand_df['Longitude'], demand_df['Latitude']),
        crs="EPSG:4326"
    )
    
    # Ensure same reference system
    zones_gdf = zones_gdf.to_crs("EPSG:4326")
    
    # Spatial join
    joined_gdf = gpd.sjoin(
        demand_gdf, 
        zones_gdf[['z1', 'geometry']],  
        how='left',  
        predicate='within'  
    )
    
    # Handle points outside of zones
    unassigned_mask = joined_gdf['z1'].isna()
    unassigned_gdf = demand_gdf[unassigned_mask].copy()
    
    if not unassigned_gdf.empty:  
        for idx, point in unassigned_gdf.geometry.items():
            nearest_zone = None
            min_distance = float('inf')
            
            for _, zone in zones_gdf.iterrows():
                distance = point.distance(zone.geometry)
                if distance < min_distance:
                    min_distance = distance
                    nearest_zone = zone['z1']
            
            joined_gdf.loc[idx, 'z1'] = nearest_zone
    
        demand_df['zone'] = joined_gdf['z1']
        demand_df.to_csv(output_csv, index=False)



assign_zones(snakemake.input[0], snakemake.input[1], snakemake.output[0])

# Scenario 2
assign_zones(snakemake.input[2], snakemake.input[1], snakemake.output[1])