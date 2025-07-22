import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
from shapely.ops import nearest_points

offshore_df = pd.read_csv(snakemake.input[0])
zones_gdf = gpd.read_file(snakemake.input[1])

# Adding missing values in data (mainly coordinates)
missing_data = offshore_df[
    (offshore_df['Latitude'].isnull() | (offshore_df['Latitude'] == 0)) |
    (offshore_df['Longitude'].isnull() | (offshore_df['Longitude'] == 0)) |
    (offshore_df['Installed Capacity (MWelec)'].isnull() | (offshore_df['Installed Capacity (MWelec)'] == 0))
]

print(missing_data[['Site Name', 'Latitude', 'Longitude', 'Installed Capacity (MWelec)', 'Development Status (short)']])

manual_data = {
    "North Falls Offshore Wind Farm": { # https://www.northfallsoffshore.com/project-description/
        "Latitude": 51.872299047018764, # https://maps.app.goo.gl/vJMx4GqocF7Yuw8bA
        "Longitude": 1.9286435026614093,
        "Installed Capacity (MW)": 504
   },
   "Morgan Offshore Wind Farm": { # https://infrastructure.planninginspectorate.gov.uk/projects/north-west/morgan-offshore-wind-project-generation-assets/
        "Latitude": 53.911474,
        "Longitude": -3.969241
    },
    "Morecombe Offshore Windfarm": { # https://web.archive.org/web/20240126062242/https://morecambeoffshorewind.com/
        "Latitude": 53.9221,
        "Longitude": -3.7088
    },
    "Dogger Bank South East & South West": { # https://geohack.toolforge.org/geohack.php?pagename=Dogger_Bank_Wind_Farm&params=54_45_N_1_55_E_region:GB_type:landmark
        "Latitude": 54.75,
        "Longitude": 1.916667
    },
    "Ossian": { # https://www.gem.wiki/Ossian_Offshore_wind_farm
        "Latitude": 56.5150,
        "Longitude": -1.3787 
    },
    "Muir Mhòr": { # https://www.gem.wiki/Muir_Mh%C3%B2r_Floating_Offshore_wind_farm
        "Latitude": 57.4189,
        "Longitude": -0.5178
    },
    "Cenos Offshore Wind Farm": { # https://www.gem.wiki/Cenos_floating_wind_farm 
        "Latitude": 57.1348,
        "Longitude": 0.6870 
    },
}

for project, values in manual_data.items():
    idx = offshore_df[offshore_df['Site Name'] == project].index
    if not idx.empty:
        for field, val in values.items():
            offshore_df.loc[idx, field] = val


# Assign to closest tzone and add distance (for transmission calculations)
offshore_df['geometry'] = offshore_df.apply(lambda row: Point(row['Longitude'], row['Latitude']), axis=1)
offshore_gdf = gpd.GeoDataFrame(offshore_df, geometry='geometry', crs="EPSG:4326")

# convert to get distances in m
offshore_gdf = offshore_gdf.to_crs(epsg=3857)
zones_gdf = zones_gdf.to_crs(epsg=3857)

def get_nearest_zone(project_point, zones):
    nearest_geom = zones.geometry.distance(project_point).sort_values().index[0]
    zone_row = zones.loc[nearest_geom]
    distance_km = project_point.distance(zone_row.geometry) / 1000  # meters to km
    return pd.Series([zone_row['z1'], distance_km])

offshore_gdf[['tzone', 'distance_km']] = offshore_gdf['geometry'].apply(lambda geom: get_nearest_zone(geom, zones_gdf))
offshore_df = pd.DataFrame(offshore_gdf.drop(columns='geometry'))



offshore_df.to_csv(snakemake.output[0],index=False)
offshore_df.to_csv(snakemake.output[1], index=False)
