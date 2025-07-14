import pandas as pd
import geopandas as gpd
from geopy.distance import geodesic


df = pd.read_csv(snakemake.input[0])
centroids = gpd.read_file(snakemake.input[1])

locations_dict = {
    row["z1"]: {
        "coordinates": {
            "lat": float(row.geometry.y),
            "lon": float(row.geometry.x)
        }
    }
    for _, row in centroids.iterrows()
}

# Value taken from national grid, mean figure chosen for AC transmission since GB almost all AC
# https://www.nationalgrid.com/sites/default/files/documents/13784-High%20Voltage%20Direct%20Current%20Electricity%20%E2%80%93%20technical%20information.pdf
loss_per_100km = 0.0065


def calculate_loss(zone1, zone2, loss_rate_per_100km):
    """Returns loss fraction for each onshore link"""
    z1 = zone1.strip()
    z2 = zone2.strip()

    coords1 = (locations_dict[z1]["coordinates"]["lat"], 
               locations_dict[z1]["coordinates"]["lon"])
    coords2 = (locations_dict[z2]["coordinates"]["lat"], 
                   locations_dict[z2]["coordinates"]["lon"])
        
    distance_km = geodesic(coords1, coords2).km
    return (distance_km / 100) * loss_rate_per_100km 

df["loss"] = df.apply(
    lambda row: calculate_loss(
        row["Zone 1"],
        row["Zone 2"],
        loss_per_100km
    ),
    axis=1
)


df.to_csv(snakemake.output[0])