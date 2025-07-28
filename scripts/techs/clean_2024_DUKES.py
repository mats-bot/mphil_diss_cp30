import pandas as pd
import geopandas as gpd
import pyproj

DUKES_df = pd.read_excel(snakemake.input[0], sheet_name="5.11 Full list", header=5)
ref_coords_df  = pd.read_excel(snakemake.input[1], header=0)
zones_gdf = gpd.read_file(snakemake.input[2])

# Keep only fossil fuel and primary fuel and sort by type
DUKES_df = DUKES_df[DUKES_df["Technology"].isin(["Fossil Fuel", "Nuclear"])]

# remove unnecessary colums
columns = ["Site Name", "Type", "CHP", "Primary Fuel", "InstalledCapacity (MW)", "Postcode", "X-Coordinate", "Y-Coordinate"]
DUKES_df = DUKES_df[columns]





# Convert coordinates from BNG to WGS84
transformer = pyproj.Transformer.from_crs("epsg:27700", "epsg:4326", always_xy=True)

valid_coords = DUKES_df[["X-Coordinate", "Y-Coordinate"]].dropna()
valid_coords = valid_coords.drop(DUKES_df[DUKES_df["Site Name"] == "Grain (OCGT)"].index) # Drop this 
# one manually since coords are (0,0)
lons, lats = transformer.transform(valid_coords["X-Coordinate"].values, valid_coords["Y-Coordinate"].values)

DUKES_df.loc[valid_coords.index, "Longitude"] = lons
DUKES_df.loc[valid_coords.index, "Latitude"] = lats

# give entries with no coordinates tzone reference coordinates
def assign_coords(row, zone_dict, zone_coords):
    site_name = row["Site Name"]
    if site_name in zone_dict:
        tzone = zone_dict[site_name]
        match = zone_coords[zone_coords["Tzone name"] == tzone]
        if not match.empty:
            row["Latitude"] = match["Reference Latitude"].values[0]
            row["Longitude"] = match["Reference Longitude"].values[0]
    return row

onshore_dict = {
    "Spalding Expansion OCGT ": "z17", 
    "Croydon": "z17",
    "Derby": "z10",
    "Exeter": "z15",
    "Heartlands": "z10",
    "Peterborough Power Station PPL2 Gas Engines": "z11",
    "Thornhill": "z8",
    "Viking": "z16",
    "Kings Lynn": "z11",
    "Cheshire": "z7",
    "Hythe": "z17",
    "Grimsby": "z8",
    "Keadby 2": "z8",
    "Grain (OCGT)": "z17"
    }

DUKES_df = DUKES_df.apply(lambda row: assign_coords(row, onshore_dict, ref_coords_df), axis=1)



# drop the coordinates in BNG
DUKES_df = DUKES_df.drop(columns=["X-Coordinate", "Y-Coordinate"])


# maps to link with tech names in cost processing, see aggregation assumptions in mapping file
# may need to modify to incl chp or undefined techs
TECH_MAPPING = {
    ("CCGT", "Natural Gas"): "Gas_CCGT",
    ("CCGT", "Sour Gas"): "Gas_CCGT", 
    ("Conventional steam", "Coal"): "Coal",
    ("Conventional steam", "Natural Gas"): "Gas_CCGT",
    ("Single cycle", "Natural Gas"): "Gas_OCGT",
    ("Single cycle", "Diesel/Gas Oil"): "Diesel",
    ("AGR", None): "Nuclear",
    ("PWR", None): "Nuclear",
}
    
def map_cp30_technology(row):
    type_ = row["Type"]
    fuel = row.get("Primary Fuel", None)
    return TECH_MAPPING.get((type_, fuel), TECH_MAPPING.get((type_, None), "other"))

    
DUKES_df["CP30 technology"] = DUKES_df.apply(map_cp30_technology, axis=1)


# Print to see current capacity
cp30_summary = DUKES_df.groupby("CP30 technology")["InstalledCapacity (MW)"].sum().reset_index()
print(cp30_summary)

# remove useles cols
DUKES_df = DUKES_df.drop(columns=["Type", "Primary Fuel"])


# Map the power generators to the tzones
DUKES_gdf = gpd.GeoDataFrame(
    DUKES_df,
    geometry=gpd.points_from_xy(DUKES_df["Longitude"], DUKES_df["Latitude"]),
    crs="EPSG:4326"
)
joined_gdf = gpd.sjoin(DUKES_gdf, zones_gdf[["z1", "geometry"]], how="left", predicate="within")


# Handle points outside of the zones 
unassigned_mask = joined_gdf['z1'].isna()
unassigned_gdf = DUKES_gdf[unassigned_mask].copy()

if not unassigned_gdf.empty:  
    # Iterates over unassigned points to find closest zone   
    for idx, point in unassigned_gdf.geometry.items():
        nearest_zone = None
        min_distance = 5
        

        for zone_idx, zone in zones_gdf.iterrows():
            distance = point.distance(zone.geometry)
            if distance < min_distance:
                min_distance = distance
                nearest_zone = zone['z1']
        joined_gdf.loc[idx, 'z1'] = nearest_zone



DUKES_df['zone'] = joined_gdf['z1']

# Check all power stations assgined to a zone
# unassigned_final = DUKES_df[DUKES_df["zone"].isna()]
# print("Unassigned entries after final zone mapping:")
# print(unassigned_final[["Site Name", "Latitude", "Longitude", "CP30 technology", "InstalledCapacity (MW)"]])

# Aggregate data by zone
DUKES_df = (
    DUKES_df
    .groupby(["zone", "CP30 technology"])["InstalledCapacity (MW)"]
    .sum()
    .reset_index()
    .sort_values(by=["zone", "InstalledCapacity (MW)"])
)


DUKES_df.to_csv(snakemake.output[0], index=None)