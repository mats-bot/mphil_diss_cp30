import pandas as pd
import pyproj
import geopandas as gpd

REPD_df = pd.read_excel(snakemake.input[0],sheet_name="REPD", header=0)
ref_coords_df  = pd.read_excel(snakemake.input[1], header=0)
zones_gdf = gpd.read_file(snakemake.input[2])

tech_types = REPD_df["Technology Type"].dropna().unique()
unique_df = pd.DataFrame(tech_types, columns=["Technology Type"])


# offshore_df = REPD_df[REPD_df["Technology Type"] == "Wind Offshore"]
# capacity_by_status = offshore_df.groupby("Development Status (short)")["Installed Capacity (MWelec)"].sum()
# print(capacity_by_status)

# Remove NI, Jersey, Isle of Man
REPD_df = REPD_df[REPD_df["Country"] != "Northern Ireland"]
REPD_df = REPD_df[~REPD_df["County"].isin(["Jersey", "Isle of Man"])]

# Keep only projects that could be built by 2030
allowed_statuses = ["Operational", "Under Construction", "Awaiting Construction", "Application Submitted"]
REPD_df = REPD_df[REPD_df["Development Status (short)"].isin(allowed_statuses)]


# Convert coordinates from BNG to WGS84
transformer = pyproj.Transformer.from_crs("epsg:27700", "epsg:4326", always_xy=True)

valid_coords = REPD_df[["X-coordinate", "Y-coordinate"]].dropna()
lons, lats = transformer.transform(valid_coords["X-coordinate"].values, valid_coords["Y-coordinate"].values)

REPD_df.loc[valid_coords.index, "Longitude"] = lons
REPD_df.loc[valid_coords.index, "Latitude"] = lats



# Handle the datapoints with no coordinates, using reference coordinates from each zone
# Seperate handling for offshore wind farms since need precise coordinates
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
    "Sainsbury's Stores (169 individual Stores)": "z17", # Since most stores in z17
    "First Wessex Housing Properties (multiple)": "z16",
    "Bristol City Council - Solar PV Project": "z14",
    "Balliemeanoch Pumped Hydro Project": "z3",
    "Persley Wastewater Treatment Plant - Battery Storage": "z2",
    "Mynydd Maen Solar Farm (Cil-Lonydd)": "z13",
    "Alwen Forest Project": "z7",
    "Seastar Project - Tidal Energy Farm": "z1",
    "Oceanstar Project - Tidal Energy Farm": "z1",
    "Mynydd Maen Solar Farm (Cil-Lonydd)": "z13",
    "Mynydd Maen Solar Farm (Cil-Lonydd)": "z13"
}

REPD_df = REPD_df.apply(lambda row: assign_coords(row, onshore_dict, ref_coords_df), axis=1)


# keep only necessary columns

cols = ["Site Name", "Technology Type", "CHP Enabled", "Installed Capacity (MWelec)",
        "Development Status (short)", "Longitude", "Latitude", "Planning Application Submitted", 
        "Planning Permission Granted", "Under Construction", "Operational"]
REPD_df = REPD_df[cols]


# Assign renewables to techs of CP30 output.
def map_cp30_technology(row):

    # Framework for splitting inputs into Biomass vs Biomass CHP. Currently disabled because of 
    # significant difference with CP30, check methodology.
    if row["Technology Type"] == "Biomass (co-firing)":
        if str(row["CHP Enabled"]).strip().lower() == "yes":
            return "Biomass CHP"
        elif str(row["CHP Enabled"]).strip().lower() == "no":
            return "Biomass"
        
    elif row["Technology Type"] == "Biomass (dedicated)":
        if str(row["CHP Enabled"]).strip().lower() == "yes":
            return "Biomass CHP"
        elif str(row["CHP Enabled"]).strip().lower() == "no":
            return "Biomass"
        else:
            return "Biomass" # Concerns only a few samll plants and a Drax unit known not to have CHP
    
    elif row["Technology Type"] == "EfW Incineration":
        if str(row["CHP Enabled"]).strip().lower() == "yes":
            return "Waste"
        elif str(row["CHP Enabled"]).strip().lower() == "no":
            return "Waste"
        else:
            # print(f'\nIssue with waste')
            # print(row)
            return "Waste"
    
    elif row["Technology Type"] == "Large Hydro":
        return "Hydro"
    
    elif row["Technology Type"] == "Small Hydro":
        return "Hydro"
    
    elif row["Technology Type"] == "Battery":
        return "Battery"
    
    elif row["Technology Type"] == "Solar Photovoltaics":
        return "Solar PV"

    elif row["Technology Type"] == "Wind Offshore":
        return "Offshore Wind"
    
    elif row["Technology Type"] == "Wind Onshore":
        return "Onshore Wind"
    
    elif row["Technology Type"] == "Pumped Storage Hydroelectricity":
        return "Pumped hydro (LDES)"
    
    elif row["Technology Type"] == "Liquid Air Energy Storage":
        return "LAES (LDES)"
    
    elif row["Technology Type"] == "Compressed Air Energy Storage":
        return "CAES (LDES)"
    
    elif row["Technology Type"] == "Hydrogen":
        if str(row["CHP Enabled"]).strip().lower() == "yes":
            return "Hydrogen"
        elif str(row["CHP Enabled"]).strip().lower() == "no":
            return "Hydrogen"
        else:
            return "Waste"
    else:
        return "Other Renewables"
    

REPD_df["CP30 technology"] = REPD_df.apply(map_cp30_technology, axis=1)
REPD_df = REPD_df.drop(columns=["Technology Type"])


# Find 2023 total capacities for comparison with CP30
REPD_df["Operational"] = pd.to_datetime(REPD_df["Operational"], errors="coerce")

operational_2023_df = REPD_df[
    (REPD_df["Development Status (short)"] == "Operational") &
    (REPD_df["Operational"] < "2024-01-01")
]

capacity_by_cp30 = operational_2023_df.groupby("CP30 technology")["Installed Capacity (MWelec)"].sum()
print(capacity_by_cp30)

# Remove offshore wind since this needs to be handled seperately
offshore_wind_df = operational_2023_df[operational_2023_df["CP30 technology"] == "Offshore Wind"]
offshore_wind_df.to_csv(snakemake.output[0], index=False)
total_offshore_capacity = offshore_wind_df["Installed Capacity (MWelec)"].sum()
print(f"Total installed capacity of offshore wind: {total_offshore_capacity:.2f} MW")

operational_2023_df = operational_2023_df[operational_2023_df["CP30 technology"] != "Offshore Wind"]



# Map the power generators to the tzones
operational_2023_gdf = gpd.GeoDataFrame(
    operational_2023_df,
    geometry=gpd.points_from_xy(operational_2023_df["Longitude"], operational_2023_df["Latitude"]),
    crs="EPSG:4326"  
)
joined_gdf = gpd.sjoin(operational_2023_gdf, zones_gdf[["z1", "geometry"]], how="left", predicate="within")


# Handle points outside of the zones 
unassigned_mask = joined_gdf['z1'].isna()
unassigned_gdf = operational_2023_gdf[unassigned_mask].copy()

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

operational_2023_df['zone'] = joined_gdf['z1']

# Check all entries assigned a zone
unassigned_final = operational_2023_df[operational_2023_df["zone"].isna()]
print("Unassigned entries after final zone mapping:")
print(unassigned_final[["Site Name", "Latitude", "Longitude", "CP30 technology", "Installed Capacity (MWelec)"]])


# Handle entries where capacity is empty - 40 entries out of ~2000
operational_2023_df["Installed Capacity (MWelec)"] = pd.to_numeric(
    operational_2023_df["Installed Capacity (MWelec)"], errors="coerce"
)

# Check these entries 
missing_capacity = operational_2023_df[operational_2023_df["Installed Capacity (MWelec)"].isna()]
print("Entries with missing Installed Capacity:")
print(missing_capacity[["Site Name", "CP30 technology", "zone", "Latitude", "Longitude"]])
print(f"Total entries with missing capacity: {len(missing_capacity)}")
print(f"Entries with missing capacity and technology 'Battery': {len(missing_capacity[missing_capacity['CP30 technology'] == 'Battery'])}")



# Aggregate data by zone
operational_2023_df = (
    operational_2023_df
    .groupby(["zone", "CP30 technology"])["Installed Capacity (MWelec)"]
    .sum()
    .reset_index()
    .sort_values(by=["zone", "Installed Capacity (MWelec)"])
)


operational_2023_df.to_csv(snakemake.output[1])


