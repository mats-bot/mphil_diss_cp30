import pandas as pd
import pyproj
import numpy as np


REPD_df = pd.read_excel(snakemake.input[0],sheet_name="REPD", header=0)
ref_coords_df  = pd.read_excel(snakemake.input[0])

tech_types = REPD_df["Technology Type"].dropna().unique()
unique_df = pd.DataFrame(tech_types, columns=["Technology Type"])


# offshore_df = REPD_df[REPD_df["Technology Type"] == "Wind Offshore"]
# capacity_by_status = offshore_df.groupby("Development Status (short)")["Installed Capacity (MWelec)"].sum()
# print(capacity_by_status)

# Remove NI
REPD_df = REPD_df[REPD_df["Country"] != "Northern Ireland"]

# Keep only projects that could be built by 2030
allowed_statuses = ["Operational", "Under Construction", "Awaiting Construction", "Application Submitted"]
REPD_df = REPD_df[REPD_df["Development Status (short)"].isin(allowed_statuses)]


# Convert coordinates from BNG to WGS84
transformer = pyproj.Transformer.from_crs("epsg:27700", "epsg:4326", always_xy=True)

valid_coords = REPD_df[["X-coordinate", "Y-coordinate"]].dropna()
lons, lats = transformer.transform(valid_coords["X-coordinate"].values, valid_coords["Y-coordinate"].values)

REPD_df.loc[valid_coords.index, "Longitude"] = lons
REPD_df.loc[valid_coords.index, "Latitude"] = lats



# Handle the 20 datapoints with no coordinates, using reference coordinates from each zone
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


missing_x = REPD_df[REPD_df["X-coordinate"].isna()]
print(f"Rows missing X-coordinate: {len(missing_x)}")
print(missing_x[["Technology Type", "Installed Capacity (MWelec)", "Site Name", "Country", "Region", "County", "Post Code"]])



# keep only necessary columns

cols = ["Site Name", "Technology Type", "Installed Capacity (MWelec)","Development Status (short)", "Longitude", "Latitude"]
REPD_df = REPD_df[cols]


# Assign renewables to techs of CP30 output.


REPD_df.to_csv(snakemake.output[0])
