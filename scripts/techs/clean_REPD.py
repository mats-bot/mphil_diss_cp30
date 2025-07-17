import pandas as pd
import pyproj
import numpy as np


REPD_df = pd.read_excel(snakemake.input[0],sheet_name="REPD", header=0)

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

# def convert_coords(x, y):
#     if pd.notna(x) and pd.notna(y):
#         lon, lat = transformer.transform(x, y)
#         return pd.Series([lon, lat])
#     else:
#         return pd.Series([np.nan, np.nan])
    
# REPD_df[["Longitude", "Latitude"]] = REPD_df.apply(lambda row: convert_coords(row["X-coordinate"], row["Y-coordinate"]), axis=1)

valid_coords = REPD_df[["X-coordinate", "Y-coordinate"]].dropna()
lons, lats = transformer.transform(valid_coords["X-coordinate"].values, valid_coords["Y-coordinate"].values)

# Assign results back, filling NaNs where coordinates were missing
REPD_df.loc[valid_coords.index, "Longitude"] = lons
REPD_df.loc[valid_coords.index, "Latitude"] = lats


print(transformer.transform(545000, 180000))

# Add coordinates to the 20 datapoints with missing ones



missing_x = REPD_df[REPD_df["X-coordinate"].isna()]
print(f"Rows missing X-coordinate: {len(missing_x)}")




# keep only necessary columns

cols = ["Site Name", "Technology Type", "Installed Capacity (MWelec)","Development Status (short)", "Longitude", "Latitude"]
REPD_df = REPD_df[cols]


# Assign renewables to techs of CP30 output.


REPD_df.to_csv(snakemake.output[0])
