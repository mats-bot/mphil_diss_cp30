import pandas as pd

demand_df = pd.read_csv(snakemake.input["demand"])
coords_df = pd.read_csv(snakemake.input["coords"])
unmatched_coords_df = pd.read_excel(snakemake.input["unmatched_GSP_coords"])


demand_key = "GSP"
coords_key = "GSP ID"

# merge coordinates of GSPs for which we have IDs
merged_df = pd.merge(
    demand_df, coords_df, how="left", left_on=demand_key, right_on=coords_key
)

# Used to check which ones missing from GSP coordinates list. Created document using this where 
# manually matched coordinates to GSP_IDs in function of FES24 regional breakdown "MAIN DATA" sheet
# giving location information (e.g. name of a town).
unmatched_ids = set(demand_df[demand_key]) - set(coords_df[coords_key])
if unmatched_ids:
    print(f"\n Found {len(unmatched_ids)} unmatched GSP IDs in demand data (not found in coordinates):\n")

for i, row in merged_df.iterrows():
    if pd.isna(row["Latitude"]) or pd.isna(row["Longitude"]):
        gsp = row["GSP"]
        manual_match = unmatched_coords_df[unmatched_coords_df["GSP"] == gsp]
        if not manual_match.empty:
            merged_df.at[i, "Latitude"] = manual_match.iloc[0]["Latitude"]
            merged_df.at[i, "Longitude"] = manual_match.iloc[0]["Longitude"]

for idx, row in merged_df.iterrows():
    if pd.isna(row["Latitude"]) or pd.isna(row["Longitude"]):
        print(f"Missing coords at row {idx}: {row['Latitude']}, {row['Longitude']}")


# Remove empty columns
merged_df = merged_df.loc[:, ~merged_df.columns.str.contains("^Unnamed")]
merged_df = merged_df.dropna(axis=1, how="all")

merged_df = merged_df.drop(columns=["scenario", "GSP ID"])

merged_df.to_csv(snakemake.output[0], index=False)

# After this still missing ~70 GSPs, to be added by using ETYS2017 FLOP zones available in Appendix
# A. Need to manually check minor FLOP zone on FES regional 'MAIN DATA' sheet.
