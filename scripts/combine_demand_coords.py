import pandas as pd

demand_df = pd.read_csv(snakemake.input["demand"])
coords_df = pd.read_csv(snakemake.input["coords"])
  
demand_key = "GSP"
coords_key = "GSP ID"   

# Used to check which ones missing from GSP coordinates list. See final comment.
unmatched_ids = set(demand_df[demand_key]) - set(coords_df[coords_key])
# if unmatched_ids:
#     print(f"\n Found {len(unmatched_ids)} unmatched GSP IDs in demand data (not found in coordinates):\n")
#     for gsp_id in sorted(unmatched_ids):
#         print(f"  - {gsp_id}")
# else:
#     print("All GSP IDs in demand data have corresponding entries in coordinates.")

merged_df = pd.merge(
    demand_df,
    coords_df,
    how="left",
    left_on=demand_key,
    right_on=coords_key
)

# Remove empty columns
merged_df = merged_df.loc[:, ~merged_df.columns.str.contains('^Unnamed')]
merged_df = merged_df.dropna(axis=1, how='all')

merged_df.to_csv(snakemake.output[0], index=False)

# After this still missing ~70 GSPs, to be added by using ETYS2017 FLOP zones available in Appendix
# A. Need to manually check minor FLOP zone on FES regional 'MAIN DATA' sheet.