import pandas as pd

demand_df = pd.read_csv(snakemake.input["demand"])
coords_df = pd.read_csv(snakemake.input["coords"])
  
demand_key = "GSP"
coords_key = "GSP ID"   

unmatched_ids = set(demand_df[demand_key]) - set(coords_df[coords_key])
print(f'\nNUMBER OF UNMATCHED GSPS: {len(unmatched_ids)}')
# if unmatched_ids:
#     print(f"WARNING: {len(unmatched_ids)} GSP IDs in demand data not found in coordinates:")
#     for uid in unmatched_ids:
#         print(f"  - {uid}")

merged_df = pd.merge(
    demand_df,
    coords_df,
    how="left",
    left_on=demand_key,
    right_on=coords_key
)


merged_df.to_csv(snakemake.output[0], index=False)