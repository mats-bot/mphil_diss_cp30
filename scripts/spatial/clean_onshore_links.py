import pandas as pd
import yaml

GTC_df = pd.read_excel(snakemake.input[0], sheet_name="FES24 Chart Data", header=0, index_col=None)
links_df = pd.read_excel(snakemake.input[1], index_col=None)

GTC_df = GTC_df[GTC_df["Scenario"] == "Holistic Transition"]
GTC_df = GTC_df[GTC_df.iloc[:, 0] == "Capability"]

GTC_cols = GTC_df.columns.tolist()
idx_2030 = GTC_cols.index(2030)
GTC_df = GTC_df.iloc[:, :idx_2030+1]

# Add column for 2023, assuming same transmission values as 2024
idx_2024 = GTC_cols.index(2024)
GTC_df.insert(loc=idx_2024, column=2023, value=GTC_df[2024])

links_df = links_df.merge(GTC_df, how="left", left_on="ETYS boundary", right_on="Boundary")
links_df.drop(columns=["Unnamed: 0", "Boundary", "Scenario"], inplace=True)

links_df.to_csv(snakemake.output[0])
   