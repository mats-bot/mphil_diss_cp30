import pandas as pd

power_df = pd.read_excel(snakemake.input[0], sheet_name="ES1", header=9)

generation_types = power_df["SubType"].dropna().unique()
unique_df = pd.DataFrame(generation_types, columns=["SubType"])


unique_df.to_csv(snakemake.output[0], index=False)


