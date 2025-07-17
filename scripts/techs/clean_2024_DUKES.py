import pandas as pd

power_df = pd.read_excel(snakemake.input[0], sheet_name="5.11 Full list", header=5)





power_df.to_csv(snakemake.output[0])