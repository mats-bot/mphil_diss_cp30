import pandas as pd

df = pd.read_excel(snakemake.input[0], sheet_name="GSP info", skiprows=4)

# skip whitespace
df.columns = df.columns.str.strip()

df = df[["GSP ID", "Latitude", "Longitude"]]

df.to_csv(snakemake.output[0], index=False)