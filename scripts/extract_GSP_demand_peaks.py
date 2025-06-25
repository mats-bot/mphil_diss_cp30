import pandas as pd


df = pd.read_excel(snakemake.input[0], sheet_name="Active", header=None)

# skip empty rows
df = df.iloc[2:]

# keep headers if needed later (otherwise can remove)
df.columns = df.iloc[0]
df = df[1:]

# remove non-HT, after 2030.
df = df[df.iloc[:, 0].astype(str).str.contains("HT")]
df = df[pd.to_numeric(df.iloc[:, 6], errors='coerce') < 31]
df = df[df.iloc[:, 2].notna() & (df.iloc[:, 2].astype(str).str.strip() != "")] # check if empty

df.to_csv(snakemake.output[0], index=False)