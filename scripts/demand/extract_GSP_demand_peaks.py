import pandas as pd

df = pd.read_excel(snakemake.input[0], sheet_name="Active", header=None)

# skip empty rows
df = df.iloc[2:]

# keep headers for further processing
df.columns = df.iloc[0]
df = df[1:]

# remove non-HT, data after 2030.
df = df[df.iloc[:, 0].astype(str).str.contains("HT")]
df = df[pd.to_numeric(df.iloc[:, 6], errors="coerce") < 31]
# check if empty
df = df[df.iloc[:, 2].notna() & (df.iloc[:, 2].astype(str).str.strip() != "")]

# remove GSPs where no demand estimated before 2030 (often no location info)
df["DemandPk"] = pd.to_numeric(df["DemandPk"], errors="coerce")
non_zero_gsps = df.groupby(df.columns[1])["DemandPk"].transform("sum") > 0
df = df[non_zero_gsps]

df.to_csv(snakemake.output[0], index=False)