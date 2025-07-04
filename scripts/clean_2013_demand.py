import pandas as pd

df = pd.read_csv(snakemake.input[0])

df = df[["SETTLEMENT_DATE", "SETTLEMENT_PERIOD", "ND"]]
df["SETTLEMENT_DATE"] = pd.to_datetime(df["SETTLEMENT_DATE"])
df["HOUR"] = ((df["SETTLEMENT_PERIOD"] - 1) // 2) + 1
hourly_df = df.groupby(["SETTLEMENT_DATE", "HOUR"], as_index=False)["ND"].sum()
hourly_df.to_csv(snakemake.output[0], index=False)
