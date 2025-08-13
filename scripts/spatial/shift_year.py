import pandas as pd

old_year = 2013
new_year = 2030

df = pd.read_csv(snakemake.input[0], index_col=0, parse_dates=True)
new_index = df.index.map(lambda ts: ts.replace(year=new_year))
df.index = new_index
df.to_csv(snakemake.output[0])
