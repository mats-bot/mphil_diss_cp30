import pandas as pd
import yaml

df = pd.read_csv(snakemake.input[0], index_col=0)


