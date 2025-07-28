import pandas as pd
import yaml

df = pd.read_csv(snakemake.input[0], index_col=0)

renewable_techs = ["Onshore_Wind", "Solar_PV"]



techs_yaml = {}
