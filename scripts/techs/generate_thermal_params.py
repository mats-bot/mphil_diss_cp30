import pandas as pd
import yaml

df = pd.read_excel(snakemake.input[0])

df = df[~df["params"].isin(["start_time", "min_uptime"])]

df.to_csv(snakemake.output[0], index=False)


thermal_params = {
    "thermal_params": {
        "data": "data/processed/techs/thermal_params.csv",
        "rows": "parameters",
        "columns": "techs",
        'rename_dims': {
            'params': 'parameters'
        },
    }
}

with open(snakemake.output[1], "w") as f:
    yaml.dump({"data_tables": thermal_params}, f, sort_keys=False)