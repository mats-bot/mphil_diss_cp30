import pandas as pd
import numpy as np
import yaml

# For renewables max capacities (from REPD queue, entries with missing capacities omitted)
df3 = pd.read_csv(snakemake.input[0])

ren_flow_cap_max_table = {
    "data": "data/processed/techs/renewables_2030_queue.csv",
    "rows": "techs",
    "columns": "nodes",
    "add_dims": {
        "parameters": "flow_cap_max"
        },
}


# Set infinite flow cap for other techs so model can build as much as it wants
nodes = df3.columns.tolist()
nodes.remove('techs')

unconstrained_techs = ["gas_ccgt_new", "gas_ocgt_new", "beccs", "gas_ccs"]

unconstrained_df = pd.DataFrame(
    data=np.inf,
    index=unconstrained_techs,
    columns=nodes
)

unconstrained_df.index.name = "techs"
unconstrained_df.columns.name = "nodes"
unconstrained_df.to_csv(snakemake.output[1])

unconstrained_caps_table = {
    "data": "data/processed/techs/unconstrained_techs_2030.csv",
    "rows": "techs",
    "columns": "nodes",
    "add_dims": {
        "parameters": ["flow_cap_max"]
    },
}

output_data = {
    "data_tables": {
        "renewables_flow_cap_max": ren_flow_cap_max_table,
        "unconstrained_flow_cap_max": unconstrained_caps_table,
    }
}

with open(snakemake.output[0], "w") as f:
    yaml.dump(output_data, f, sort_keys=False)