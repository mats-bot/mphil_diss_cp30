import pandas as pd
import yaml

# For renewables max capacities (from REPD queue, entries with missing capacities omitted)
df3 = pd.read_csv(snakemake.input[0])

ren_flow_cap_max_table = {
    "data": "data/processed/techs/renewables_2030_queue.csv",
    "rows": "techs",
    "columns": "nodes",
    "add_dims": {"parameters": "flow_cap_max"},
}


# Set infinite flow cap for other techs so model can build as much as it wants
unconstrained_caps = {
    "gas_ccgt_new": {
        "flow_cap_max": ".inf"
    },
    "gas_ocgt_new": {
        "flow_cap_max": ".inf"
    },
    "beccs": {
        "flow_cap_max": ".inf"
    },
    "gas_ccs": {
        "flow_cap_max": ".inf"
    },
}


output_data = {
    "data_tables": {
        "renewables_flow_cap_max": ren_flow_cap_max_table
    }
}

with open(snakemake.output[0], "w") as f:
    yaml.dump(output_data, f, sort_keys=False)