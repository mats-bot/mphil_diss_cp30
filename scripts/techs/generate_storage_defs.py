import pandas as pd
import yaml

df = pd.read_csv(snakemake.input[0], index_col=0)

storage_techs = {}

for tech in df.index:
    row = df.loc[tech]

    # Default to 0 if missing (NaN)
    def get(val):
        return float(val) if pd.notna(val) else 0

    storage_techs[tech.lower()] = {
        "category": "storage",
        "cp30_category": "thermal",

        "essentials": {
            "name": tech,
            "carrier_in": "electricity",
            "carrier_out": "electricity",
            "storage": True
        },
        "constraints": {
            "energy_eff": float(row["energy_eff"]),
            "storage_loss": float(row["storage_loss"]),
            "charge_rate": float(row["charge_rate"]),
            "energy_cap_equals": float(row["energy_cap_equals"]),
            "storage_cap_equals": float(row["storage_cap_equals"]),
            "lifetime": int(row["lifetime"])
        },
        "costs": {
            "cost_energy_cap": get(row["cost_energy_cap"]),
            "cost_storage_cap": get(row["cost_storage_cap"]),
            "om_annual": get(row["om_annual"]),
            "om_prod": get(row["om_prod"]),
            "om_con": get(row["om_con"])
        }
    }

with open(snakemake.output[0], "w") as f:
    yaml.dump(storage_techs, f, sort_keys=False)