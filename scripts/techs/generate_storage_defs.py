import pandas as pd
import yaml

df = pd.read_csv(snakemake.input[0])

storage_techs = {}

disallow_investment = {"pumped_hydro"}


for _, row in df.iterrows():
    tech = row["technology"]

    # Grouping by cp30_category
    category = "battery" if tech == "battery" else "ldes"

    # Parent tech with shared constraints and base category
    storage_techs[tech] = {
        "category": category,
        "cp30_category": "storage",
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
            "lifetime": int(row["lifetime"])
        }
    }

    # Existing version (fixed capacity, no capex)
    storage_techs[f"{tech}_existing"] = {
        "parent": tech,
        "costs": {
            "cost_energy_cap": 0,
            "cost_storage_cap": 0,
            "om_annual": float(row["om_annual"]) if pd.notna(row["om_annual"]) else 0,
            "om_prod": float(row["om_prod"]) if pd.notna(row["om_prod"]) else 0,
            "om_con": float(row["om_con"]) if pd.notna(row["om_con"]) else 0,
        },
        "constraints": {
            "energy_cap_equals": float(row["energy_cap_equals"]),
            "storage_cap_equals": float(row["storage_cap_equals"])
        }
    }

    # New version (investable only if allowed)
    if tech not in disallow_investment:
        storage_techs[f"{tech}_new"] = {
            "parent": tech,
            "costs": {
                "cost_energy_cap": float(row["cost_energy_cap"]) if pd.notna(row["cost_energy_cap"]) else 0,
                "cost_storage_cap": float(row["cost_storage_cap"]) if pd.notna(row["cost_storage_cap"]) else 0,
                "om_annual": float(row["om_annual"]) if pd.notna(row["om_annual"]) else 0,
                "om_prod": float(row["om_prod"]) if pd.notna(row["om_prod"]) else 0,
                "om_con": float(row["om_con"]) if pd.notna(row["om_con"]) else 0,
            },
            "constraints": {
                "energy_cap_per_unit": float(row["energy_cap_equals"]),
                "storage_cap_per_unit": float(row["storage_cap_equals"]),
                "integer": True,
                "units_max": 1
            }
        }

with open(snakemake.output[0], "w") as f:
    yaml.dump(storage_techs, f, sort_keys=False)