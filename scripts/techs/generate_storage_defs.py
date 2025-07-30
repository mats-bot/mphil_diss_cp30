import pandas as pd
import yaml

df = pd.read_csv(snakemake.input[0])

capacity_df = pd.read_csv(snakemake.input[1])
zones = sorted(capacity_df['zone'].unique())

storage_techs = {}

disallow_investment = {"pumped_hydro"}

tech_zone_capacity = {}
for tech in capacity_df['CP30 technology'].unique():
    tech_zone_capacity[tech.lower()] = [0.0] * len(zones)

for _, row in capacity_df.iterrows():
    zone = row['zone']
    tech = row['CP30 technology'].lower()  # lowercased to match keys
    capacity = float(row['Installed Capacity (MWelec)'])
    if tech in tech_zone_capacity and zone in zones:
        idx = zones.index(zone)
        tech_zone_capacity[tech][idx] = capacity


for _, row in df.iterrows():
    tech = row["technology"]

    # Grouping by cp30_category
    category = "battery" if tech == "battery" else "ldes"

    # Parent tech with shared constraints and base category
    storage_techs[tech] = {
        "category": category,
        "cp30_category": "storage",
        "base_tech": "storage",
        "name": tech,
        "carrier_in": "electricity",
        "carrier_out": "electricity",
        "energy_eff": float(row["energy_eff"]),
        "storage_loss": float(row["storage_loss"]),
        "charge_rate": float(row["charge_rate"]),
        "lifetime": int(row["lifetime"]),
    }


    flow_cap_data = tech_zone_capacity.get(tech.lower(), [0.0] * len(zones))
    # Calculate storage_cap_max = flow_cap_max * (1 / charge_rate)
    charge_rate = float(row["charge_rate"])
    if charge_rate == 0:
        # Avoid division by zero; fallback to zeros or raise error as needed
        storage_cap_data = [0.0] * len(zones)
    else:
        storage_cap_data = [fc / charge_rate for fc in flow_cap_data]

    flow_cap_max = {
        "data": flow_cap_data,
        "dims": ["carriers"],
        "index": zones,
    }
    storage_cap_max = {
        "data": storage_cap_data,
        "dims": ["carriers"],
        "index": zones,
    }

    # Existing version (fixed capacity, no capex)
    storage_techs[f"{tech}_existing"] = {
        "parent": tech,
        "base_tech": "storage",
        "cost_energy_cap": 0,
        "cost_storage_cap": 0,
        "om_annual": float(row["om_annual"]) if pd.notna(row["om_annual"]) else 0,
        "om_prod": float(row["om_prod"]) if pd.notna(row["om_prod"]) else 0,
        "om_con": float(row["om_con"]) if pd.notna(row["om_con"]) else 0,
        "flow_cap_max": flow_cap_max,
        "storage_cap_max": storage_cap_max,
    }

    # New version (investable only if allowed)
    if tech not in disallow_investment:
        storage_techs[f"{tech}_new"] = {
            "parent": tech,
            "cost_energy_cap": float(row["cost_energy_cap"]) if pd.notna(row["cost_energy_cap"]) else 0,
            "cost_storage_cap": float(row["cost_storage_cap"]) if pd.notna(row["cost_storage_cap"]) else 0,
            "om_annual": float(row["om_annual"]) if pd.notna(row["om_annual"]) else 0,
            "om_prod": float(row["om_prod"]) if pd.notna(row["om_prod"]) else 0,
            "om_con": float(row["om_con"]) if pd.notna(row["om_con"]) else 0,
            "flow_cap_min": float(row["energy_cap_equals"]),
            "flow_cap_max": float(row["energy_cap_equals"]),
            "storage_cap_min": float(row["storage_cap_equals"]),
            "storage_cap_max": float(row["storage_cap_equals"]),
            "flow_cap_per_unit": float(row["energy_cap_equals"]),
            "storage_cap_per_unit": float(row["storage_cap_equals"]),
            "integer": True,
            "units_max": 1
        }

with open(snakemake.output[0], "w") as f:
    yaml.dump({"techs": storage_techs}, f, sort_keys=False)