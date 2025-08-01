import pandas as pd
import yaml
import numpy as np

df = pd.read_csv(snakemake.input[0])

capacity_df = pd.read_csv(snakemake.input[1])
zones = sorted(capacity_df["zone"].unique())

storage_techs = {"battery", "pumped_hydro", "caes", "laes"}

disallow_investment = {"pumped_hydro"}

tech_zone_capacity = {}
# for tech in capacity_df["CP30 technology"].unique():
for tech in storage_techs:
    tech_zone_capacity[tech.lower()] = [0.0] * len(zones)

# TODO: turn this into a per-node definition
# TODO: add storage caps
for _, row in capacity_df.iterrows():
    zone = row["zone"]
    tech = row["CP30 technology"].lower()  # lowercased to match keys
    capacity = float(row["Installed Capacity (MWelec)"])
    if tech in tech_zone_capacity and zone in zones:
        idx = zones.index(zone)
        tech_zone_capacity[tech][idx] = capacity

capacities = pd.DataFrame.from_dict(
    tech_zone_capacity, orient="index", columns=zones
)
capacities.to_csv(snakemake.output[1])




# discharge_time = {
#     "battery_existing": 4.0,
#     "pumped_hydro_existing": 8,
#     "caes_existing": 6,
#     "laes_existing": 14,
# }       

# zonal def
# nodes = {}
# for zone in capacity_df.columns:
#     nodes[zone] = {'techs': {}}
#     for tech in storage_techs:
#         flow_cap = float(capacity_df.loc[(capacity_df['zone'] == zone) & (capacity_df['CP30 technology'].lower() == tech), 'Installed Capacity (MWelec)'].values[0])
#         storage_cap = flow_cap * storage_techs[tech]
#         flow_cap_min = {"data": flow_cap, "dims": ["carriers"], "index": zones}
#         storage_cap_min = {"data": storage_cap, "dims": ["carriers"], "index": zones}

#         nodes[zone]['techs'][tech] = {
#                 tech: {
#                     'flow_cap_min': flow_cap,
#                     'storage_cap_min': storage_cap
#             }
#         }
storage_techs = {}

for _, row in df.iterrows():
    tech = row["technology"]

    # Grouping by cp30_category
    category = "battery" if tech == "battery" else "ldes"

    flow_cap_data = tech_zone_capacity.get(tech.lower(), [0.0] * len(zones))
    # Calculate storage_cap_max = flow_cap_max * (1 / charge_rate)
    charge_rate = float(row["charge_rate"])
    if charge_rate == 0:
        storage_cap_data = [0.0] * len(zones)
    else:
        storage_cap_data = [fc / charge_rate for fc in flow_cap_data]


    def make_cost_dict(value):
        return {"data": value, "index": "monetary", "dims": ["costs"]}

    # Parent tech with shared constraints and base category
    storage_techs[tech] = {
        "category": category,
        "cp30_category": "storage",
        "base_tech": "storage",
        "name": tech,
        "carrier_in": "electricity",
        "carrier_out": "electricity",
        "flow_out_eff": float(row["energy_eff"]),
        "storage_loss": float(row["storage_loss"]),
        "flow_cap_per_storage_cap_max": float(row["charge_rate"]),
        "lifetime": int(row["lifetime"]),
        "cost_om_annual": make_cost_dict(
            float(row["om_annual"]) if pd.notna(row["om_annual"]) else 0
        ),
        "cost_flow_out": make_cost_dict(
            float(row["om_prod"]) if pd.notna(row["om_prod"]) else 0
        ),
        "cost_flow_in": make_cost_dict(
            float(row["om_con"]) if pd.notna(row["om_con"]) else 0
        ),
        "cost_flow_cap": make_cost_dict(
            float(row["cost_energy_cap"]) if pd.notna(row["cost_energy_cap"]) else 0
        ),
        "cost_storage_cap": make_cost_dict(
            float(row["cost_storage_cap"]) if pd.notna(row["cost_storage_cap"]) else 0
        ),
        # "flow_cap_min": float(row["energy_cap_equals"]),
        # "flow_cap_max": float(row["energy_cap_equals"]),
        # "storage_cap_min": float(row["storage_cap_equals"]),
        # "storage_cap_max": float(row["storage_cap_equals"]),
    }


# data_table = {
#     'flow_cap_min': {
#         'data': 'data/processed/techs/storage_capacities.csv',
#         'rows': 'techs',
#         'columns': 'nodes',
#         'add_dims': {
#             'parameters': 'flow_cap_min'
#             }
#         }
#     }

with open(snakemake.output[0], "w") as f:
    yaml.dump({"techs": storage_techs}, f, sort_keys=False)
