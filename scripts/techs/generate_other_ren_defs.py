import pandas as pd
import yaml

df = pd.read_csv(snakemake.input[0], index_col=0)
capacity_df = pd.read_csv(snakemake.input[1])
zones = sorted(capacity_df["zone"].unique())

techs = ["Hydro"]


# build csv with per zone capacity
flow_caps = [0.0] * len(zones)
for _, row in capacity_df.iterrows():
    if row["CP30 technology"].lower() == "hydro":
        zone = row["zone"]
        capacity = float(row["Installed Capacity (MWelec)"])
        if zone in zones:
            idx = zones.index(zone)
            flow_caps[idx] = capacity

cap_df = pd.DataFrame(
    [flow_caps, flow_caps],
    index=["flow_cap_min", "flow_cap_max"],
    columns=zones)

cap_df.to_csv(snakemake.output[1])

techs_yaml = {}

for tech in techs:
    tech_name = tech.lower()

    # flow_cap_data = tech_zone_capacity.get(tech_name, [0.0] * len(zones))
    # flow_cap_max = {"data": flow_cap_data, "dims": ["carriers"], "index": zones}
    techs_yaml[tech_name] = {
        "category": "renewable",
        "cp30_category": "renewable",
        "base_tech": "supply",
        "name": tech,
        "carrier_out": "electricity",
        "flow_out_eff": float(df.loc["efficiency", tech]),  # unitless (fraction)
        "lifetime": int(df.loc["lifetime", tech]),  # years
        "capacity_factor_max": 0.45,  # from Electricity Generation costs 2023 to avoid full load
        "cost_om_annual": {
            "data": float(df.loc["om_annual", tech]),
            "index": "monetary",
            "dims": ["costs"],
        },
        "cost_flow_out": {
            "data": float(df.loc["om_prod", tech]),
            "index": "monetary",
            "dims": ["costs"],
        },
        "cost_flow_in": {
            "data": float(df.loc["fuel_cost", tech]),
            "index": "monetary",
            "dims": ["costs"],
        },
        # Fix hydro capacity to today's values
        # "flow_cap_min": flow_cap_max,
        # "flow_cap_max": flow_cap_max,  
    }

data_tables_yaml = {
    "hydro_capacities": {
        "data": "data/processed/techs/hydro_capacities.csv",
        "rows": "parameters",
        "columns": "nodes",
        "add_dims": {
            "techs": tech_name
        }
    }
}

# Write output YAML
with open(snakemake.output[0], "w") as f:
    yaml.dump({"techs": techs_yaml, "data_tables": data_tables_yaml}, f, sort_keys=False)
