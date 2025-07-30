import pandas as pd
import yaml

df = pd.read_csv(snakemake.input[0], index_col=0)
capacity_df = pd.read_csv(snakemake.input[1])
zones = sorted(capacity_df['zone'].unique())

techs = ["Hydro"]

tech_zone_capacity = {}
for tech in capacity_df['CP30 technology'].unique():
    tech_zone_capacity[tech.lower()] = [0.0] * len(zones)

for _, row in capacity_df.iterrows():
    zone = row['zone']
    tech = row['CP30 technology'].lower()
    capacity = float(row['Installed Capacity (MWelec)'])
    if tech in tech_zone_capacity and zone in zones:
        idx = zones.index(zone)
        tech_zone_capacity[tech][idx] = capacity

techs_yaml = {}

for tech in techs:
    tech_name = tech.lower()

    techs_yaml[tech_name] = {
        "category": "renewable",
        "cp30_category": "renewable",
        "base_tech": "supply",
        "name": tech,
        "carrier_out": "electricity",
        "energy_eff": float(df.loc["efficiency", tech]),   # unitless (fraction)
        "lifetime": int(df.loc["lifetime", tech]),         # years
        "capacity_factor_max": 0.45  # from Electricity Generation costs 2023 to avoid full load
    }

    flow_cap_data = tech_zone_capacity.get(tech_name, [0.0] * len(zones))
    flow_cap_max = {
        "data": flow_cap_data,
        "dims": ["carriers"],
        "index": zones,
    }

    techs_yaml[f"{tech_name}_existing"] = {
        "parent": tech_name,
        "base_tech": "supply",
        "cost_om_annual": {
            "data": float(df.loc["om_annual", tech]),
            "index": "monetary",
            "dims": ["costs"]
        },
        "cost_flow_out": {
            "data": float(df.loc["om_prod", tech]),
            "index": "monetary",
            "dims": ["costs"]
        },
        "cost_source": {
            "data": float(df.loc["fuel_cost", tech]),
            "index": "monetary",
            "dims": ["costs"]
        },
        "flow_cap_max": flow_cap_max,
    }


# Write output YAML
with open(snakemake.output[0], "w") as f:
    yaml.dump({"techs": techs_yaml}, f, sort_keys=False)
