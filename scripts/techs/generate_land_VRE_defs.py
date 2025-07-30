import pandas as pd
import yaml

df = pd.read_csv(snakemake.input[0], index_col=0)

renewable_techs = ["Onshore_Wind", "Solar_PV"]

solar_cfs = pd.read_csv(snakemake.input[1])
onshore_cfs = pd.read_csv(snakemake.input[2])

capacity_df = pd.read_csv(snakemake.input[3])
zones = sorted(capacity_df['zone'].unique())

year = int(snakemake.config["weather_year"])

resources = {
    "Solar_PV": f"data/processed/spatial/solar_cf_{year}.csv",
    "Onshore_Wind": f"data/processed/spatial/onshore_cf_{year}.csv"
}

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

for tech in renewable_techs:
    tech_base = tech.lower()

    techs_yaml[tech_base] = {
        "category": "renewable",
        "cp30_category": "renewable",
        "name": tech,
        "base_tech": "supply",
        "carrier_out": "electricity",
        "energy_eff": float(df.loc["efficiency", tech]),
        "resource": resources[tech],
        "resource_unit": "per_unit",
        "lifetime": int(df.loc["lifetime", tech])
    }

    flow_cap_max = {
        "data": tech_zone_capacity.get(tech_base, [0.0] * len(zones)),
        "dims": ["carriers"],
        "index": zones
    }

    techs_yaml[f"{tech_base}_existing"] = {
        "parent": tech_base,
        "om_annual": float(df.loc["om_annual", tech]),
        "om_prod": float(df.loc["om_prod", tech])
    }

    techs_yaml[f"{tech_base}_new"] = {
        "parent": tech_base,
        "cost_energy_cap": float(df.loc["capex", tech]),
        "om_annual": float(df.loc["om_annual", tech]),
        "om_prod": float(df.loc["om_prod", tech])
    }

with open(snakemake.output[0], "w") as f:
    yaml.dump({"techs": techs_yaml}, f, sort_keys=False)