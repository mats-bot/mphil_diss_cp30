import pandas as pd
import yaml

df = pd.read_csv(snakemake.input[0], index_col=0)
capacity_df = pd.read_csv(snakemake.input[1]) 

fossil_techs = ["Gas_CCGT", "Gas_CCGT_CHP", "Gas_OCGT", "Diesel", "Coal"]
techs_yaml = {}

fuels = {
    "Gas_CCGT": "gas",
    "Gas_CCGT_CHP": "gas",
    "Gas_OCGT": "gas",
    "Diesel": "diesel",
    "Coal": "coal"
}

# Choose which techs can invest in
allow_investment = {"Gas_CCGT", "Gas_CCGT_CHP", "Gas_OCGT"}

zones = sorted(capacity_df['zone'].unique())
tech_zone_capacity = {}
for tech in capacity_df['CP30 technology'].unique():
    tech_zone_capacity[tech.lower()] = [0.0] * len(zones)

for _, row in capacity_df.iterrows():
    zone = row['zone']
    tech = row['CP30 technology'].lower()
    capacity = float(row['InstalledCapacity (MW)'])
    if tech in tech_zone_capacity and zone in zones:
        idx = zones.index(zone)
        tech_zone_capacity[tech][idx] = capacity


for tech in fossil_techs:
    tech_name = tech.lower()

    techs_yaml[tech_name] = {
        "category": "fossil",
        "cp30_category": "thermal",
        "base_tech": "conversion",
        "name": tech,
        "carrier_in": fuels[tech],
        "carrier_out": "electricity",
        "energy_eff": float(df.loc["efficiency", tech]),
        "lifetime": int(df.loc["lifetime", tech]),
    }

    flow_cap_max = {
        "data": tech_zone_capacity.get(tech_name, [0.0] * len(zones)),
        "dims": ["carriers"],
        "index": zones
    }

    # Existing tech: no capex, only operational/fuel costs
    techs_yaml[f"{tech_name}_existing"] = {
        "parent": tech_name,
        "base_tech": "conversion",
        "om_annual": float(df.loc["om_annual", tech]),
        "om_prod": float(df.loc["om_prod", tech]),
        "fuel": float(df.loc["fuel_cost", tech]),
        "flow_cap_max": flow_cap_max
    }

    # New tech: full cost structure
    if tech in allow_investment:
        techs_yaml[f"{tech_name}_new"] = {
            "parent": tech_name,
            "base_tech": "conversion",
            "cost_energy_cap": float(df.loc["capex", tech]),
            "om_annual": float(df.loc["om_annual", tech]),
            "om_prod": float(df.loc["om_prod", tech]),
            "fuel": float(df.loc["fuel_cost", tech]),
        }


with open(snakemake.output[0], "w") as f:
    yaml.dump({"techs": techs_yaml}, f, sort_keys=False)