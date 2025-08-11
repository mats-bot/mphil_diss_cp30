import pandas as pd
import yaml

df = pd.read_csv(snakemake.input[0], index_col=0)
capacity_df = pd.read_csv(snakemake.input[1])

# remvoed coal since not present in 2030
fossil_techs = ["Gas_CCGT", "Gas_CCGT_CHP", "Gas_OCGT", "Diesel"]
techs_yaml = {}
templates_yaml = {}

# removal of fuels since defining as supply and not conversion
fuels = {
    "Gas_CCGT": "gas",
    "Gas_CCGT_CHP": "gas",
    "Gas_OCGT": "gas",
    "Diesel": "diesel",
    "Coal": "coal",
}

# Choose which techs can invest in
allow_investment = {"Gas_CCGT", "Gas_CCGT_CHP", "Gas_OCGT"}

zones = sorted(capacity_df["zone"].unique())
tech_zone_capacity = {}
for tech in capacity_df["CP30 technology"].unique():
    tech_zone_capacity[tech.lower()] = [0.0] * len(zones)

for _, row in capacity_df.iterrows():
    zone = row["zone"]
    tech = row["CP30 technology"].lower()
    capacity = float(row["InstalledCapacity (MW)"])
    if tech in tech_zone_capacity and zone in zones:
        idx = zones.index(zone)
        tech_zone_capacity[tech][idx] = capacity


for tech in fossil_techs:
    tech_name = tech.lower()

    templates_yaml[tech_name] = {
        "category": "fossil",
        "cp30_category": "thermal",
        "base_tech": "supply",
#        "carrier_in": fuels[tech],
        "carrier_out": "electricity",
        "flow_out_eff": float(df.loc["efficiency", tech]),
        "lifetime": int(df.loc["lifetime", tech]),
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
        "cost_flow_in": {  # fuel cost
            "data": float(df.loc["fuel_cost", tech]),
            "index": "monetary",
            "dims": ["costs"],
        },
    }

    # flow_cap_max = {
    #     "data": tech_zone_capacity.get(tech_name, [0.0] * len(zones)),
    #     "dims": ["carriers"],
    #     "index": zones,
    # }

    # Existing tech: no capex, only operational/fuel costs
    techs_yaml[f"{tech_name}_existing"] = {
        "template": tech_name,
        # "flow_cap_max": flow_cap_max,
    }

    # New tech: full cost structure
    if tech in allow_investment:
        techs_yaml[f"{tech_name}_new"] = {
            "template": tech_name,
            "cost_flow_cap": {
                "data": float(df.loc["capex", tech]),
                "index": "monetary",
                "dims": ["costs"],
            },
        }

with open(snakemake.output[0], "w") as f:
    yaml.dump({"techs": techs_yaml, "templates": templates_yaml}, f, sort_keys=False)
