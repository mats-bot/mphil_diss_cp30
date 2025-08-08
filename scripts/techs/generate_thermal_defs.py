import pandas as pd
import yaml

df = pd.read_csv(snakemake.input[0], index_col=0)

techs = ["Biomass", "Biomass_CHP", "Waste", "Waste_CHP", "Hydrogen", "Nuclear"]

fuels = {
    "Biomass": "biomass",
    "Biomass_CHP": "biomass",
    "Waste": "waste",
    "Waste_CHP": "waste",
    "Hydrogen": "hydrogen",
    "Nuclear": "nuclear_fuel",
}

# Per cp30 outputs
classif = {
    "Biomass": "renewable",
    "Biomass_CHP": "renewable",
    "Waste": "renewable",
    "Waste_CHP": "renewable",
    "Hydrogen": "thermal",
    "Nuclear": "low_carbon",
}

disallow_investment = {"nuclear"}

capacity_df = pd.read_csv(snakemake.input[1])  # renewable capacities CSV
nuclear_capacity_df = pd.read_csv(snakemake.input[2])  # nuclear/fossil capacities CSV
zones = sorted(capacity_df["zone"].unique())


def build_capacity_dict(df, cap_col_name):
    tech_zone_capacity = {}
    for tech in df["CP30 technology"].unique():
        tech_zone_capacity[tech.lower()] = [0.0] * len(zones)
    for _, row in df.iterrows():
        zone = row["zone"]
        tech = row["CP30 technology"].lower()
        capacity = float(row[cap_col_name])
        if tech in tech_zone_capacity and zone in zones:
            idx = zones.index(zone)
            tech_zone_capacity[tech][idx] = capacity
    return tech_zone_capacity


main_capacity = build_capacity_dict(capacity_df, "Installed Capacity (MWelec)")
nuclear_capacity = build_capacity_dict(nuclear_capacity_df, "InstalledCapacity (MW)")

techs_yaml = {}

for tech in techs:
    tech_name = tech.lower()

    if tech_name == "nuclear":
        flow_cap_data = nuclear_capacity.get(tech_name, [0.0] * len(zones))
    else:
        flow_cap_data = main_capacity.get(tech_name, [0.0] * len(zones))

    flow_cap_max = {"data": flow_cap_data, "dims": ["nodes"], "index": zones}

    # Set CAPEX for nuclear to 0 since all new build, and set min flow to 0.5 of cap
    if tech_name == "nuclear":
        flow_cap_data = nuclear_capacity.get(tech_name, [0.0] * len(zones))
        cost_flow_cap_val = 0.0
        flow_out_min_relative = 0.5
    else:
        flow_cap_data = main_capacity.get(tech_name, [0.0] * len(zones))
        cost_flow_cap_val = float(df.loc["capex", tech])
        flow_out_min_relative = ".null"

    techs_yaml[tech_name] = {
        "category": "thermal",
        "cp30_category": classif[tech],
        "base_tech": "conversion",
        "name": tech,
        "carrier_in": fuels[tech],
        "carrier_out": "electricity",
        "flow_out_eff": float(df.loc["efficiency", tech]),
        "lifetime": int(df.loc["lifetime", tech]),
        "flow_out_min_relative": flow_out_min_relative,
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
        "cost_flow_cap": {
            "data": cost_flow_cap_val,
            "index": "monetary",
            "dims": ["costs"],
        },
    }

with open(snakemake.output[0], "w") as f:
    yaml.dump({"techs": techs_yaml}, f, sort_keys=False)
