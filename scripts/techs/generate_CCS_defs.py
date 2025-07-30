import pandas as pd
import yaml

df = pd.read_csv(snakemake.input[0], index_col=0)

techs = ["BECCS", "Gas_CCS"]

fuels = {
    "BECCS": "biomass",
    "Gas_CCS": "gas"
}

techs_yaml = {}

for tech in techs:
    tech_name = tech.lower()

    techs_yaml[tech_name] = {
        "category": "ccs",
        "cp30_category": "low_carbon",
        "base_tech": "conversion",
        "name": tech,
        "carrier_in": fuels[tech],
        "carrier_out": "electricity",
        "energy_eff": float(df.loc["efficiency", tech]),   # unitless (fraction)
        "lifetime": int(df.loc["lifetime", tech])          # years
    }

    techs_yaml[f"{tech_name}_new"] = {
        "parent": tech_name,
        "base_tech": "conversion",
        "cost_flow_cap": {  # changed from cost_energy_cap
            "data": float(df.loc["capex", tech]),
            "index": "monetary",
            "dims": ["costs"]
        },
        "cost_om_annual": {
            "data": float(df.loc["om_annual", tech]),
            "index": "monetary",
            "dims": ["costs"]
        },
        "cost_flow_out": {  # changed from cost_om_prod
            "data": float(df.loc["om_prod", tech]),
            "index": "monetary",
            "dims": ["costs"]
        },
        "cost_source": {
            "data": float(df.loc["fuel_cost", tech]),
            "index": "monetary",
            "dims": ["costs"]
        }
    }

# Write output YAML
with open(snakemake.output[0], "w") as f:
    yaml.dump({"techs": techs_yaml}, f, sort_keys=False)
