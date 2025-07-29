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

        "essentials": {
            "name": tech,
            "carrier_in": fuels[tech],
            "carrier_out": "electricity"
        },
        "constraints": {
            "energy_eff": float(df.loc["efficiency", tech]),   # unitless (fraction)
            "lifetime": int(df.loc["lifetime", tech])          # years
        },
    }
    techs_yaml[f"{tech_name}_new"] = {
            "parent": tech_name,
            "costs": {
                "cost_energy_cap": float(df.loc["capex", tech]),
                "om_cost": float(df.loc["om_annual", tech]),
                "om_prod": float(df.loc["om_prod", tech]),
                "fuel": float(df.loc["fuel_cost", tech])
        }
    }

# Write output YAML
with open(snakemake.output[0], "w") as f:
    yaml.dump(techs_yaml, f, sort_keys=False)