import pandas as pd
import yaml

df = pd.read_csv(snakemake.input[0], index_col=0)

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


for tech in fossil_techs:
    tech_name = tech.lower()

    techs_yaml[tech_name] = {
        "category": "fossil",         
        "cp30_category": "thermal", 
        
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

    # Existing tech: no capex, only operational/fuel costs
    techs_yaml[f"{tech_name}_existing"] = {
        "parent": tech_name,
        "costs": {
            "om_cost": float(df.loc["om_annual", tech]),
            "om_prod": float(df.loc["om_prod", tech]),
            "fuel": float(df.loc["fuel_cost", tech])
        }
    }

    # New tech: full cost structure
    if tech in allow_investment:
        techs_yaml[f"{tech_name}_new"] = {
            "parent": tech_name,
            "costs": {
                "cost_energy_cap": float(df.loc["capex", tech]),
                "om_cost": float(df.loc["om_annual", tech]),
                "om_prod": float(df.loc["om_prod", tech]),
                "fuel": float(df.loc["fuel_cost", tech])
        }
    }


with open(snakemake.output[0], "w") as f:
    yaml.dump(techs_yaml, f, sort_keys=False)