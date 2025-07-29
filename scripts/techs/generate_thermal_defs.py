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
    "Nuclear": "nuclear_fuel"
}

# Per cp30 outputs
classif = {
    "Biomass": "renewable",
    "Biomass_CHP": "renewable",
    "Waste": "renewable",
    "Waste_CHP": "renewable",
    "Hydrogen": "thermal",
    "Nuclear": "low_carbon"
}


disallow_investment = {"nuclear"}


# Initialize YAML dict
techs_yaml = {}

for tech in techs:
    tech_name = tech.lower()

    techs_yaml[tech_name] = {
        "category": "thermal",
        "cp30_category": classif[tech],  

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
    if tech not in disallow_investment:
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