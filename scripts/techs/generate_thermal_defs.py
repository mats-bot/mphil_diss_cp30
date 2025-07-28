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


new_build_cap = {
    "Nuclear": 0
}

# Initialize YAML dict
techs_yaml = {}

for tech in techs:
    techs_yaml[tech.lower()] = {
        "category": "thermal",
        "cp30_category": classif[tech],  

        "essentials": {
            "name": tech,
            "carrier_in": fuels[tech],
            "carrier_out": "electricity"
        },
        "constraints": {
            "energy_cap_max": new_build_cap.get(tech, float("nan")),
            "energy_eff": float(df.loc["efficiency", tech]),   # unitless (fraction)
            "lifetime": int(df.loc["lifetime", tech])          # years
        },
        "costs": {
            "cost_energy_cap": float(df.loc["capex", tech]),        # £/kW installed
            "om_cost": float(df.loc["om_annual", tech]),     # £/kW/year
            "om_prod": float(df.loc["om_prod", tech]),         # £/kWh
            "fuel": float(df.loc["fuel_cost", tech])           # £/kWh fuel
        }
    }

# Write output YAML
with open(snakemake.output[0], "w") as f:
    yaml.dump(techs_yaml, f, sort_keys=False)