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
    techs_yaml[tech.lower()] = {
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
        "costs": {
            "energy_cap": float(df.loc["capex", tech]),        # £/kW installed
            "om_annual": float(df.loc["om_annual", tech]),     # £/kW/year
            "om_prod": float(df.loc["om_prod", tech]),         # £/kWh
            "fuel": float(df.loc["fuel_cost", tech])           # £/kWh fuel
        }
    }

# Write output YAML
with open(snakemake.output[0], "w") as f:
    yaml.dump(techs_yaml, f, sort_keys=False)