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


for tech in fossil_techs:
    techs_yaml[tech.lower()] = {
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
        "costs": {
            "cost_energy_cap": float(df.loc["capex", tech]),        # £/kW installed
            "om_cost": float(df.loc["om_annual", tech]),     # £/kW/year
            "om_prod": float(df.loc["om_prod", tech]),         # £/kWh
            "fuel": float(df.loc["fuel_cost", tech])           # £/kWh fuel
        }
    }



with open(snakemake.output[0], "w") as f:
    yaml.dump(techs_yaml, f, sort_keys=False)