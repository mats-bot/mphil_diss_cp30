import pandas as pd 
import pandas as pd
import yaml

df = pd.read_csv(snakemake.input[0], index_col=0)

techs = ["Hydro"]

techs_yaml = {}

for tech in techs:
    techs_yaml[tech.lower()] = {
        "category": "renewable",
        "cp30_category": "renewable",  

        "essentials": {
            "name": tech,
            "carrier_out": "electricity"
        },
        "constraints": {
            "energy_eff": float(df.loc["efficiency", tech]),   # unitless (fraction)
            "lifetime": int(df.loc["lifetime", tech]),         # years
            "capacity_factor_max": 0.45 # from Electricity Generation costs 2023 to avoid full load
        },
        "costs": {
            "cost_energy_cap": float(df.loc["capex", tech]),        # £/kW installed
            "om_cost": float(df.loc["om_annual", tech]),     # £/kW/year
            "om_prod": float(df.loc["om_prod", tech]),         # £/kWh
        }
    }

# Write output YAML
with open(snakemake.output[0], "w") as f:
    yaml.dump(techs_yaml, f, sort_keys=False)