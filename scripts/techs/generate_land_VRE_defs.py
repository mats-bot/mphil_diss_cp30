import pandas as pd
import yaml

df = pd.read_csv(snakemake.input[0], index_col=0)

renewable_techs = ["Onshore_Wind", "Solar_PV"]

solar_cfs = pd.read_csv(snakemake.input[1])
onshore_cfs = pd.read_csv(snakemake.input[2])


year = int(snakemake.config["weather_year"])

resources = {
    "Solar_PV": f"data/processed/spatial/solar_cf_{year}.csv",
    "Onshore_Wind": f"data/processed/spatial/onshore_cf_{year}.csv"
}


techs_yaml = {}

for tech in renewable_techs:
    tech_base = tech.lower()

    techs_yaml[tech_base] = {
        "category": "renewable",
        "cp30_category": "renewable",
        "essentials": {
            "name": tech,
            "carrier_out": "electricity"
        },
        "constraints": {
            "energy_eff": float(df.loc["efficiency", tech]),
            "resource": resources[tech],
            "resource_unit": "per_unit",
            "lifetime": int(df.loc["lifetime", tech])
        }
    }

    # Existing sites 
    techs_yaml[f"{tech_base}_existing"] = {
        "parent": tech_base,
        "costs": {
            "om_annual": float(df.loc["om_annual", tech]),
            "om_prod": float(df.loc["om_prod", tech])
        }
    }

    # New projects
    techs_yaml[f"{tech_base}_new"] = {
            "parent": tech_base,
            "costs": {
                "cost_energy_cap": float(df.loc["capex", tech]),
                "om_annual": float(df.loc["om_annual", tech]),
                "om_prod": float(df.loc["om_prod", tech])
            }
        }

with open(snakemake.output[0], "w") as f:
    yaml.dump(techs_yaml, f, sort_keys=False)