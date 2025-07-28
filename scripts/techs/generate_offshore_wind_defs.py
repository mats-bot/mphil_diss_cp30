import pandas as pd
import yaml

costs_df  = pd.read_csv(snakemake.input[0], index_col=0)
projects_df = pd.read_csv(snakemake.input[1])

cfs_path = str(snakemake.input[2])

projects_df["Operational"] = pd.to_datetime(projects_df["Operational"], errors="coerce")
cutoff_date = pd.Timestamp("2024-01-01")

capex = float(costs_df.at["capex", "Offshore_Wind"])
om_annual = float(costs_df.at["om_annual", "Offshore_Wind"])
om_prod = float(costs_df.at["om_prod", "Offshore_Wind"])
lifetime = int(float(costs_df.at["lifetime", "Offshore_Wind"]))


techs = {}

for _, row in projects_df.iterrows():
    site = row["Site Name"]
    slug = site.replace(" ", "_").replace("-", "_")
    tech_name = f"OffshoreWind_{slug}"
    installed_cap = row["Installed Capacity (MW)"]

    status = row["Development Status (short)"]
    op_date = row["Operational"]

    constraints = {
        "resource": cfs_path,
        "resource_column": site,
        "resource_unit": "per_unit",
        "resource_area": site,
        "lifetime": lifetime,
        "energy_cap_per_unit": installed_cap,
        "integer": True,
        "units_max": 1
    }

    # CATEGORY 1: Operational by end of 2023
    if pd.notnull(op_date) and op_date < cutoff_date:
        constraints["energy_cap_equals"] = installed_cap
        constraints.pop("energy_cap_per_unit", None)
        constraints.pop("integer", None)
        constraints.pop("units_max", None)

    # CATEGORY 2: Operational between 2024-2025, in/awaiting construction
    elif status in {"Operational", "Under Construction", "Awaiting Construction"}:
        constraints["units_min"] = 1

    # CATEGORY 3: Projects in planning phase
    else:
        pass

    techs[tech_name] = {
        "category": "renewable",
        "cp30_category": "renewable",
        "essentials": {
            "name": tech_name,
            "carrier_out": "electricity"
        },
        "constraints": constraints,
        "costs": {
            "cost_energy_cap": capex,
            "om_annual": om_annual,
            "om_prod": om_prod
        }
    }

with open(snakemake.output[0], "w") as f:
    yaml.dump(techs, f, sort_keys=False)