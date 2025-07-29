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

# Parent tech
techs = {
    "offshore_wind": {
        "category": "renewable",
        "cp30_category": "renewable",
        "essentials": {
            "name": "offshore_wind",
            "carrier_out": "electricity"
        },
        "constraints": {
            "resource_unit": "per_unit",
            "lifetime": lifetime
        },
        "costs": {
            "cost_energy_cap": capex,
            "om_annual": om_annual,
            "om_prod": om_prod
        }
    }
}

cf_columns = pd.read_csv(cfs_path, nrows=1).columns
# Site specific
for _, row in projects_df.iterrows():
    site = row["Site Name"]
    slug = site.replace(" ", "_").replace("-", "_")
    tech_name = f"OffshoreWind_{slug}"
    tech_name = tech_name.lower()
    installed_cap = row["Installed Capacity (MWelec)"]

    status = row["Development Status (short)"]
    op_date = row["Operational"]

    # Match site names explicitly since funky characters
    matched_column = next((col for col in cf_columns if col.strip() == site.strip()), None)
    if matched_column is None:
        raise ValueError(f"Could not match '{site}' to any column in {cfs_path}")

    # Start with a copy of parent constraints and update with site-specific
    costs = techs["offshore_wind"]["costs"].copy()
    constraints = techs["offshore_wind"]["constraints"].copy()
    
    constraints.update({
        "resource": cfs_path,
        "resource_column": matched_column,
        "resource_area": matched_column,
        "energy_cap_per_unit": installed_cap,
        "integer": True,
        "units_max": 1
    })

    # Category 1 override: Operational by end of 2023
    if pd.notnull(op_date) and op_date < cutoff_date:
        constraints = {
            "resource": cfs_path,
            "resource_column": matched_column,
            "resource_area": matched_column,
            "energy_cap_equals": installed_cap
        }
        costs["cost_energy_cap"] = 0  # Override capex for built projects

    # Category 2 override: Will be built by 2030
    elif status in {"Operational", "Under Construction", "Awaiting Construction"}:
        constraints["units_min"] = 1

    # Category 3 (all others): use constraints as is (model decides)
    tech = {
        "parent_category": "offshore_wind",
        "essentials": {
            "name": tech_name,
            "location": row["tzone"]
        },
        "constraints": constraints,
        "costs": costs
    }


    techs[tech_name] = tech

with open(snakemake.output[0], "w") as f:
    yaml.dump(techs, f, sort_keys=False)