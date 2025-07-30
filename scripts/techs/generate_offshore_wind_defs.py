import re

import pandas as pd
import yaml

costs_df = pd.read_csv(snakemake.input[0], index_col=0)
projects_df = pd.read_csv(snakemake.input[1])

cfs_path = str(snakemake.input[2])

projects_df["Operational"] = pd.to_datetime(projects_df["Operational"], errors="coerce")
cutoff_date = pd.Timestamp("2024-01-01")

capex = float(costs_df.at["capex", "Offshore_Wind"])
om_annual = float(costs_df.at["om_annual", "Offshore_Wind"])
om_prod = float(costs_df.at["om_prod", "Offshore_Wind"])
lifetime = int(float(costs_df.at["lifetime", "Offshore_Wind"]))

# Parent tech
template = {
    "offshore_wind": {
        "category": "renewable",
        "cp30_category": "renewable",
        "base_tech": "supply",
        "name": "offshore_wind",
        "carrier_out": "electricity",
        "source_unit": "per_cap",
        "lifetime": lifetime,
        # "resource": cfs_path, # TODO: use data table
        "cost_flow_cap": {"data": capex, "index": "monetary", "dims": ["costs"]},
        "cost_om_annual": {"data": om_annual, "index": "monetary", "dims": ["costs"]},
        "cost_flow_out": {"data": om_prod, "index": "monetary", "dims": ["costs"]},
    }
}
techs = {}


# Func to clean names for calliope
def sanitize_tech_name(name):
    name = re.sub(r"\([^)]*\)", "", name)
    name = re.sub(r"[ \-]", "_", name)
    name = re.sub(r"[^\w]", "", name)
    name = name.lower().strip("_")
    if not name or not name[0].isalpha():
        name = "t_" + name
    return name


cf_columns = pd.read_csv(cfs_path, nrows=1).columns

# Site specific
for _, row in projects_df.iterrows():
    site = row["Site Name"]
    slug = sanitize_tech_name(site)
    tech_name = f"offshorewind_{slug}"

    installed_cap = row["Installed Capacity (MWelec)"]

    status = row["Development Status (short)"]
    op_date = row["Operational"]

    # Match site names explicitly since funky characters
    matched_column = next(
        (col for col in cf_columns if col.strip() == site.strip()), None
    )
    if matched_column is None:
        raise ValueError(f"Could not match '{site}' to any column in {cfs_path}")

    tech = {
        "template": "offshore_wind",
        "resource_column": matched_column,
        "flow_cap_per_unit": installed_cap,
        "units_max": 1,
    }

    # Category 1 override: Operational by end of 2023
    if pd.notnull(op_date) and op_date < cutoff_date:
        tech.update(
            {
                "flow_cap_min": installed_cap,
                "flow_cap_max": installed_cap,
                "cost_flow_cap": 0,  # keep scalar zero here
            }
        )
        # remove keys not needed when energy_cap_equals is used
        tech.pop("flow_cap_per_unit", None)
        tech.pop("units_max", None)
        tech.pop("integer", None)

    # Category 2 override: Will be built by 2030
    elif status in {"Operational", "Under Construction", "Awaiting Construction"}:
        tech["units_min"] = 1

    # Category 3 (all others): use constraints as is (model decides)
    techs[tech_name] = tech

with open(snakemake.output[0], "w") as f:
    yaml.dump({"techs": techs, "templates": template}, f, sort_keys=False)
