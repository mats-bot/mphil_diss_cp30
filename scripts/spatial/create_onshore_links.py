import yaml
import pandas as pd

df = pd.read_csv(snakemake.input[0])

# Set default year for GTCs
with open("config/default.yaml") as f:
    config = yaml.safe_load(f)
target_year = str(config["default_year"])

techs = {}
for _, row in df.iterrows():
    zone_1 = row["Zone 1"].strip()
    zone_2 = row["Zone 2"].strip()

    tech_name = "_".join(sorted([zone_1, zone_2]))
    
    tech = {
        "parent": "ac_transmission",
        "link_from": zone_1,
        "link_to": zone_2,
    }

    if not pd.isna(row.get(target_year)):
        tech["flow_cap_max"] = float(row[target_year])

    if not pd.isna(row.get("loss")):
        tech["loss"] = float(row["loss"])

    techs[tech_name] = tech

with open(snakemake.output[0], 'w') as f:
    yaml.dump({"techs": techs}, f, sort_keys=False, indent=2)