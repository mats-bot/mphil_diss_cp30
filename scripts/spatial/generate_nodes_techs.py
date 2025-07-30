import re

import pandas as pd
import yaml


# Func to clean names for calliope
def sanitize_tech_name(name):
    name = re.sub(r"\([^)]*\)", "", name)
    name = re.sub(r"[ \-]", "_", name)
    name = re.sub(r"[^\w]", "", name)
    name = name.lower().strip("_")
    if not name or not name[0].isalpha():
        name = "t_" + name
    return name


with open(snakemake.input.nodes_file) as f:
    locations = yaml.safe_load(f).get("nodes", {})

nodes = list(locations.keys())


all_techs = {}
for file_path in snakemake.input.tech_files:
    with open(file_path) as f:
        tech_data = yaml.safe_load(f)
        if tech_data and "techs" in tech_data:
            all_techs.update(tech_data["techs"])

offshore_parent = "offshorewind"
offshore_children = {k for k in all_techs.keys() if k.startswith(offshore_parent + "_")}
other_techs = set(all_techs.keys()) - offshore_children - {offshore_parent}

df_sites = pd.read_csv(snakemake.input.offshore_df)
df_sites["tzone"] = df_sites["tzone"].str.strip()

site_to_zone = {}
for _, row in df_sites.iterrows():
    slug = sanitize_tech_name(row["Site Name"])
    tech_name = f"{offshore_parent}_{slug}"
    site_to_zone[tech_name] = row["tzone"]

techs_per_node = {node: {} for node in nodes}

for tech in other_techs:
    if tech.endswith("_existing") or tech.endswith("_new"):
        continue
    for node in nodes:
        techs_per_node[node][tech] = {}

for tech in offshore_children:
    node = site_to_zone.get(tech)
    if node in nodes:
        techs_per_node[node][tech] = {}


nodes_yaml = {}
for node in nodes:
    nodes_yaml[node] = {
        "coordinates": locations[node]["coordinates"],
        "techs": techs_per_node.get(node, {}),
    }

with open(snakemake.output[0], "w") as f:
    yaml.dump({"nodes": nodes_yaml}, f, sort_keys=True)
