import pandas as pd
import yaml

df = pd.read_csv(snakemake.input[0])

df["zone"] = df["zone"].str.lower()

locations_yaml = {"locations": {}}

for _, row in df.iterrows():
    zone = row["zone"]
    tech = row["CP30 technology"].lower() + "_existing"
    cap = float(row["InstalledCapacity (MW)"])

    if zone not in locations_yaml["locations"]:
        locations_yaml["locations"][zone] = {"techs": {}}
    
    locations_yaml["locations"][zone]["techs"][tech] = {
        "constraints": {
            "energy_cap_max": cap
        }
    }

with open(snakemake.output[0], "w") as f:
    yaml.dump(locations_yaml, f, sort_keys=False, default_flow_style=False)