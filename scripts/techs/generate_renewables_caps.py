import pandas as pd
import yaml

df = pd.read_csv(snakemake.input[0])

df = df[df["CP30 technology"] != "Other_Renewables"]

tech_capacities = {"locations": {}}

for _, row in df.iterrows():
    tech = row["CP30 technology"].lower() + "_existing"
    zone = row["zone"]
    cap = float(row["Installed Capacity (MWelec)"])

    if zone not in tech_capacities["locations"]:
        tech_capacities["locations"][zone] = {"techs": {}}
    
    tech_capacities["locations"][zone]["techs"][tech] = {
        "constraints": {
            "energy_cap_max": cap
        }
    }

with open(snakemake.output[0], "w") as f:
    yaml.dump(tech_capacities, f, sort_keys=False)