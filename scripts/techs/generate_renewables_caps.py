import pandas as pd
import yaml

df = pd.read_csv(snakemake.input[0])

df = df[df["CP30 technology"] != "Other_Renewables"]

tech_capacities = {"nodes": {}}

# To assign storage caps, with charging time in hrs
storage_techs = {
    "battery_existing": 4.0,
    "pumped_hydro_existing": 8,
    "caes_existing": 6,
    "laes_existing": 14,
}

for _, row in df.iterrows():
    tech = row["CP30 technology"].lower() + "_existing"
    zone = row["zone"]
    cap = float(row["Installed Capacity (MWelec)"])

    if zone not in tech_capacities["nodes"]:
        tech_capacities["nodes"][zone] = {"techs": {}}

    tech_dict = {"flow_cap_max": cap}

    if tech in storage_techs:
        multiplier = storage_techs[tech]
        storage_cap = cap * multiplier
        tech_dict["storage_cap_max"] = storage_cap

    tech_capacities["nodes"][zone]["techs"][tech] = tech_dict


with open(snakemake.output[0], "w") as f:
    yaml.dump(tech_capacities, f, sort_keys=False)
