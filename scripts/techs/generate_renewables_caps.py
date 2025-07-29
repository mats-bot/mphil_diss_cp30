import pandas as pd
import yaml

df = pd.read_csv(snakemake.input[0])

df = df[df["CP30 technology"] != "Other_Renewables"]

# Choose which techs cannot invest further
barred = {"Hydro", "Pumped_Hydro"}

tech_capacities = {}

for _, row in df.iterrows():
    tech = row["CP30 technology"].lower()
    zone = row["zone"]
    cap = float(row["Installed Capacity (MWelec)"])

    if tech not in tech_capacities:
        tech_capacities[tech] = {
            "constraints": {
                "energy_cap": {},
                "energy_cap_min": {},
                "energy_cap_max": {}
            }
        }

    if tech in barred:
        # Fix maximum cap (no new build)
        tech_capacities[tech]["constraints"]["energy_cap_max"][zone] = cap
    else:
        # Set minimum capacity, allow new investment (no max limit)
        tech_capacities[tech]["constraints"]["energy_cap_min"][zone] = cap

with open(snakemake.output[0], "w") as f:
    yaml.dump(tech_capacities, f, sort_keys=False)