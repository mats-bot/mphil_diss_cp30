import pandas as pd
import yaml

df = pd.read_csv(snakemake.input[0])

df = df[df["CP30 technology"] != "Other_Renewables"]

tech_capacities = {}

for _, row in df.iterrows():
    tech = row["CP30 technology"]
    zone = row["zone"]
    cap = float(row["Installed Capacity (MWelec)"])

    if tech not in tech_capacities:
        tech_capacities[tech] = {"constraints": {"energy_cap": {}}}
    tech_capacities[tech]["constraints"]["energy_cap"][zone] = cap

with open(snakemake.output[0], "w") as f:
    yaml.dump(tech_capacities, f, sort_keys=False)