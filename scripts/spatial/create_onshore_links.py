import yaml
import pandas as pd

df = pd.read_csv(snakemake.input[0])

# Set default year for GTCs
with open("config/default.yaml") as f:
    config = yaml.safe_load(f)
target_year = str(config["default_year"])

links = {}
for _, row in df.iterrows():
    zone_pair = "-".join(sorted([row["Zone 1"], row["Zone 2"]]))
    
    constraints = {
        "techs": ["ac_transmission"],
        "constraints": {}
    }
    
    if not pd.isna(row.loc[target_year]):
        constraints["constraints"]["energy_cap_max"] = float(row.loc[target_year]) / 1000  # convert to GW
    
    if not pd.isna(row["loss"]):
        constraints["constraints"]["loss"] = float(row["loss"])  
    
    links[zone_pair] = constraints


with open(snakemake.output[0], 'w') as f:
    yaml.dump({"links": links}, f, sort_keys=False, indent=2)