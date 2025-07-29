import yaml

# Load techs.yaml, spatial.yaml, demand.yaml
with open(snakemake.input.techs) as f:
    techs = yaml.safe_load(f)
    # If no top-level key "techs", wrap it
    if "techs" not in techs:
        techs = {"techs": techs}

with open(snakemake.input.spatial) as f:
    spatial = yaml.safe_load(f)
    # Wrap if no top-level key (assumed 'spatial' here)
    if "spatial" not in spatial:
        spatial = {"spatial": spatial}

with open(snakemake.input.demand) as f:
    demand = yaml.safe_load(f)
    # demand already has 'techs' top-level key, so probably it's ok

# Build full model dictionary
model = {
    "model": {
        "name": "CP30-main-rep",
        "years": [2030],
        "timestep": "1h",
        "calliope_version": "0.7.0.dev6"
    },
    "techs": techs["techs"],
    "spatial": spatial["spatial"],
    "demand": demand.get("techs", demand),  # use demand['techs'] or whole demand dict
    "constraints": {
        "limit_dirty_generation": {
            "description": "Limit dirty generation to 5% of total generation",
            "foreach": ["timesteps"],
            "expression": (
                "sum(model_output:flow, where={\"techs__category\": \"fossil\"}, dims=[\"techs\", \"nodes\"]) "
                "<= 0.05 * sum(model_output:flow, dims=[\"techs\", \"nodes\"])"
            )
        }
    },
    "LDES": {
        "essentials": {
            "name": "LDES",
            "storage": True
        },
        "children": [
            "Pumped_Hydro",
            "CAES",
            "LAES"
        ]
    }
}

# Write out combined model.yaml
with open(snakemake.output[0], "w") as f:
    yaml.dump(model, f, sort_keys=False)
