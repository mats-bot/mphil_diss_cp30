import os
import yaml

# From FES24 regional breakdown
type_labels = {
    "C": "Commercial",
    "I": "Industrial",
    "R": "Residential",
}

def generate_demand_yaml(input_dir, output_path):
    techs = {}
    data_tables = {}
    nodes = {}

    # Parent tech
    techs["demand"] = {
        "name": "Overall demand technology parent",
        "carrier_in": "electricity"
    }

    for filename in os.listdir(input_dir):
        if not filename.endswith(".csv"):
            continue

        parts = filename.replace(".csv", "").split("_")
        if len(parts) != 3:
            continue

        _, zone, dtype = parts
        dtype_label = type_labels.get(dtype, dtype)
        tech_name = f"demand_{dtype}"

        # Add tech if not already
        if tech_name not in techs:
            techs[tech_name] = {
                "name": f"{dtype_label} electricity demand",
                "base_tech": "demand",
                "carrier_in": "electricity"
            }

        filepath = os.path.relpath(os.path.join(input_dir, filename), os.getcwd())

        dt_key = f"{tech_name}_{zone}"
        data_tables[dt_key] = {
            "data": filepath.replace("\\", "/"),  # Use forward slashes in YAML
            "rows": "timesteps",
            "columns": "nodes",
            "add_dims": {
                "techs": tech_name,
                "parameters": "sink_use_equals"
            }
        }

        # Add node entry for zone and tech
        node_key = f"{zone}.techs.{tech_name}"
        nodes[node_key] = None

    # Compose full dict to dump
    out_dict = {
        "techs": techs,
        "data_tables": data_tables,
        "nodes": nodes
    }

    with open(output_path, "w") as f:
        yaml.dump(out_dict, f, sort_keys=False, default_flow_style=False, indent=4)

generate_demand_yaml(snakemake.input[0], snakemake.output[0])