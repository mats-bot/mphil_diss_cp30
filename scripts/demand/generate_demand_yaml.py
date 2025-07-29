import os
import yaml

# From FES24 regional breakdown
type_labels = {
    "R": "Residential",
    "C": "Commercial",
    "I": "Industrial",
    "T": "Transmission direct connects",
    "D": "District heating",
    "E": "Electric vehicles",
    "H": "Heat pumps",
    "Z": "Electrolyzers"
}

def generate_demand_yaml(input_dir, output_path):
    demand_data = {}

    for filename in os.listdir(input_dir):
        parts = filename.replace(".csv", "").split("_")
        if len(parts) != 3:
            continue

        _, zone, dtype = parts
        label = type_labels.get(dtype, dtype)
        tech_name = f"demand_{dtype}"

        if tech_name not in demand_data:
            demand_data[tech_name] = {
                "description": f"{label} electricity demand",
                "essentials": {
                    "carrier": "electricity",
                    "parent": "demand"
                },
                "constraints": {}
            }

        demand_data[tech_name]["constraints"][f"sink_use_equals::{zone}"] = f"file={filename}"

    yaml_output = {"techs": demand_data}

    with open(output_path, "w") as f:
        yaml.dump(yaml_output, f, sort_keys=False)


generate_demand_yaml(snakemake.input[0], snakemake.output[0])