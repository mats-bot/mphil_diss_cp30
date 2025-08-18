import os
import yaml
import re

type_labels = {
    "C": "Commercial",
    "I": "Industrial",
    "R": "Residential",
    "E": "Electric Vehicles",
    "H": "Heat Pumps",
    "D": "District heat",
    "T": "Transmission direct connects",
    "Z": "Electrolyzers"
}

def sanitize_name(name):
    return re.sub(r'[^a-zA-Z0-9]', '_', name)

def generate_demand_yaml(input_dir, output_path):
    techs = {}
    data_tables = {}

    # Base parent tech
    techs["demand"] = {
        "name": "Overall demand technology parent",
        "carrier_in": "electricity"
    }

    # We'll track which demand types we've seen flex and inflex for
    demand_types_seen = set()

    for filename in os.listdir(input_dir):
        if not filename.endswith(".csv"):
            continue

        parts = filename.replace(".csv", "").split("_")

        if parts[0] != "demand":
            continue

        # Non-flexible standard demand: demand_T.csv
        # If only two parts: demand_{type}.csv
        if len(parts) == 2:
            demand_type = parts[1]
            dtype_label = type_labels.get(demand_type, demand_type)
            tech_name = f"demand_{demand_type}"

            if tech_name not in techs:
                techs[tech_name] = {
                    "name": f"{dtype_label} electricity demand",
                    "base_tech": "demand",
                    "carrier_in": "electricity"
                }

            filepath = os.path.relpath(os.path.join(input_dir, filename), os.getcwd())

            data_tables[tech_name] = {
                "data": filepath.replace("\\", "/"),
                "rows": "timesteps",
                "columns": "nodes",
                "add_dims": {
                    "techs": tech_name,
                    "parameters": "sink_use_equals"
                }
            }
            demand_types_seen.add(demand_type)

        # Flex or inflex split CSVs with at least 3 parts:
        # e.g. demand_I_flex.csv, demand_I_inflex.csv
        elif len(parts) >= 3:
            demand_type = parts[1]
            flexflag = parts[-1]  # 'flex' or 'inflex'
            dtype_label = type_labels.get(demand_type, demand_type)

            # Tech name simplified to just flex or inflex per demand type
            tech_name = f"demand_{demand_type}_{flexflag}"

            if tech_name not in techs:
                name_desc = f"{dtype_label} electricity demand ({flexflag})"
                techs[tech_name] = {
                    "name": name_desc,
                    "base_tech": "demand",
                    "carrier_in": "electricity"
                }

            param = "sink_use_dsr" if flexflag == "flex" else "sink_use_equals"

            filepath = os.path.relpath(os.path.join(input_dir, filename), os.getcwd())

            data_tables[tech_name] = {
                "data": filepath.replace("\\", "/"),
                "rows": "timesteps",
                "columns": "nodes",
                "add_dims": {
                    "techs": tech_name,
                    "parameters": param
                }
            }
            demand_types_seen.add(demand_type)

        else:
            print(f"Warning: skipping unexpected filename format: {filename}")

    out_dict = {
        "techs": techs,
        "data_tables": data_tables
    }

    with open(output_path, "w") as f:
        yaml.dump(out_dict, f, sort_keys=False, default_flow_style=False, indent=4)

generate_demand_yaml(snakemake.input[0], snakemake.output[0])

# S2 
generate_demand_yaml(snakemake.input[1], snakemake.output[1])

# S5,
scenario = snakemake.params.scenario
s5_input = os.path.join("data/processed/demand/split_by_flex", f"S5_{scenario}")
generate_demand_yaml(s5_input, snakemake.output[2])