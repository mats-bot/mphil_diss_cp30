import os
import pandas as pd
import yaml

def gen_yaml(base_csv, capacity_csv, solar_cf_csv, onshore_cf_csv, output_yaml_path, weather_year):
    df = pd.read_csv(base_csv, index_col=0)

    renewable_techs = ["Onshore_Wind", "Solar_PV"]

    solar_cf_path = os.path.relpath(solar_cf_csv, os.getcwd()).replace("\\", "/")
    onshore_cf_path = os.path.relpath(onshore_cf_csv, os.getcwd()).replace("\\", "/")

    capacity_df = pd.read_csv(capacity_csv)
    zones = sorted(capacity_df["zone"].unique())

    year = int(weather_year)

    resources = {
        "Solar_PV": f"data/processed/spatial/solar_cf_{year}.csv",
        "Onshore_Wind": f"data/processed/spatial/onshore_cf_{year}.csv",
    }

    tech_zone_capacity = {}
    for tech in capacity_df["CP30 technology"].unique():
        tech_zone_capacity[tech.lower()] = [0.0] * len(zones)

    for _, row in capacity_df.iterrows():
        zone = row["zone"]
        tech = row["CP30 technology"].lower()
        capacity = float(row["Installed Capacity (MWelec)"])
        if tech in tech_zone_capacity and zone in zones:
            idx = zones.index(zone)
            tech_zone_capacity[tech][idx] = capacity

    techs_yaml = {}
    data_tables_yaml = {}

    cf_paths = {"solar_pv": solar_cf_path, "onshore_wind": onshore_cf_path}

    for tech in renewable_techs:
        tech_base = tech.lower()

        techs_yaml[tech_base] = {
            "category": "renewable",
            "cp30_category": "renewable",
            "name": tech,
            "base_tech": "supply",
            "carrier_out": "electricity",
            "source_unit": "per_cap",
            "lifetime": int(df.loc["lifetime", tech]),
            "cost_om_annual": {
                "data": float(df.loc["om_annual", tech]),
                "index": "monetary",
                "dims": ["costs"],
            },
            "cost_flow_out": {
                "data": float(df.loc["om_prod", tech]),
                "index": "monetary",
                "dims": ["costs"],
            },
            "cost_flow_cap": {
                "data": float(df.loc["capex", tech]),
                "index": "monetary",
                "dims": ["costs"],
            },
        }

        data_tables_yaml[f"{tech_base}_cf"] = {
            "data": cf_paths[tech_base],
            "rows": "timesteps",
            "columns": "nodes",
            "add_dims": {"techs": tech_base, "parameters": "source_use_max"},
            'rename_dims': {
                    'time': 'timesteps'
            },
        }

    output_yaml = {"techs": techs_yaml, "data_tables": data_tables_yaml}

    with open(output_yaml_path, "w") as f:
        yaml.dump(output_yaml, f, sort_keys=False)


gen_yaml(
    base_csv=snakemake.input[0],
    capacity_csv=snakemake.input[3],
    solar_cf_csv=snakemake.input[2],
    onshore_cf_csv=snakemake.input[1],
    output_yaml_path=snakemake.output[0],
    weather_year=snakemake.config["weather_year"]
)


# S2 config
gen_yaml(
    base_csv=snakemake.input[0],
    capacity_csv=snakemake.input[3],
    solar_cf_csv=snakemake.input[5],
    onshore_cf_csv=snakemake.input[4],
    output_yaml_path=snakemake.output[1],
    weather_year="2017"
)