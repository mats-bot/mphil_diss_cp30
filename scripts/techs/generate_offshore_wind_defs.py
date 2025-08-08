import re

import pandas as pd
import yaml

costs_df = pd.read_csv(snakemake.input[0], index_col=0)
projects_df = pd.read_csv(snakemake.input[1])
cfs = pd.read_csv(snakemake.input[2])

# Get min and max capacity per zone (max = full queue, min = 2023 operational)
zones = [f"z{i}" for i in range(1, 18)]

projects_df["Operational"] = pd.to_datetime(projects_df["Operational"], errors="coerce")
cutoff_date = pd.Timestamp("2024-01-01")

df_pre2024 = projects_df[projects_df['Operational'] < cutoff_date]
flow_cap_min = df_pre2024.groupby('tzone')['Installed Capacity (MWelec)'].sum().reindex(zones, fill_value=0).rename('flow_cap_min')
flow_cap_max = projects_df.groupby('tzone')['Installed Capacity (MWelec)'].sum().reindex(zones, fill_value=0).rename('flow_cap_max')
capacities_df = pd.concat([flow_cap_min, flow_cap_max], axis=1).T
capacities_df.index = ['flow_cap_min', 'flow_cap_max']

capacities_df.to_csv(snakemake.output[1])


# Avergae capacity factors by zone
time = cfs.iloc[:, 0]
cfs_data = cfs.iloc[:, 1:]

cfs_data_T = cfs_data.T
cfs_data_T.index.name = 'Site Name'

cfs_with_zone = cfs_data_T.merge(projects_df[['Site Name', 'tzone']], left_index=True, right_on='Site Name')
numeric_cols = cfs_with_zone.select_dtypes(include='number').columns
cfs_zone_avg = cfs_with_zone.groupby('tzone')[numeric_cols].mean()

cf_zone_avg_T = cfs_zone_avg.T
cf_zone_avg_T.insert(0, 'time', time)

cf_zone_avg_T.to_csv(snakemake.output[2], index=None)


# Generate yaml file
cost_flow_cap = float(costs_df.at["capex", "Offshore_Wind"])
cost_om_annual = float(costs_df.at["om_annual", "Offshore_Wind"])
cost_flow_out = float(costs_df.at["om_prod", "Offshore_Wind"])
lifetime = int(float(costs_df.at["lifetime", "Offshore_Wind"]))

# nodes = {}
# for zone in capacities_df.columns:
#     nodes[zone] = {
#         'techs': {
#             'wind_offshore': {
#                 'flow_cap_min': float(capacities_df.loc['flow_cap_min', zone]),
#                 'flow_cap_max': float(capacities_df.loc['flow_cap_max', zone])
#             }
#         }
#     }

# BRYNS BIT
# Parent tech
# template = {
#     "offshore_wind": {
#         "category": "renewable",
#         "cp30_category": "renewable",
#         "base_tech": "supply",
#         "name": "offshore_wind",
#         "carrier_out": "electricity",
#         "source_unit": "per_cap",
#         "lifetime": lifetime,
#         # "resource": cfs_path, # TODO: use data table
#         "cost_flow_cap": {"data": capex, "index": "monetary", "dims": ["costs"]},
#         "cost_om_annual": {"data": om_annual, "index": "monetary", "dims": ["costs"]},
#         "cost_flow_out": {"data": om_prod, "index": "monetary", "dims": ["costs"]},
#     }
# }
# techs = {}


# # Func to clean names for calliope
# def sanitize_tech_name(name):
#     name = re.sub(r"\([^)]*\)", "", name)
#     name = re.sub(r"[ \-]", "_", name)
#     name = re.sub(r"[^\w]", "", name)
#     name = name.lower().strip("_")
#     if not name or not name[0].isalpha():
#         name = "t_" + name
#     return name


# cf_columns = pd.read_csv(cfs_path, nrows=1).columns

# # Site specific
# for _, row in projects_df.iterrows():
#     site = row["Site Name"]
#     slug = sanitize_tech_name(site)
#     tech_name = f"offshorewind_{slug}"

#     installed_cap = row["Installed Capacity (MWelec)"]

#     status = row["Development Status (short)"]
#     op_date = row["Operational"]

#     # Match site names explicitly since funky characters
#     matched_column = next(
#         (col for col in cf_columns if col.strip() == site.strip()), None
#     )
#     if matched_column is None:
#         raise ValueError(f"Could not match '{site}' to any column in {cfs_path}")

#     tech = {
#         "template": "offshore_wind",
#         "resource_column": matched_column,
#         "flow_cap_per_unit": installed_cap,
#         "units_max": 1,
#     }

#     # Category 1 override: Operational by end of 2023
#     if pd.notnull(op_date) and op_date < cutoff_date:
#         tech.update(
#             {
#                 "flow_cap_min": installed_cap,
#                 "flow_cap_max": installed_cap,
#                 "cost_flow_cap": 0,  # keep scalar zero here
#             }
#         )
#         # remove keys not needed when energy_cap_equals is used
#         tech.pop("flow_cap_per_unit", None)
#         tech.pop("units_max", None)
#         tech.pop("integer", None)

#     # Category 2 override: Will be built by 2030
#     elif status in {"Operational", "Under Construction", "Awaiting Construction"}:
#         tech["units_min"] = 1

#     # Category 3 (all others): use constraints as is (model decides)
#     techs[tech_name] = tech

yaml_data = {
    'techs': {
        'offshore_wind': {
            "category": "renewable",
            "cp30_category": "renewable",
            "base_tech": "supply",
            "name": "offshore_wind",
            "carrier_out": "electricity",
            "source_unit": "per_cap", 
            "lifetime": lifetime,
            'cost_flow_cap': {
                'data': cost_flow_cap,
                'index': 'monetary',
                'dims': ['costs']
            },
            'cost_om_annual': {
                'data': cost_om_annual,  
                'index': 'monetary',
                'dims': ['costs']
            },
            'cost_flow_out': {
                'data': cost_flow_out,  
                'index': 'monetary',
                'dims': ['costs']
            }
        }    
    },
    'data_tables': {
        'offshore_wind_cf': {
            'data': 'data/processed/spatial/offshore_cfs_2013.csv',
            'rows': 'timesteps',
            'columns': 'nodes',
            'add_dims': {
                'techs': 'offshore_wind',
                'parameters': 'source_use_max',
            },
            'rename_dims': {
                'time': 'timesteps'
            },
        },
        'offshore_wind_capacities': {
            'data': 'data/processed/techs/offshore_wind_projects_aggregated.csv',
            'rows': 'parameters',
            'columns': 'nodes',
            'add_dims': {
                'techs': 'offshore_wind'
                }
        }
    },
} 

# with open(snakemake.output[0], "w") as f:
#     yaml.dump({"techs": techs, "templates": template}, f, sort_keys=False)

with open(snakemake.output[0], 'w') as f:
   yaml.dump(yaml_data, f, sort_keys=False)

