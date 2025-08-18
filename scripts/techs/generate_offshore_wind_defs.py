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
# cf_zone_avg_T.loc[:, cf_zone_avg_T.columns != 'time'] *= 100

cf_zone_avg_T.to_csv(snakemake.output[2], index=None)
print("\nStatistics across zones (per time step):")
print(cfs_zone_avg.agg(['mean', 'min', 'max'], axis=1))

print("\nStatistics across zones (per time step):")
print(cfs_zone_avg.agg(['mean', 'min', 'max'], axis=1))


# Generate yaml file
cost_flow_cap = float(costs_df.at["capex", "Offshore_Wind"])
cost_om_annual = float(costs_df.at["om_annual", "Offshore_Wind"])
cost_flow_out = float(costs_df.at["om_prod", "Offshore_Wind"])
lifetime = int(float(costs_df.at["lifetime", "Offshore_Wind"]))


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
            'data': 'data/processed/spatial/offshore_cfs_2013_aggregated.csv',
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


with open(snakemake.output[0], 'w') as f:
   yaml.dump(yaml_data, f, sort_keys=False)



# S2 config

# Avergae capacity factors by zone
S2_cfs = pd.read_csv(snakemake.input[3])
S2_projects_df = pd.read_csv(snakemake.input[1])

S2_time = S2_cfs.iloc[:, 0]
S2_cfs_data = S2_cfs.iloc[:, 1:]

S2_cfs_data_T = S2_cfs_data.T
S2_cfs_data_T.index.name = 'Site Name'

S2_cfs_with_zone = S2_cfs_data_T.merge(S2_projects_df[['Site Name', 'tzone']], left_index=True, right_on='Site Name')
S2_numeric_cols = S2_cfs_with_zone.select_dtypes(include='number').columns
S2_cfs_zone_avg = S2_cfs_with_zone.groupby('tzone')[S2_numeric_cols].mean()

S2_cf_zone_avg_T = S2_cfs_zone_avg.T
S2_cf_zone_avg_T.insert(0, 'time', S2_time)
S2_cf_zone_avg_T.loc[:, S2_cf_zone_avg_T.columns != 'time'] *= 100

S2_cf_zone_avg_T.to_csv(snakemake.output[3], index=None)



S2_yaml_data = {
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
            'data': 'data/processed/spatial/offshore_cfs_2017_aggregated.csv',
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

with open(snakemake.output[4], 'w') as f:
   yaml.dump(S2_yaml_data, f, sort_keys=False)