import pandas as pd
import yaml
import re

costs_df  = pd.read_csv(snakemake.input[0], index_col=0)
projects_df = pd.read_csv(snakemake.input[1])
cfs = pd.read_csv(snakemake.input[2])

# Get min and max capacity per zone (max = full queue, min = 2023 operational)
projects_df["Operational"] = pd.to_datetime(projects_df["Operational"], errors="coerce")
cutoff_date = pd.Timestamp("2024-01-01")

df_pre2024 = projects_df[projects_df['Operational'] < cutoff_date]
flow_cap_min = df_pre2024.groupby('tzone')['Installed Capacity (MWelec)'].sum().rename('flow_cap_min')
flow_cap_max = projects_df.groupby('tzone')['Installed Capacity (MWelec)'].sum().rename('flow_cap_max')
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

cf_zone_avg_T.to_csv(snakemake.output[2])


# Generate yaml file
cost_flow_cap = float(costs_df.at["capex", "Offshore_Wind"])
cost_om_annual = float(costs_df.at["om_annual", "Offshore_Wind"])
cost_flow_out = float(costs_df.at["om_prod", "Offshore_Wind"])
lifetime = int(float(costs_df.at["lifetime", "Offshore_Wind"]))


nodes = {}
for zone in capacities_df.columns:
    nodes[zone] = {
        'techs': {
            'wind_offshore': {
                'flow_cap_min': float(capacities_df.loc['flow_cap_min', zone]),
                'flow_cap_max': float(capacities_df.loc['flow_cap_max', zone])
            }
        }
    }


yaml_data = {
    'techs': {
        'offshore_wind': {
            "category": "renewable",
            "cp30_category": "renewable",
            "base_tech": "supply",
            "name": "offshore_wind",
            "carrier_out": "electricity",
            "resource_unit": "per_unit",
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
        'capacityfactors_wind_offshore': {
            'data': 'data/processed/spatial/offshore_cfs_2030_aggregated.csv',
            'rows': 'time',
            'columns': 'nodes',
            'add_dims': {
                'techs': 'wind_offshore',
                'parameters': 'source_use_max'
            }
        }
    },
    'nodes': nodes
} 


with open(snakemake.output[0], 'w') as f:
    yaml.dump(yaml_data, f, sort_keys=False)

