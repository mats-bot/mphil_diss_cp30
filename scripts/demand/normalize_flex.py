import pandas as pd

df = pd.read_excel(
    snakemake.input[0],
    sheet_name="EC.05",
    usecols="AB,AG",
    skiprows=8,  
    nrows=6,
    header=None    
)


df.columns = ['2030', '2035']

flex_types = [
    "I&C DSR",
    "I&C Heat",
    "Residential DSR",
    "Residential Heat",
    "Road Transport Smart Charging",
    "Road Transport V2G at peak"
]

df['flex_type'] = flex_types

category_to_types = {
    "I&C DSR": ["C", "I"],
    "Residential DSR": ["R"],
    "I&C Heat": ["C", "I", "H", "D"],
    "Residential Heat": ["R", "H", "D"],
    "Road Transport Smart Charging": ["E"],
    "Road Transport V2G at peak": ["E"]
}


df_melted = df.melt(id_vars='flex_type', var_name='year', value_name='flex_gw')
df_melted['demand_type'] = df_melted['flex_type'].map(category_to_types)

# Flex split proportionally by type since no data on supply by type
df_melted['num_types'] = df_melted['demand_type'].apply(len)
df_melted = df_melted.explode('demand_type')
df_melted['flex_gw'] = df_melted['flex_gw'] / df_melted['num_types']

df_final = df_melted[['flex_type', 'demand_type', 'year', 'flex_gw']]

df_final.to_csv(snakemake.output[0], index=False)