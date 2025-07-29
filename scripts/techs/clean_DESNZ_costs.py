import pandas as pd
import numpy as np

DESNZ_df = pd.read_excel(snakemake.input[0], sheet_name="Technical and cost assumptions", header=0)

# Sort for techs to keep
technologies_to_include = [
    "Large-scale Solar", # REPD database minimum capacity is 150 kW
    "Onshore Windr", 
    "Offshore Wind", # assumed to be same for floating
    "CCGT H Class", # Assume all CCGT (natural/gas) are H Class
    "CCGT CHP mode", # unsure how to define CHP profiles but may be useful then
    "OCGT 100MW 500 hr", # most single cycles <100 MW or not much over and not expected to run 
                        # year-round, mapped conventioal steam gas to this as well
    "Recip Diesel 500 hrs", # mapped to single cycle disel/gas oil, could also do 90/2000 hrs
    "EfW",
    "EfW CHP",
    "Dedicated Biomass",
    "Biomass CHP",
    "Hydro Large Storage", # (not pumped storage) no option for smaller so used for all hydro
    "Hydrogen CCGT (FOAK)" # NOAK data only available for 2040
]
# Also available: tidal stream, wave, ACT, AD, Landfill, Sewage, Geothermal

# missing techs (from intial mapping):
# nuclear
# coal (conventional steam, only 1 plant)
# LDES (all types)
# Batteries
# BECCS
# CCS gas
 
# Keep only data for projects commissioning in 2025
cols_to_keep = DESNZ_df.columns[:3].tolist()
matching_cols = [col for col in DESNZ_df.columns[3:] if col in technologies_to_include]
final_cols = cols_to_keep + matching_cols
DESNZ_df = DESNZ_df[final_cols]


# fill in 1st column so can call variable names
DESNZ_df.iloc[:, 0].replace('', pd.NA, inplace=True)
DESNZ_df.iloc[:, 0] = DESNZ_df.iloc[:, 0].fillna(method='ffill')

# keep only the medium estimates
DESNZ_df = DESNZ_df[~DESNZ_df.iloc[:, 1].str.contains('Low|High', na=False)]


parameters = ['Pre-development cost', 'Construction cost', 'Infrastructure', 'Reference plant size', 
              'Fixed O&M', 'Insurance', 'Connection and use of system charges', 'Variable O&M', 
              'Average fuel efficiency', 'Operating lifetime']

def get_param_values(df, param):
    row = df[df.iloc[:, 0] == param]
    if row.empty:
        return pd.Series(np.nan, index=df.columns[3:])
    else:
        return row.iloc[0, 3:]

pre_dev_cost = get_param_values(DESNZ_df, 'Pre-development cost').astype(float)
construction_cost = get_param_values(DESNZ_df, 'Construction cost').astype(float)
infrastructure = get_param_values(DESNZ_df, 'Infrastructure').astype(float)
plant_size = get_param_values(DESNZ_df, 'Reference plant size').astype(float)

fixed_om = get_param_values(DESNZ_df, 'Fixed O&M').astype(float)
insurance = get_param_values(DESNZ_df, 'Insurance').astype(float)
connection = get_param_values(DESNZ_df, 'Connection and use of system charges').astype(float)
variable_om = get_param_values(DESNZ_df, 'Variable O&M').astype(float)

fuel_efficiency = get_param_values(DESNZ_df, 'Average fuel efficiency').astype(float)
lifetime_years = get_param_values(DESNZ_df, 'Operating lifetime').astype(float)



capex = pre_dev_cost + construction_cost + (infrastructure * 1000 / (plant_size * 1000)) * 1000 # £/MW
om_annual = (fixed_om + insurance + connection)  # £/MW/year
om_prod = variable_om # £/MWh 
efficiency = fuel_efficiency  # fraction
lifetime = lifetime_years  # years

calc_params = {
    'capex': capex,
    'om_annual': om_annual,
    'om_prod': om_prod,
    'efficiency': efficiency,
    'lifetime': lifetime
}

new_rows = []
for param_name, values in calc_params.items():
    original_row = DESNZ_df[DESNZ_df.iloc[:, 0] == param_name]
    if not original_row.empty:
        first_3_cols = original_row.iloc[0, :3].values
    else:
        first_3_cols = [param_name] + [np.nan]*2

    new_row = list(first_3_cols) + values.tolist()
    new_rows.append(new_row)

new_df = pd.DataFrame(new_rows, columns=DESNZ_df.columns)
new_df.drop(DESNZ_df.columns[[1, 2]], axis=1, inplace=True)

# rename techs to match ones from capacities
Cp30_techs = {
    "Large-scale Solar": "Solar_PV",
    "Onshore Windr": "Onshore_Wind",
    "Offshore Wind": "Offshore_Wind",
    "CCGT H Class": "Gas_CCGT",
    "CCGT CHP mode": "Gas_CCGT_CHP",
    "OCGT 100MW 500 hr": "Gas_OCGT",
    "Recip Diesel 500 hrs": "Diesel",
    "EfW": "Waste",
    "EfW CHP": "Waste_CHP",
    "Dedicated Biomass": "Biomass",
    "Biomass CHP": "Biomass_CHP",
    "Hydro Large Storage": "Hydro",
    "Hydrogen CCGT (FOAK)": "Hydrogen"
}

new_df.rename(columns=Cp30_techs, inplace=True)



new_df.to_csv(snakemake.output[0], index=False)