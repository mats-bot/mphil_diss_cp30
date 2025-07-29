import pandas as pd

df = pd.read_csv(snakemake.input[0], index_col=0)


## Variables to add


# CPI to rebase all costs to £2024 per C30 assumptions (yearly avgs)
# https://www.ons.gov.uk/economy/inflationandpriceindices/datasets/consumerpriceinflation
c14 = 133.9/100.0
c15 = 133.9/100.0
c19 = 133.9/107.8
c21 = 133.9/111.6
c22 = 133.9/121.7

# assuming 30 MW (avg value REPD) and 4h charge/discharge time (common)
# eff, loss, life from https://www.sciencedirect.com/science/article/pii/B9780444626165000164
battery_eff = 0.95 # frac
battery_power_cap = 30  # MW
battery_energy_cap_max= 120  # MWh
battery_life = 25 # yrs, own assumption
battery_storage_loss = 0.02 / (7*24) # frac per hour

# costs from https://www.sciencedirect.com/science/article/pii/S0306261916312740#s0015 2015 euro, mean values
# 2015 costs 
battery_capital_cost = 65 / 0.7263 * c15 * 1000 # euro/MW * pound/euro|2015 https://www.exchangerates.org.uk/EUR-GBP-spot-exchange-rates-history-2015.html
battery_om_prod = 420 * 1 / 0.7263 * c15 * 1000 # euro/MWh * pound/euro|2015
battery_om_annual = battery_capital_cost * 0.02 / 20 # euro/MW * 1/yrs

# pumped hydro  technical from https://www.sciencedirect.com/science/article/pii/B9780444626165000073
pumped_hydro_eff =  0.825 # frac
pumped_hydro_power_cap = 700 # MW 
pumped_hydro_energy_cap_max = pumped_hydro_power_cap * 8 # MWh, 700 is avg pumped hydro cap, assumed 8 hrs but large variation since one 1.7 GW and rest <0.5 GW
pumped_hydro_life =  80 # yrs
pumped_hydro_storage_loss = 0.025 / 2 / 24 # frac per hour

# costs from https://www.sciencedirect.com/science/article/pii/S0306261916312740#s0015
# 2015 costs 
pumped_hydro_capital_cost = 750 / 0.7263 * c15 * 1000 # euro/MW * pound/euro|2015  https://www.exchangerates.org.uk/EUR-GBP-spot-exchange-rates-history-2015.html
pumped_hydro_energy_storage_cost = 15 * 1 /0.7263 * c15 * 1000 # euro/MWh * pound/euro|2015
pumped_hydro_om_fixed = 11 / 80 * 1 / 0.7263 * c15 * 1000 # euro/MW/yr * pound/euro|2015 
pumped_hydro_om_prod = 0.0005 * 1 / 0.7263 * c15 * 1000 # euro/MWh * pound/euro|2015


# CAES technical from https://www.sciencedirect.com/science/article/pii/B9780444626165000073
# assume 5 MW from REPD and 6 hrs from https://assets.publishing.service.gov.uk/media/659bde4dd7737c000ef3351a/long-duration-electricity-storage-policy-framework-consultation.pdf#:~:text=Pumped%20hydro%20storage%204%20hours,of%20projects%20that%20provide%20the
CAES_eff = 0.75 # frac
CAES_power_cap = 5 # MW
CAES_energy_cap_max = CAES_power_cap * 6 # MWh
CAES_life = 25 # yrs
CAES_storage_loss = 0.75 / 24 # frac per hour

# costs from https://www.sciencedirect.com/science/article/pii/S0306261916312740#s0015
# 2015 costs
CAES_capital_cost = 795 * 1.25 * 1 / 0.7263 * c15 * 1000 # euro/MW * pound/euro|2015  https://www.exchangerates.org.uk/EUR-GBP-spot-exchange-rates-history-2015.html
CAES_energy_storage_cost = 25 * 1 / 0.7263 * c15 * 1000 # euro/MWh * pound/euro|2015
CAES_om_fixed = 9 / 35 * 1 / 0.7263 * c15 * 1000 # euro/MW/yr * pound/euro|2015 
CAES_om_prod = 0.0033 * 1 / 0.7263 * c15 * 1000 # euro/MWh * pound/euro|2015


# LAES 14 hrs, 0.55 eff from https://assets.publishing.service.gov.uk/media/659bde4dd7737c000ef3351a/long-duration-electricity-storage-policy-framework-consultation.pdf#:~:text=Pumped%20hydro%20storage%204%20hours,of%20projects%20that%20provide%20the
# assumed 50 MW from REPD
# lifespan from https://assets.publishing.service.gov.uk/media/659bde4dd7737c000ef3351a/long-duration-electricity-storage-policy-framework-consultation.pdf
# storage loss from https://www.sciencedirect.com/science/article/pii/S1364032120308571
LAES_eff = 0.55 # frac
LAES_power_cap = 50 # MW
LAES_energy_cap_max = LAES_power_cap * 14 # MWh
LAES_life = 35 # yrs
LAES_storage_loss = 0.15 / 24 # frac per hr

# LAES capital costs, median values from https://www.sciencedirect.com/science/article/pii/S0360544219323758?via%3Dihub
# OPEX costs from https://www.sciencedirect.com/science/article/pii/S0360544220303820#bib30
# 2019 costs 
LAES_capital_cost = 1365 * c19 * 1000 # £/MW
LAES_storage_cap = 330 * c19 * 1000 # £/MWh
# 2015 costs 
LAES_om_fixed = 11.2 * 1 / 0.7263 * c15 * 1000 # euro/MW/yr * pound/euro|2015 #https://www.exchangerates.org.uk/EUR-GBP-spot-exchange-rates-history-2015.html
LAES_om_con = 30 * 0.00264 * 1 / 0.7263 * c15 * 1000 # euro/MWh * pound/euro|2015

storage_data = {
    "technology": ["battery", "pumped_hydro", "caes", "laes"],
    
    "energy_eff": [battery_eff, pumped_hydro_eff, CAES_eff, LAES_eff],
    "storage_loss": [battery_storage_loss, pumped_hydro_storage_loss, CAES_storage_loss, LAES_storage_loss],
    "charge_rate": [
        battery_power_cap / battery_energy_cap_max,
        pumped_hydro_power_cap / pumped_hydro_energy_cap_max,
        CAES_power_cap / CAES_energy_cap_max,
        LAES_power_cap / LAES_energy_cap_max
    ],
    "energy_cap_equals": [battery_power_cap, pumped_hydro_power_cap, CAES_power_cap, LAES_power_cap],
    "storage_cap_equals": [battery_energy_cap_max, pumped_hydro_energy_cap_max, CAES_energy_cap_max, LAES_energy_cap_max],
    "lifetime": [battery_life, pumped_hydro_life, CAES_life, LAES_life],
    
    "cost_energy_cap": [battery_capital_cost, pumped_hydro_capital_cost, CAES_capital_cost, LAES_capital_cost],
    "cost_storage_cap": [None, pumped_hydro_energy_storage_cost, CAES_energy_storage_cost, LAES_storage_cap],
    "om_annual": [battery_om_annual, pumped_hydro_om_fixed, CAES_om_fixed, LAES_om_fixed],
    "om_prod": [battery_om_prod, pumped_hydro_om_prod, CAES_om_prod, None],
    "om_con": [None, None, None, LAES_om_con],
}

storage_df = pd.DataFrame(storage_data)

storage_df.to_csv(snakemake.output[0], index=None)



## Modification of generation (DESNZ) df 

# Convert existing costs to GBP2024 from GBP2021
monetary_rows = ['capex', 'om_annual', 'om_prod']
df.loc[monetary_rows] = df.loc[monetary_rows] * c21


# Adding in column for undefined coal and nuclear (and BECCS in future)
# coal values relatively unimportant since production cease, just need to make tech unattractive so 
# using gas ccgt as ref
df['Coal'] = df['Gas_CCGT']
df.loc['capex', 'Coal'] *= 1.5
df.loc['om_annual', 'Coal'] *= 1.2
df.loc['om_prod', 'Coal'] *= 2
df.loc['efficiency', 'Coal'] = 0.25 # own assumption

# Adding nuclear (only need operational costs since no new can be built)
# costs from https://www.gov.uk/government/publications/beis-electricity-generation-costs-november-2016 (both report and spreadhseet)
# while for FOAK PWR commissioning 2025, is CP30 ref. took medium estimates.
# assumed insurance meant to be per MW and not MWh since value otherwie uncreasonable
# costs in 2014 GBP
nuclear_capex = (240 + 4100 + 11500000/3300000) * c14  * 1000 #  £/MW includes insurance
nuclear_om_annual = (72900+500+10000) * c14 # £/MW/yr  includes connection/use of system and insurance
nuclear_om_prod = (5 + 5 + 2) * c14 # £/MWh includes fuel and decommissioning  

df['Nuclear'] = None
df.loc['capex', 'Nuclear'] = nuclear_capex
df.loc['om_annual', 'Nuclear'] = nuclear_om_annual
df.loc['om_prod', 'Nuclear'] = nuclear_om_prod
df.loc['efficiency', 'Nuclear'] = 1.0
df.loc['lifetime', 'Nuclear'] = 60


# CCS cost (transport and storage) from CP30 - provide 78/tco2 as hihg cost for sensitivity
co2ts = 18 * c15 # £/tCO2, 2015 cost https://ukerc.rl.ac.uk/publications/technical_report/CCS_CC1026_14.pdf
# Converting to add to BECCS and Gas CCS costs 

# sawmill residue CCS at 88% capture rate: https://www.sciencedirect.com/science/article/pii/S0961953421002002
BECCS_co2 = 1067 * co2ts/1000  # kgCO2/MWh * £/tCO2 * tCO2/kgCO2  = £/MWh

# CCGT + CCS source for UK, but considers co2e 
gas_ccs_co2 = 676.8 * co2ts/1000  # kgCO2/MWh * £/tCO2 * tCO2/kgCO2 = £/MWh



# Adding BECCS, assuming 2025 commissioning and FOAK, medium estimates
# costs from spreadsheet https://www.gov.uk/government/publications/beis-electricity-generation-costs-2020
# costs in 2018 basis
BECCS_capex = (100 + 3400 + 29600/396) * c14 * 1000 #  £/MW
BECCS_om_annual = (160400+29600+4100) * c14 # £/MW/yr  includes connection/use of system and insurance
BECCS_om_prod = 4 * c14 + BECCS_co2 # £/MWh includes fuel and decommissioning

df['BECCS'] = None
df.loc['capex', 'BECCS'] = BECCS_capex
df.loc['om_annual', 'BECCS'] = BECCS_om_annual
df.loc['om_prod', 'BECCS'] = BECCS_om_prod
df.loc['efficiency', 'BECCS'] = 1.0 # need to change if assign biomass cost
df.loc['lifetime', 'BECCS'] = 25 


# Adding Gas + CCS, assuming 2025 commissioning and FOAK, medium estimates
# costs from spreadsheet https://www.gov.uk/government/publications/beis-electricity-generation-costs-2020
# costs in 2018 basis
Gas_CCS_capex = (10 + 1500 + 16100/1056) * c14 * 1000 #  £/MW
Gas_CCS_om_annual = (25800+16100+3500) * c14 # £/MW/yr  includes connection/use of system and insurance
Gas_CCS_om_prod = 5 * c14 + BECCS_co2 # £/MWh includes fuel and decommissioning

df['Gas_CCS'] = None
df.loc['capex', 'Gas_CCS'] = Gas_CCS_capex
df.loc['om_annual', 'Gas_CCS'] = Gas_CCS_om_annual
df.loc['om_prod', 'Gas_CCS'] = Gas_CCS_om_prod
df.loc['efficiency', 'Gas_CCS'] = 0.47 
df.loc['lifetime', 'Gas_CCS'] = 25 


# Gas cost said to be ~100 p/therm in CP30 assumptions, close to Oct 23/24 values. took 23 yr avg.
# 2022 cost
gas_cost = 102/100  * 1/29.31 * c22 * 1000 # GBP/therm * therm/MWh https://www.metric-conversions.org/energy-and-power/therms-uk-to-kilowatt-hours.htm
hydrogen_cost = 1.2 * gas_cost

# Own assumption for now. 
waste_cost = 0
biomass_cost = 0

# giving coal random high values since cannot be used in 2030
coal_cost = 1.5 * gas_cost

# own assumption for diesel
diesel_fuel_oil_cost = gas_cost

# Overwriting diesel om_annual since in data is negative (receives subsidies to operate), set to zero for now
df.loc["om_annual", "Diesel"] = 0

# Overwriting efficiency=0 for biomasss (since no cost)
df.loc["efficiency", "Biomass"] = 1.0


fuel_costs = {
    "Gas_CCGT": gas_cost,
    "Gas_CCGT_CHP": gas_cost,
    "Gas_OCGT": gas_cost,
    "Gas_CCS": gas_cost,
    "Diesel": diesel_fuel_oil_cost,
    "Coal": coal_cost,
    "Hydrogen": hydrogen_cost, # Assumption is 1.2x gas price per CP30 assumptions
    "Waste": waste_cost,
    "Waste_CHP": waste_cost,
    "Biomass_CHP": biomass_cost,
    "Biomass": biomass_cost
}

fuel_cost_row = pd.Series({tech: fuel_costs.get(tech, 0) for tech in df.columns}, name="fuel_cost")
df.loc["fuel_cost"] = fuel_cost_row

df.to_csv(snakemake.output[1])

