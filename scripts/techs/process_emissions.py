import pandas as pd
import yaml 

emissions_df = pd.read_csv(snakemake.input[0])

# From IPCC 1.5 aligned scenarios (in line with demand proj) for all of EU
# https://1p5ndc-pathways.climateanalytics.org/countries/european-union/v6/sectors/power
import_emissions = 50 # gCO2 / kWh


# 175 euro/tonne from cp30 assumptions
# https://www.exchangerates.org.uk/EUR-GBP-spot-exchange-rates-history-2024.html
carbon_price = 175 * 1 / 0.8468 # EUR/t * GBP/EURO | 2024

mapping = {
    "beccs": "CCS Biomass",
    "diesel_existing": "Diesel",
    "gas_ocgt_existing": "OCGT",
    "biomass": "Biomass (dedicated)",
    "gas_ccgt_existing": "CCGT",
    "gas_ccgt_new": "CCGT",
    "gas_ccs": "CCS Gas",
    "waste": "EfW Incineration",
    "hydrogen": "Hydrogen",
    "nuclear": "Nuclear",
    "offshore_wind": "Wind Offshore",
    "hydro": "Large Hydro",
    "solar_pv": "Solar Photovoltaics",
    "onshore_wind": "Wind Onshore",
    "battery": "Battery",
    "pumped_hydro": "Pumped Storage Hydroelectric",
    "caes": "Compressed Air",
    "laes": "Liquid Air"
}

import_techs = [
    "import_bel_electricity",
    "import_deu_electricity",
    "import_dnk_electricity",
    "import_fra_electricity",
    "import_irl_electricity",
    "import_nor_electricity",
]

tech_emissions_yaml = {}

for tech, source in mapping.items():
    row = emissions_df.loc[emissions_df["Type"].str.strip() == source,
                           "Emissions Intensity [gCO2/kWh]"]
    if not row.empty:
        g_per_kWh = row.values[0]
        tech_emissions_yaml[tech] = float(g_per_kWh / 1e3)
    else:
        tech_emissions_yaml[tech] = 0.0

for tech in import_techs:
    tech_emissions_yaml[tech] = float(import_emissions / 1e3)

yaml_dict = {"overrides": {"emissions": {"techs": {}}}}

for tech, value in tech_emissions_yaml.items():
    yaml_dict["overrides"]["emissions"]["techs"][tech] = {
        "cost_flow_out": {
            "data": value,
            "index": "emissions",
            "dims": "costs"
        }
    }


with open(snakemake.output[0], "w") as f:
    yaml.dump(yaml_dict, f, sort_keys=False)

