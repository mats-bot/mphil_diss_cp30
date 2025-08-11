import pandas as pd
import yaml

import_prices = pd.read_csv(snakemake.input[0])
import_prices["timesteps"] = pd.to_datetime(import_prices["timesteps"], format="%m/%d/%Y %H:%M").dt.strftime("%Y-%m-%d %H:%M:%S")

# Costs are in  2014 bn€ / 100 GWh so need to convert to 2024 £/MWh 
# https://euro-calliope.readthedocs.io/en/stable/model/customisatio

# https://www.ons.gov.uk/economy/inflationandpriceindices/datasets/consumerpriceinflation
CPI_14_24 = 133.9 / 100.0

# https://www.exchangerates.org.uk/GBP-EUR-spot-exchange-rates-history-2014.html
eur14_gbp14 = 1.2411111

units = 10000 # bn/100GWh -> £/Mwh

import_prices.iloc[:, 1:] = import_prices.iloc[:, 1:] * CPI_14_24 * eur14_gbp14 * units

templates = {
    "import_electricity": {
        "name": "Electricity import",
        "base_tech": "supply",
        "carrier_out": "electricity"
    },
    "export_electricity": {
        "name": "Electricity export",
        "base_tech": "demand",
        "carrier_in": "electricity"
    }
}

# # Subsea HVDC: 0.35%/100km, from https://www.eia.gov/analysis/studies/electricity/hvdctransmission/pdf/transmission.pdf#page=18
loss_per_100 = 0.0035

interconnectors = pd.read_excel(snakemake.input[1])
interconnectors['loss'] = loss_per_100 * (interconnectors['Distance (km)'] / 100)

grouped = interconnectors.groupby(['Country', 'Zone']).apply(
    lambda df: pd.Series({
        'total_capacity': df['Capacity (MW)'].sum(),
        'weighted_loss': (df['loss'] * df['Capacity (MW)']).sum() / df['Capacity (MW)'].sum()
    })
).reset_index()

techs = {}
nodes = {}

for country in grouped['Country'].unique():
    tech_import = f"import_{country.lower()}_electricity"
    tech_export = f"export_{country.lower()}_electricity"
    techs[tech_import] = {"template": "import_electricity"}
    techs[tech_export] = {"template": "export_electricity", "category": "export"}

    if country in import_prices.columns:
        import_prices[tech_import] = import_prices[country]

    zones_for_country = grouped[grouped['Country'] == country]
    for _, row in zones_for_country.iterrows():
        zone = row['Zone']
        cap = float(row['total_capacity'])
        loss = float(row['weighted_loss'])
        if zone not in nodes:
            nodes[zone] = {"techs": {}}
        nodes[zone]["techs"][tech_import] = {
            "flow_cap_min": cap,
            "flow_cap_max": cap,
            "flow_out_eff": 1 - loss,
        }
        nodes[zone]["techs"][tech_export] = {
            "flow_cap_min": cap,
            "flow_cap_max": cap,
            "flow_in_eff": 1 - loss,
        }

data = {
    "import_prices": {
        "data": "data/processed/demand/import_prices.csv",
        "rows": "timesteps",
        "columns": "techs",
        "add_dims": {
            "parameters": "cost_flow_out"
        }
    },
    "export_prices": {
        "data": "data/processed/demand/export_prices.csv",
        "rows": "timesteps",
        "columns": "techs",
        "add_dims": {
            "parameters": "cost_flow_in"
        }
    },
}

import_prices.drop(columns=interconnectors["Country"].unique(), inplace=True)
export_prices = import_prices.copy()
export_prices.iloc[:, 1:] = -export_prices.iloc[:, 1:]

import_prices.to_csv(snakemake.output[0], index=False)
export_prices.to_csv(snakemake.output[1], index=False)

with open(snakemake.output[2], "w") as f:
    yaml.dump({"templates": templates, "techs": techs, "nodes": nodes, "data_tables": data}, f, sort_keys=False, indent=2)
