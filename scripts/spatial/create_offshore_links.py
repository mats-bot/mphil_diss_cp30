import pandas as pd
import yaml

df = pd.read_excel(snakemake.input[0])

# # Subsea HVDC: 0.35%/100km, from https://www.eia.gov/analysis/studies/electricity/hvdctransmission/pdf/transmission.pdf#page=18
loss_per_100 = 0.0035

techs = {}

for index, row in df.iterrows():
    link_name = f"{row['Zone 1']}_{row['Zone 2']}_offshore"
    flow_cap_max = float(row['2030 Capacity (MW)'])
    distance = float(row['Distance (km)'])
    loss = 1 - loss_per_100 * (distance / 100)

    techs[link_name] = {
        "template": "subsea_dc_transmission",
        "link_from": row['Zone 1'],
        "link_to": row['Zone 2'],
        "flow_cap_max": flow_cap_max,
        "flow_out_eff": loss,
    }

with open(snakemake.output[0], "w") as f:
    yaml.dump({"techs": techs}, f, sort_keys=False, indent=2)