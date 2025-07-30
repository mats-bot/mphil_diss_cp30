import pandas as pd
import yaml

df = pd.read_excel(snakemake.input[0])

techs = {}

# All modelled as subsea HVDC
loss_per_100_km = 0.0035 # from https://www.eia.gov/analysis/studies/electricity/hvdctransmission/pdf/transmission.pdf#page=18

for index, row in df.iterrows():
    link_name = f"{row['from_zone']}_{row['to_zone']}"
    flow_cap_max = float(row['capacity'])
    distance = float(row['distance_km'])
    loss = loss_per_100_km * (distance / 100)

    techs[link_name] = {
        "parent": "ac_transmission",
        "link_from": row['from_zone'],
        "link_to": row['to_zone'],
        "flow_cap_max": flow_cap_max,
        "loss": round(loss, 16),
    }