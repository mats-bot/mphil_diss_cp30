import pandas as pd
import yaml

# Note that losses on the transmission lines are calculated directly from the distances of each link
# in other processing steps for compatibility with Calliope. The values used to do so are:

# Onshore HVAC: median value of 0.65%/100km from https://www.nationalgrid.com/sites/default/files/documents/13784-High%20Voltage%20Direct%20Current%20Electricity%20%E2%80%93%20technical%20information.pdf
# Subsea HVDC: 0.35%/100km, from https://www.eia.gov/analysis/studies/electricity/hvdctransmission/pdf/transmission.pdf#page=18
# values for HVAC in 2nd source aligned with 1st source.

# Costs for the transmission system are neglected. This is because the capital cost estimations
# for upgrades/new projects are well defined by NESO and use proprietary data, and have much more
# information as to the cost of the new projects including their acceleration. Additionally, the 
# operational/maintenance costs are negligible (~0.01 p/kWh for a 1000 MW line) compared to 
# generation as per figures from https://www.pure.ed.ac.uk/ws/portalfiles/portal/21980985/Grid_Carbon_Footprint_Paper.pdf#:~:text=transmission%20network%20are%20around%2010,mix%20is%20largely%20outside%20the

techs = {
    "subsea_dc_transmission": {
        "name": "subsea_dc_transmission",
        "base_tech": "transmission",
        "carrier_in": "electricity",
        "carrier_out": "electricity",
        "link_to": "null",
        "link_from": "null"
    },
    "ac_transmission": {
        "name": "ac_transmission",
        "base_tech": "transmission",
        "carrier_in": "electricity",
        "carrier_out": "electricity",
        "link_to": "null",
        "link_from": "null"
    }
}

with open(snakemake.output[0], "w") as f:
    yaml.dump({"techs": techs}, f, sort_keys=False)


