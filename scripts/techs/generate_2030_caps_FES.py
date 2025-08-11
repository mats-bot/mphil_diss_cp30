import pandas as pd
import yaml

df = pd.read_excel(snakemake.input[0], sheet_name="ES1", skiprows=9)
year_col = 2030

# Omitted hereL marine and other renewables (other is ~1.3GW, marine ~250MW in 2030)
tech_map = {
    "gas_ccgt": [
        {"pathway": "Holistic Transition", "type": "Gas", "subtype": "CCGT"},
    ],
    "gas_ocgt": [
        {"pathway": "Holistic Transition", "type": "OCGT", "subtype": "existing"},
        {"pathway": "Holistic Transition", "type": "Gas", "subtype": "Gas Reciprocating"},
    ],
    "diesel": [
        {"pathway": "Holistic Transition", "type": "Other Thermal", "subtype": "Fuel Oil"},
        {"pathway": "Holistic Transition", "type": "Other Thermal", "subtype": "OCGT"},
        {"pathway": "Holistic Transition", "type": "Other Thermal", "subtype": "Diesel"},
    ],
    "hydrogen": [
        {"pathway": "Net Zero Pathway", "type": "Hydrogen", "subtype": "Hydrogen"},
        {"pathway": "Net Zero Pathway", "type": "Hydrogen", "subtype": "Hydrogen CHP"},
    ],
    # "coal": [
    #     {"pathway": "Holistic Transition", "type": "Coal", "subtype": "Coal"},
    # ],
    "battery": [
        {"pathway": "Holistic Transition", "type": "Storage", "subtype": "Battery"},
    ],
    "pumped_hydro": [
        {"pathway": "Holistic Transition", "type": "Storage", "subtype": "Pumped Hydro"},
    ],
    "caes": [
        {"pathway": "Holistic Transition", "type": "Storage", "subtype": "Compressed Air"},
    ],
    "laes": [
        {"pathway": "Holistic Transition", "type": "Storage", "subtype": "Liquid Air"},
    ],
    "nuclear": [
        {"pathway": "Holistic Transition", "type": "Low Carbon", "subtype": "Nuclear - Large"},
    ],
    "beccs": [
        {"pathway": "Holistic Transition", "type": "CCS", "subtype": "CCS Biomass"},
    ],
    "gas_ccs": [
        {"pathway": "Holistic Transition", "type": "CCS", "subtype": "CCS Gas"},
    ],
    "biomass": [
        {"pathway": "Holistic Transition", "type": "Biomass", "subtype": "Biomass"},
        {"pathway": "Holistic Transition", "type": "Biomass", "subtype": "Biomass CHP"},
    ],
    "hydro": [
        {"pathway": "Holistic Transition", "type": "Hydro", "subtype": "Hydro"},
    ],
    "waste": [
        {"pathway": "Holistic Transition", "type": "Waste", "subtype": "Waste"},
        {"pathway": "Holistic Transition", "type": "Waste", "subtype": "Waste CHP"},
    ],
    "offshore_wind": [
        {"pathway": "Holistic Transition", "type": "Offshore Wind", "subtype": "Offshore Wind"},
    ],
    "onshore_wind": [
        {"pathway": "Holistic Transition", "type": "Onshore Wind", "subtype": "Onshore Wind"},
    ],
    "solar_pv": [
        {"pathway": "Holistic Transition", "type": "Solar PV", "subtype": "Solar PV"},
    ],   
}

def build_yaml(tech_map, filename):
    out = {"techs": {}} 
    for tech, combos in tech_map.items():
        total_val = 0
        for crit in combos:
            mask = (
                (df["Pathway"] == crit["pathway"]) &
                (df["Type"] == crit["type"]) &
                (df["SubType"] == crit["subtype"])
            )
            total_val += df.loc[mask, year_col].sum()
        out[tech] = {"flow_cap_max_systemwide": float(total_val)}
    
    with open(filename, "w") as f:
        yaml.safe_dump(out, f, sort_keys=False)

build_yaml(tech_map, snakemake.output[0])

tech_map_h2 = tech_map.copy()
tech_map_h2 = {k: v[:] for k, v in tech_map.items()}  # deep-ish copy of lists
tech_map_h2["battery"] = [
    {"pathway": "Hydrogen Evolution", "type": "Storage", "subtype": "Battery"},
]
tech_map_h2["offshore_wind"] = [
    {"pathway": "Hydrogen Evolution", "type": "Offshore Wind", "subtype": "Offshore Wind"},
]

build_yaml(tech_map_h2, snakemake.output[1])