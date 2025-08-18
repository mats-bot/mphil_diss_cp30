import pandas as pd
import yaml

df = pd.read_excel(snakemake.input[0], sheet_name="ES1", skiprows=9)
year_col = 2030

# Omitted here: marine and other renewables (other is ~1.3GW, marine ~250MW in 2030)
# commented for difference between cp30 and FES24 if need to grab those
def make_tech_map(pathway):
    return {
        "gas_ccgt_new": [
            {"pathway": pathway, "type": "Gas", "subtype": "Gas"},
            {"pathway": pathway, "type": "Gas", "subtype": "Gas CHP"},
        ],
        # "gas_ocgt": [
        #     {"pathway": pathway, "type": "Gas", "subtype": "OCGT"},
        #     {"pathway": pathway, "type": "Gas", "subtype": "Gas Reciprocating Engines"},
        # ],
        "diesel_existing": [
            {"pathway": pathway, "type": "Other Thermal", "subtype": "Fuel Oil"},
            # {"pathway": pathway, "type": "Other Thermal", "subtype": "OCGT"},
            {"pathway": pathway, "type": "Other Thermal", "subtype": "Diesel"},
        ],
        "hydrogen": [
            {"pathway": pathway, "type": "Hydrogen", "subtype": "Hydrogen"},
            {"pathway": pathway, "type": "Hydrogen", "subtype": "Hydrogen CHP"},
        ],
        # "coal": [
        #     {"pathway": pathway, "type": "Coal", "subtype": "Coal"},
        # ],
        "battery": [
            {"pathway": pathway, "type": "Storage", "subtype": "Battery"},
        ],
        # "pumped_hydro": [
        #     {"pathway": pathway, "type": "Storage", "subtype": "Pumped Hydro"},
        # ],
        # "caes": [
        #     {"pathway": pathway, "type": "Storage", "subtype": "Compressed Air"},
        # ],
        # "laes": [
        #     {"pathway": pathway, "type": "Storage", "subtype": "Liquid Air"},
        # ],
        # "ldes": [
        #     {"pathway": pathway, "type": "Storage", "subtype": "LDES"},
        # ],
        "nuclear": [
            {"pathway": pathway, "type": "Nuclear", "subtype": "Nuclear - Large"},
        ],
        "beccs": [
            {"pathway": pathway, "type": "CCS", "subtype": "CCS Biomass"},
        ],
        "gas_ccs": [
            {"pathway": pathway, "type": "CCS", "subtype": "CCS Gas"},
        ],
        "biomass": [
            {"pathway": pathway, "type": "Biomass", "subtype": "Biomass"},
            {"pathway": pathway, "type": "Biomass", "subtype": "Biomass CHP"},
        ],
        "hydro": [
            {"pathway": pathway, "type": "Hydro", "subtype": "Hydro"},
        ],
        "waste": [
            {"pathway": pathway, "type": "Waste", "subtype": "Waste"},
            {"pathway": pathway, "type": "Waste", "subtype": "Waste CHP"},
        ],
        "offshore_wind": [
            {"pathway": pathway, "type": "Offshore Wind", "subtype": "Offshore Wind"},
        ],
        "onshore_wind": [
            {"pathway": pathway, "type": "Onshore Wind", "subtype": "Onshore Wind"},
        ],
        "solar_pv": [
            {"pathway": pathway, "type": "Solar", "subtype": "Solar PV"},
        ],   
    }

def build_yaml(tech_map):
    out = {"techs": {}} 

    # gas offset to account for existing capacity
    # 35 GW - ccgt_existing - ocgt_existing = 35 - 33 = 2
    offsets = {
        "gas_ccgt_new": -33000,
#        "battery": 10000,
    }

    for tech, combos in tech_map.items():
        total_val = 0
        for crit in combos:
            mask = (
                (df["Pathway"] == crit["pathway"]) &
                (df["Type"] == crit["type"]) &
                (df["SubType"] == crit["subtype"]) &
                (df["Variable"] == "Capacity (MW)")
            )
            total_val += df.loc[mask, year_col].sum()

        total_val += offsets.get(tech, 0)

        out["techs"][tech] = {"flow_cap_max_systemwide": float(total_val)}
    return out
    
def save_yaml(data, filename):
    with open(filename, "w") as f:
        yaml.safe_dump(data, f, sort_keys=False, default_flow_style=False)


FFR_full = build_yaml(make_tech_map("Further Flex and Renewables"))
save_yaml(FFR_full, snakemake.output[0])

ND_full = build_yaml(make_tech_map("New Dispatch"))
save_yaml(ND_full, snakemake.output[1])

# csvs of caps for postprocessing
def build_cap_df(full_yaml, pathway):
    base_df = pd.DataFrame([
        {"tech": k, "flow_cap_max_systemwide": v["flow_cap_max_systemwide"]}
        for k, v in full_yaml["techs"].items()
    ])

    extra_types = ["Other Renewables", "Marine"]
    extra_df = (
        df[(df["Type"].isin(extra_types)) & (df["Pathway"] == pathway) & (df["Variable"] == "Capacity (MW)")]
        .groupby("Type")[year_col]
        .sum()
        .reset_index()
        .rename(columns={"Type": "tech", year_col: "flow_cap_max_systemwide"})
    )
    extra_df["tech"] = extra_df["tech"].str.lower().str.replace(" ", "_")

    return pd.concat([base_df, extra_df], ignore_index=True)

ffr_with_extra = build_cap_df(FFR_full, "Further Flex and Renewables")
ffr_with_extra.to_csv(snakemake.output[2], index=False)

nd_with_extra = build_cap_df(ND_full, "New Dispatch")
nd_with_extra.to_csv(snakemake.output[3], index=False)