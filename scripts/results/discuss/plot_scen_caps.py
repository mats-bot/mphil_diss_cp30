import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import LogLocator
import numpy as np

import calliope

tech_group_map = {
    "solar_pv": "Solar PV",
    "onshore_wind": "Onshore wind",
    "offshore_wind": "Offshore wind",
    "hydro": "Hydro",
    "nuclear": "Nuclear",
    "biomass": "Biomass",
    "beccs": "BECCS",
    "gas_ccs": "Gas CCS",
    "waste": "Waste",
    "hydrogen": "Hydrogen",
    "battery": "Battery",
    "caes": "LDES",
    "laes": "LDES",
    "pumped_hydro": "LDES",
    "gas_ccgt_existing": "Unabated gas",
    "gas_ccgt_new": "Unabated gas",
    "gas_ocgt_existing": "Unabated gas",
    "diesel_existing": "Unabated gas",
    # "import_bel_electricity": "Interconnectors",
    # "import_deu_electricity": "Interconnectors",
    # "import_dnk_electricity": "Interconnectors",
    # "import_fra_electricity": "Interconnectors",
    # "import_irl_electricity": "Interconnectors",
    # "import_nor_electricity": "Interconnectors",
}

tech_colors = {
    "Solar PV":       "#FDBF6F", 
    "Onshore wind":   "#A6D854",  
    "Offshore wind":  "#66C2A5",  
    "Hydro":          "#1F78B4",  
    "Nuclear":        "#E31A1C",  
    "Biomass":        "#33A02C",  
    "Waste":          "#B15928", 
    "Hydrogen":       "#FB9A99",  
    "BECCS":          "#CAB2D6", 
    "Gas CCS":        "#B3DEEB",  
    "Battery":        "#FF69B4",  
    "LDES":           "#984EA3",  
    "Unabated gas":   "#4D4D4D",  
    # "Interconnectors":"#FFFF99",  
}


def collect_capacity_grouped(inputs, tech_group_map):
    records = []
    for scen in ["ND", "FFR"]:
        filepath = inputs.get(scen)
        if filepath is None:
            continue
        model = calliope.read_netcdf(filepath)

        for tech, group in tech_group_map.items():
            try:
                cap_val = model.results.flow_cap.loc[{"techs": tech}].sum().item()
            except KeyError:
                cap_val = 0.0
            records.append({
                "scenario": scen,
                "pathway": "ND" if "ND" in scen else "FFR",
                "tech_group": group,
                "capacity_GW": cap_val
            })
    df = pd.DataFrame.from_records(records)
    return df.groupby(["scenario", "pathway", "tech_group"], as_index=False).sum()


def collect_cp30_capacities(filepath, pathway):
    df = pd.read_excel(filepath, sheet_name="ES1", header=9)  
    df = df[["Pathway", "Variable", "SubType", "Type", 2030]]  
    df = df.dropna(subset=[2030])

    tech_order = [
        "Solar PV", "Onshore wind", "Offshore wind", "Hydro",
        "Nuclear", "Biomass", "BECCS", "Waste", "Gas CCS", "Hydrogen",
        "Battery", "LDES", "Unabated gas", 
    ]

    if pathway == "ND":
        df_path = df[df["Pathway"] == "New Dispatch"]
    else:
        df_path = df[df["Pathway"] == "Further Flex and Renewables"]

    records = []
    for tech in tech_order:  
        if tech == "LDES":
            value = df_path[(df_path["SubType"] == "LDES") & (df_path["Variable"] == "Capacity (MW)")][2030].sum()
        elif tech == "Battery":
            value = df_path[(df_path["SubType"] == "Battery") & (df_path["Variable"] == "Capacity (MW)")][2030].sum()
        elif tech == "Unabated gas":
            value = df_path[(df_path["Type"].isin(["Other Thermal", "Gas"])) & (df_path["Variable"] == "Capacity (MW)")][2030].sum()
        elif tech == "Gas CCS":
            value = df_path[(df_path["SubType"] == "CCS Gas") & (df_path["Variable"] == "Capacity (MW)")][2030].sum()
        elif tech == "Hydrogen":
            value = df_path[(df_path["Type"] == "Hydrogen") & (df_path["Variable"] == "Capacity (MW)")][2030].sum()
        elif tech == "Solar PV":
            value = df_path[(df_path["Type"] == "Solar") & (df_path["Variable"] == "Capacity (MW)")][2030].sum()             
        elif tech == "Onshore wind":
            value = df_path[(df_path["Type"] == "Onshore Wind") & (df_path["Variable"] == "Capacity (MW)")][2030].sum()          
        elif tech == "Offshore wind":
            value = df_path[(df_path["Type"] == "Offshore Wind") & (df_path["Variable"] == "Capacity (MW)")][2030].sum()
        elif tech == "Hydro":
            value = df_path[(df_path["Type"] == "Hydro") & (df_path["Variable"] == "Capacity (MW)")][2030].sum()
        elif tech == "Nuclear":
            value = df_path[(df_path["Type"] == "Nuclear") & (df_path["Variable"] == "Capacity (MW)")][2030].sum()
        elif tech == "Biomass":
            value = df_path[(df_path["Type"] == "Biomass") & (df_path["Variable"] == "Capacity (MW)")][2030].sum()
        elif tech == "Waste":
            value = df_path[(df_path["Type"] == "Waste") & (df_path["Variable"] == "Capacity (MW)")][2030].sum()
        elif tech == "BECCS":
            value = df_path[(df_path["SubType"] == "CCS Biomass") & (df_path["Variable"] == "Capacity (MW)")][2030].sum()
        else:
            value = 0

        records.append({"pathway": pathway, "tech_group": tech, "capacity_GW": value / 1000})

    return pd.DataFrame(records)


def plot_capacity_baselines(df_grouped, df_cp30, tech_colors, out_nd, out_ffr):
    major_techs = ["Solar PV", "Onshore wind", "Offshore wind", "Battery", "Unabated gas"]
    all_tech_order = [
        "Solar PV", "Onshore wind", "Offshore wind", "Hydro",
        "Nuclear", "Biomass", "BECCS", "Waste", "Gas CCS", "Hydrogen",
        "Battery", "LDES", "Unabated gas"
    ]

    df_grouped["tech_group"] = pd.Categorical(df_grouped["tech_group"], categories=all_tech_order, ordered=True)
    df_grouped = df_grouped.sort_values("tech_group")
    df_grouped["capacity_GW"] = df_grouped["capacity_GW"].replace(0, 1e-3)

    other_techs = [t for t in all_tech_order if t not in major_techs]

    tech_sets = [("major", major_techs), ("minor", other_techs)]
    pathways = [("ND", out_nd), ("FFR", out_ffr)]

    for pathway_name, out_file in pathways:
        fig, axes = plt.subplots(1, 2, figsize=(12, 8), sharey=False)

        for ax_idx, (flow_type, tech_set) in enumerate(tech_sets):
            ax = axes[ax_idx]
            subset = df_grouped[
                (df_grouped["pathway"] == pathway_name) & 
                (df_grouped["tech_group"].isin(tech_set))
            ].dropna(subset=["tech_group"])

            ax.bar(
                subset["tech_group"],
                subset["capacity_GW"],
                color=[tech_colors.get(t, "gray") for t in subset["tech_group"]],
                alpha=0.85,
                width=0.3,
                zorder=2
            )

            cp30_subset = df_cp30[
                (df_cp30["pathway"] == pathway_name) & 
                (df_cp30["tech_group"].isin(tech_set))
            ].copy()
            cp30_subset["capacity_GW"] = cp30_subset["capacity_GW"].replace(0, 1e-3)

            ax.scatter(
                cp30_subset["tech_group"],
                cp30_subset["capacity_GW"],
                marker="D",
                color="black",
                zorder=3,
                s=80
            )

            title = "New Dispatch" if pathway_name == "ND" else "Further Flex & Renewables"
            ax.set_title(f"{title} – {flow_type} capacities", fontsize=16)
            ax.set_xticks(range(len(subset)))
            ax.set_xticklabels(subset["tech_group"], rotation=45, ha="right", fontsize=16)
            ax.tick_params(axis="y", labelsize=16)
            ax.tick_params(axis="x", labelsize=16)

            ax.grid(axis="y", which="major", linestyle="-", linewidth=0.8, alpha=0.8, zorder=1)
            ax.grid(axis="y", which="minor", linestyle="--", linewidth=0.5, alpha=0.5, zorder=1)
            ax.yaxis.set_minor_locator(AutoMinorLocator(5))

            combined_max = pd.concat([subset["capacity_GW"], cp30_subset["capacity_GW"]]).max() * 1.1
            ax.set_ylim(0, combined_max)
            ax.set_ylabel("Installed capacity [GW]", fontsize=16)

        fig.tight_layout()
        fig.savefig(out_file, dpi=300)
        plt.close(fig)


inputs = snakemake.input
outputs = snakemake.output

df_grouped = collect_capacity_grouped(inputs, tech_group_map)

cp30_nd = collect_cp30_capacities(inputs.cp30, "ND")
cp30_ffr = collect_cp30_capacities(inputs.cp30, "FFR")
df_cp30 = pd.concat([cp30_nd, cp30_ffr], ignore_index=True)

for scen in [k for k in inputs.keys() if k != "cp30"]:
    df_sub = df_grouped[df_grouped["scenario"] == scen]
    pathway = "ND" if "ND" in scen else "FFR"
    cp30_sub = df_cp30[df_cp30["pathway"] == pathway]
    plot_capacity_baselines(
        df_sub,
        cp30_sub,
        tech_colors,
        outputs[scen]
    )