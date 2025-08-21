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


def plot_capacity_baselines(df_grouped, df_cp30, tech_colors, output_path):
    high_value_techs = ["Solar PV", "Onshore wind", "Offshore wind", "Battery", "Unabated gas"]
    other_techs = [t for t in df_grouped["tech_group"].unique() if t not in high_value_techs]

    fig, axes = plt.subplots(2, 2, figsize=(14, 10), sharey=False)

    pathway_list = ["ND", "FFR"]
    max_high = df_grouped[df_grouped["tech_group"].isin(high_value_techs)]["capacity_GW"].max() * 1.1
    max_other = df_grouped[df_grouped["tech_group"].isin(other_techs)]["capacity_GW"].max() * 1.1

    for row_idx, tech_group in enumerate([high_value_techs, other_techs]):
        for col_idx, pathway in enumerate(pathway_list):
            ax = axes[row_idx, col_idx]
            subset = df_grouped[df_grouped["pathway"] == pathway]
            subset = subset[subset["tech_group"].isin(tech_group)]
            subset = subset.dropna(subset=["tech_group"])

            ax.bar(
                subset["tech_group"],
                subset["capacity_GW"].replace(0, 1e-3),
                color=[tech_colors.get(t, "gray") for t in subset["tech_group"]],
                alpha=0.85,
                width=0.6,
                zorder=2
            )

            cp30_subset = df_cp30[df_cp30["pathway"] == pathway]
            cp30_subset = cp30_subset[cp30_subset["tech_group"].isin(tech_group)]
            ax.scatter(
                cp30_subset["tech_group"],
                cp30_subset["capacity_GW"].replace(0, 1e-3),
                marker='D',
                color='black',
                zorder=3,
                s=30
            )

            title = "New Dispatch" if pathway == "ND" else "Further Flex & Renewables"
            ax.set_title(title + (" major capacities" if row_idx==0 else " minor capacities"), fontsize=13)
            ax.set_xticks(range(len(subset)))
            ax.set_xticklabels(subset["tech_group"], rotation=45, ha="right", fontsize=13)
            ax.tick_params(axis="y", labelsize=13)

            ax.grid(axis='y', which="major", linestyle="-", linewidth=0.8, alpha=0.8, zorder=1)
            ax.grid(axis='y', which="minor", linestyle="--", linewidth=0.5, alpha=0.5, zorder=1)
            ax.yaxis.set_minor_locator(AutoMinorLocator(5))

            high_value_techs = ["Solar PV", "Onshore wind", "Offshore wind", "Battery", "Interconnectors", "Demand"]
            other_techs = [t for t in df_grouped["tech_group"].unique() if t not in high_value_techs]

            if row_idx == 0:  
                ax.set_ylim(0, max_high)
            else:  
                ax.set_ylim(0, max_other)

    axes[0, 0].set_ylabel("Installed capacity [GW]", fontsize=13)
    axes[1, 0].set_ylabel("Installed capacity [GW]", fontsize=13)
    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)





inputs = snakemake.input
output = snakemake.output.caps_comp
df_grouped = collect_capacity_grouped(inputs, tech_group_map)
cp30_nd = collect_cp30_capacities(inputs.cp30, "ND")
cp30_ffr = collect_cp30_capacities(inputs.cp30, "FFR")
df_cp30 = pd.concat([cp30_nd, cp30_ffr], ignore_index=True)

plot_capacity_baselines(df_grouped, df_cp30, tech_colors, output)
