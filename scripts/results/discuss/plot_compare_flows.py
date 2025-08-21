import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MaxNLocator
from matplotlib.ticker import LogLocator, LogFormatter
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
    "import_bel_electricity": "Interconnectors",
    "import_deu_electricity": "Interconnectors",
    "import_dnk_electricity": "Interconnectors",
    "import_fra_electricity": "Interconnectors",
    "import_irl_electricity": "Interconnectors",
    "import_nor_electricity": "Interconnectors",
    "export_bel_electricity": "Interconnectors",
    "export_deu_electricity": "Interconnectors",
    "export_dnk_electricity": "Interconnectors",
    "export_fra_electricity": "Interconnectors",
    "export_irl_electricity": "Interconnectors",
    "export_nor_electricity": "Interconnectors",
    "demand_C_flex": "Demand",
    "demand_C_inflex": "Demand",
    "demand_D_flex": "Demand",
    "demand_D_inflex": "Demand",
    "demand_E_flex": "Demand",
    "demand_E_inflex": "Demand",
    "demand_H_flex": "Demand",
    "demand_H_inflex": "Demand",
    "demand_I_flex": "Demand",
    "demand_I_inflex": "Demand",
    "demand_R_flex": "Demand",
    "demand_R_inflex": "Demand",
    "demand_T": "Demand",
    "demand_Z": "Demand",
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
    "Interconnectors":"#FFFF99",  
    "Demand":         "#7F7F7F",
}



def collect_flows_grouped(inputs, tech_group_map):
    records = []
    for scen in ["ND", "FFR"]:
        filepath = inputs.get(scen)
        if filepath is None:
            continue
        print(f"Loading {filepath} for {scen}...")
        model = calliope.read_netcdf(filepath)

        for tech, group in tech_group_map.items():
            flow_out_val = 0.0
            flow_in_val = 0.0

            try:
                if tech.startswith("import_"):
                    flow_out_val = model.results.flow_out.loc[{"techs": tech}].sum().item()
                elif tech.startswith("demand_"):
                    flow_in_val = model.results.flow_in.loc[{"techs": tech}].sum().item()
                elif tech.startswith("export_"):
                    flow_in_val = model.results.flow_in.loc[{"techs": tech}].sum().item()
                else:
                    flow_out_val = model.results.flow_out.loc[{"techs": tech}].sum().item()
            except KeyError:
                pass

            records.append({
                "scenario": scen,
                "pathway": "ND" if "ND" in scen else "FFR",
                "tech_group": group,
                "flow_out_GWh": flow_out_val,
                "flow_in_GWh": flow_in_val
            })

    df = pd.DataFrame.from_records(records)
    df_grouped = (
        df.groupby(["scenario", "pathway", "tech_group"], as_index=False)
          .agg({"flow_out_GWh": "sum", "flow_in_GWh": "sum"})
    )
    df_grouped["net_flow_GWh"] = df_grouped["flow_out_GWh"] - df_grouped["flow_in_GWh"]

    interconnectors = df_grouped[df_grouped["tech_group"] == "Interconnectors"]
    print(f"{scen} Calliope interconnectors net flow values:")
    print(interconnectors)


    return df_grouped

def collect_cp30_capacities(filepath, pathway):
    df = pd.read_excel(filepath, sheet_name="ES1", header=9)  
    df = df[["Pathway", "Variable", "SubType", "Type", 2030]]  
    df = df.dropna(subset=[2030])

    tech_order = [
        "Solar PV", "Onshore wind", "Offshore wind", "Hydro",
        "Nuclear", "Biomass", "BECCS", "Waste", "Gas CCS", "Hydrogen",
        "Battery", "LDES", "Unabated gas", 
        "Interconnectors", "Demand"
    ]

    if pathway == "ND":
        df_path = df[df["Pathway"] == "New Dispatch"]
    else:
        df_path = df[df["Pathway"] == "Further Flex and Renewables"]

    records = []
    for tech in tech_order:  
        if tech == "LDES":
            value = df_path[(df_path["SubType"] == "LDES") & (df_path["Variable"] == "Generation (GWh)")][2030].sum()
        elif tech == "Battery":
            value = df_path[(df_path["SubType"] == "Battery") & (df_path["Variable"] == "Generation (GWh)")][2030].sum()
        elif tech == "Unabated gas":
            value = df_path[(df_path["Type"].isin(["Other Thermal", "Gas"])) & (df_path["Variable"] == "Generation (GWh)")][2030].sum()
        elif tech == "Gas CCS":
            value = df_path[(df_path["SubType"] == "CCS Gas") & (df_path["Variable"] == "Generation (GWh)")][2030].sum()
        elif tech == "Hydrogen":
            value = df_path[(df_path["Type"] == "Hydrogen") & (df_path["Variable"] == "Generation (GWh)")][2030].sum()
        elif tech == "Solar PV":
            value = df_path[(df_path["Type"] == "Solar") & (df_path["Variable"] == "Generation (GWh)")][2030].sum()             
        elif tech == "Onshore wind":
            value = df_path[(df_path["Type"] == "Onshore Wind") & (df_path["Variable"] == "Generation (GWh)")][2030].sum()          
        elif tech == "Offshore wind":
            value = df_path[(df_path["Type"] == "Offshore Wind") & (df_path["Variable"] == "Generation (GWh)")][2030].sum()
        elif tech == "Hydro":
            value = df_path[(df_path["Type"] == "Hydro") & (df_path["Variable"] == "Generation (GWh)")][2030].sum()
        elif tech == "Nuclear":
            value = df_path[(df_path["Type"] == "Nuclear") & (df_path["Variable"] == "Generation (GWh)")][2030].sum()
        elif tech == "Biomass":
            value = df_path[(df_path["Type"] == "Biomass") & (df_path["Variable"] == "Generation (GWh)")][2030].sum()
        elif tech == "Waste":
            value = df_path[(df_path["Type"] == "Waste") & (df_path["Variable"] == "Generation (GWh)")][2030].sum()
        elif tech == "BECCS":
            value = df_path[(df_path["SubType"] == "CCS Biomass") & (df_path["Variable"] == "Generation (GWh)")][2030].sum()
        elif tech == "Interconnectors":
            value = abs(df_path[(df_path["SubType"] == "Interconnectors") & (df_path["Variable"] == "Netflows (GWh)")][2030].sum())
        elif tech == "Demand":
            value = df_path[(df_path["Variable"] == "Electrolysis Demand (GWh)") | (df_path["Variable"] == "Demand Total (GWh) - Excluding Electrolysis")][2030].sum()
        else:
            value = 0

        records.append({"pathway": pathway, "tech_group": tech, "net_flow_GWh": value })

    return pd.DataFrame(records)


def plot_capacity_baselines(df_grouped, df_cp30, tech_colors, output_path):
    major_techs = ["Solar PV", "Onshore wind", "Offshore wind", "Battery", "Interconnectors", "Demand"]
    all_tech_order = [
        "Solar PV", "Onshore wind", "Offshore wind", "Hydro",
        "Nuclear", "Biomass", "BECCS", "Waste", "Gas CCS", "Hydrogen",
        "Battery", "LDES", "Unabated gas", "Interconnectors", "Demand"
    ]

    df_grouped["tech_group"] = pd.Categorical(df_grouped["tech_group"], categories=all_tech_order, ordered=True)
    df_grouped = df_grouped.sort_values("tech_group")
    df_grouped["net_flow_GWh"] = df_grouped["net_flow_GWh"].abs().replace(0, 1e0)

    other_techs = [t for t in all_tech_order if t not in major_techs]
    max_other = df_grouped[df_grouped["tech_group"].isin(other_techs)]["net_flow_GWh"].max() * 1.1

    fig, axes = plt.subplots(2, 2, figsize=(16, 12), sharey=False)

    pathways = ["ND", "FFR"]
    for row_idx, tech_set in enumerate([major_techs, other_techs]):
        for col_idx, pathway in enumerate(pathways):
            ax = axes[row_idx, col_idx]
            subset = df_grouped[(df_grouped["pathway"] == pathway) & (df_grouped["tech_group"].isin(tech_set))].dropna(subset=["tech_group"])
            ax.bar(
                subset["tech_group"],
                subset["net_flow_GWh"],
                color=[tech_colors.get(t, "gray") for t in subset["tech_group"]],
                alpha=0.85,
                width=0.6,
                zorder=2
            )

            cp30_subset = df_cp30[(df_cp30["pathway"] == pathway) & (df_cp30["tech_group"].isin(tech_set))]
            y_values = cp30_subset["net_flow_GWh"].replace(0, 1e-1)
            ax.scatter(
                cp30_subset["tech_group"],
                y_values,
                marker='D',
                color='black',
                zorder=3,
                s=60
            )

            title = "New Dispatch" if pathway == "ND" else "Further Flex & Renewables"
            ax.set_title(title + (" major flows" if row_idx==0 else " minor flows"), fontsize=13)
            ax.set_xticks(range(len(subset)))
            ax.set_xticklabels(subset["tech_group"], rotation=45, ha="right", fontsize=13)
            ax.tick_params(axis="y", labelsize=13)
            ax.tick_params(axis="x", labelsize=13)

            ax.grid(axis='y', which="major", linestyle="-", linewidth=0.8, alpha=0.8, zorder=1)
            ax.grid(axis='y', which="minor", linestyle="--", linewidth=0.5, alpha=0.5, zorder=1)
            ax.yaxis.set_minor_locator(AutoMinorLocator(4))
            ax.yaxis.set_major_locator(MaxNLocator(integer=True))
            ax.set_axisbelow(True)

            if row_idx == 1:  
                ax.set_ylim(0, max_other)
            else:
                y_max = subset["net_flow_GWh"].max() * 1.1
                ax.set_ylim(0, y_max)

            if col_idx == 0:
                ax.set_ylabel("Total net flows [GWh]", fontsize=13)
            else:
                ax.set_ylabel("")

    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)



inputs = snakemake.input
output = snakemake.output.caps_comp
df_grouped = collect_flows_grouped(inputs, tech_group_map)
cp30_nd = collect_cp30_capacities(inputs.cp30, "ND")
cp30_ffr = collect_cp30_capacities(inputs.cp30, "FFR")
df_cp30 = pd.concat([cp30_nd, cp30_ffr], ignore_index=True)

plot_capacity_baselines(df_grouped, df_cp30, tech_colors, output)
