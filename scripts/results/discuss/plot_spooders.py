import calliope
import pandas as pd
import numpy as np
import yaml
import matplotlib.pyplot as plt
import seaborn as sns

tech_group_map = {
    "solar_pv": "RES", "onshore_wind": "RES", "offshore_wind": "RES",
    "hydro": "LowC", "nuclear": "LowC", "biomass": "LowC", "beccs": "LowC",
    "gas_ccs": "LowC", "waste": "LowC", "hydrogen": "LowC",
    "battery": "Stor", "caes": "Stor", "laes": "Stor", "pumped_hydro": "Stor",
    "gas_ccgt_existing": "Fossil", "gas_ccgt_new": "Fossil", 
    "gas_ocgt_existing": "Fossil", "diesel_existing": "Fossil",
    "import_bel_electricity": "Import", "import_deu_electricity": "Import",
    "import_dnk_electricity": "Import", "import_fra_electricity": "Import",
    "import_irl_electricity": "Import", "import_nor_electricity": "Import"
}

def collect_and_scale_metrics(scenarios_dict, emissions_yaml):
    with open(emissions_yaml) as f:
        emissions_data = yaml.safe_load(f)["overrides"]["emissions"]["techs"]

    records = []

    for scen_name, filepath in scenarios_dict.items():
        model = calliope.read_netcdf(filepath)

        totals = {"RES":0, "LowC":0, "Fossil":0, "Stor":0, "Import":0}

        for tech, group in tech_group_map.items():
            try:
                flow = model.results.flow_out.loc[{"techs": tech}].sum().item()
            except KeyError:
                flow = 0
            totals[group] += flow

        co2 = 0
        for tech in tech_group_map:
            try:
                flow = model.results.flow_out.loc[{"techs": tech}].sum().item()
                ef = emissions_data.get(tech, {}).get("cost_flow_out", {}).get("data", 0)
                co2 += flow * ef
            except KeyError:
                continue

        cost = model.results.systemwide_levelised_cost.sum().item() 

        unmet_demand = model.results.systemwide_levelised_cost.sum().item() 

        if scen_name.endswith("_ND") or scen_name == "B1":
            pathway = "ND"
        elif scen_name.endswith("_FFR") or scen_name == "B2":
            pathway = "FFR"

        records.append({
            "scenario": scen_name,
            "pathway": pathway,
            **totals,
            "CO2": co2,
            "Cost": cost,
            "Unmet demand": unmet_demand,
        })

    df = pd.DataFrame(records)

    scaled_df_list = []
    metrics = ["RES", "LowC", "Stor", "Fossil", "Import", "CO2", "Cost"]
    for pathway in df["pathway"].unique():
        df_path = df[df["pathway"]==pathway].copy()
        max_per_metric = df_path[metrics].max().replace(0, 1) 
        df_path_scaled = df_path.copy()
        df_path_scaled[metrics] = df_path_scaled[metrics] / max_per_metric
        scaled_df_list.append(df_path_scaled)

    scaled_df = pd.concat(scaled_df_list, ignore_index=True)
    return scaled_df


# def plot_radar(metrics_df, output_paths):
#     metrics = ["RES", "LowC", "Stor", "Fossil", "Import"]
#     pathways = ["ND", "FFR"]

#     for pathway in pathways:
#         df = metrics_df[metrics_df.pathway == pathway].set_index("scenario")

#         labels = metrics
#         num_vars = len(labels)
#         angles = np.linspace(0, 2 * np.pi, num_vars, endpoint=False).tolist()
#         angles += angles[:1]  

#         fig, ax = plt.subplots(figsize=(8, 8), subplot_kw=dict(polar=True))
#         ax.set_theta_offset(np.pi / 2)  
#         ax.set_theta_direction(-1)
#         ax.set_rlabel_position(30)

#         cmap = plt.cm.get_cmap("tab10", len(df.index))

#         for i, scenario in enumerate(df.index):
#             values = df.loc[scenario, metrics].tolist()
#             values += values[:1]  
#             ax.plot(
#                 angles, values, label=scenario,
#                 linewidth=2, marker='o', markersize=5,
#                 color=cmap(i), alpha=0.8
#             )
#             ax.fill(angles, values, color=cmap(i), alpha=0.1)

#         ax.set_xticks(angles[:-1])
#         ax.set_xticklabels(labels, fontsize=12)
#         ax.set_yticks([0.25, 0.5, 0.75, 1.0])
#         ax.set_yticklabels(["0.25", "0.5", "0.75", "1.0"], fontsize=10)
#         ax.set_ylim(0, 1)
#         ax.grid(True)

#         ax.legend(loc="upper right", bbox_to_anchor=(1.3, 1.1), fontsize=10)
#         fig.tight_layout()
#         fig.savefig(output_paths[pathway], dpi=300)
#         plt.close(fig)


scenarios = {k: v for k, v in snakemake.input.items() if k != "emissions"}
metrics_df = collect_and_scale_metrics(scenarios, snakemake.input.emissions)
print(metrics_df)

#plot_radar(metrics_df, {"ND": snakemake.output.ND, "FFR": snakemake.output.FFR})


metrics = ["RES", "LowC", "Fossil", "Stor", "Import"]

for pathway in ["ND", "FFR"]:
    df = metrics_df[metrics_df.pathway == pathway].copy()
    
    if pathway == "ND":
        baseline = "B1"
    else:
        baseline = "B2"
    
    df = df[(df["scenario"] == baseline) | ((df["scenario"].str.startswith("S")) & (df["scenario"].str.endswith(pathway)))]
    
    df_plot = df.set_index("scenario")[metrics]
    
    plt.figure(figsize=(8,6))
    sns.heatmap(df_plot, annot=True, cmap="YlGnBu", cbar_kws={'label':'Scaled value (0-1)'}, annot_kws={"size": 14} )
    plt.title(f"{pathway} system metrics", fontsize=14)
    plt.tight_layout()
    plt.ylabel('')
    plt.yticks(fontsize=14)
    plt.xticks(fontsize=14)
    output_path = snakemake.output.ND if pathway == "ND" else snakemake.output.FFR
    plt.savefig(output_path, dpi=300)
    plt.close()