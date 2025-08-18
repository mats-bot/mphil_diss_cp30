import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import calliope

inputs = snakemake.input
output = snakemake.output

# Define variables
flow_groups = {
    "imports": {"category": "import"},
    "LDES": {"category": "LDES"},
    "batteries": {"category": "Battery"},
    "fossil fuels": {"techs": ["gas_ccgt_existing", "gas_ccgt_new", "gas_ocgt_existing", "diesel_existing"]},
    "low carbon": {"techs": ["beccs", "ccs", "biomass", "hydrogen", "waste", "hydro"]},
}

tech_order = ["imports", "LDES", "batteries", "fossil fuels", "low carbon"]
tech_colors = {
    "imports": "#1f77b4",
    "LDES": "#ff7f0e",
    "batteries": "#2ca02c",
    "fossil fuels": "#d62728",
    "low carbon": "#9467bd",
}

scenario_labels = {
    "B1": "Base",
    "B2": "Base",
    "S1_ND": "S1",
    "S1_FFR": "S1",
    "S2_ND": "S2",
    "S2_FFR": "S2",
    "S3_ND": "S3",
    "S3_FFR": "S3",
    "S4_ND": "S4",
    "S4_FFR": "S4",
    "S5_ND": "S5",
    "S5_FFR": "S5",
}

def collect_model_flows(inputs, pathway):
    records = []
    for scen, label in scenario_labels.items():
        # Only keep scenarios for this pathway
        if (pathway not in scen) and not (
            (pathway == "ND" and scen == "B1") or (pathway == "FFR" and scen == "B2")
        ):
            continue

        filepath = inputs[scen]
        model = calliope.read_netcdf(filepath)

        if not hasattr(model.results, "flow_out"):
            continue  # Skip models with no flow_out

        techs_da = model.inputs.techs  # xarray DataArray of techs

        for name, cond in flow_groups.items():
            if "category" in cond:
                # select techs matching the category
                tech_list = techs_da.sel(techs=techs_da["category"] == cond["category"]).techs.values.tolist()
            else:  # specific tech list
                tech_list = [t for t in techs_da.techs.values if t in cond["techs"]]

            if not tech_list:
                continue  # skip if no matching techs

            flows = model.results.flow_out.sel(techs=tech_list)
            value = flows.sum(dim=["timesteps", "nodes", "techs"]).item()
            records.append({"scenario": label, "flow_type": name, "GWh": value})

    return pd.DataFrame(records)


def collect_cp30_flows(filepath, pathway):
    df = pd.read_excel(filepath, sheet_name="ES1", header=9)  # row 10 is header
    df = df[["Pathway", "Variable", "SubType", "Type", 2030]]  # 2030 column as int
    df = df.dropna(subset=[2030])

    if pathway == "ND":
        df_path = df[df["Pathway"] == "New Dispatch"]
    else:
        df_path = df[df["Pathway"] == "Further Flex and Renewables"]

    records = []
    for name in tech_order:
        if name == "imports":
            value = df_path[df_path["Variable"] == "Import (GWh)"][2030].sum()
        elif name == "LDES":
            value = df_path[(df_path["SubType"] == "LDES") & (df_path["Variable"] == "Generation (GWh)")][2030].sum()
        elif name == "batteries":
            value = df_path[(df_path["SubType"] == "Battery") & (df_path["Variable"] == "Generation (GWh)")][2030].sum()
        elif name == "fossil fuels":
            value = df_path[(df_path["Type"].isin(["Other Thermal", "Gas"])) & (df_path["Variable"] == "Generation (GWh)")][2030].sum()
        elif name == "low carbon":
            value = df_path[(df_path["Type"].isin(["Biomass", "CCS", "Hydro", "Low Carbon", "Waste"])) & (df_path["Variable"] == "Generation (GWh)")][2030].sum()
        records.append({"scenario": "CP30", "flow_type": name, "GWh": value})

    return pd.DataFrame(records)

def plot_flows(df, pathway, output_path):
    subset = df[df["scenario"].isin(["CP30"] + [s for s in scenario_labels.values() if s != "Base"])]
    plt.figure(figsize=(10, 6))
    sns.barplot(
        data=subset,
        x="scenario",
        y="GWh",
        hue="flow_type",
        hue_order=tech_order,
        palette=tech_colors
    )
    plt.ylabel("Flow [GWh]")
    plt.xlabel("Scenario")
    plt.title(f"Flexibility flows – {pathway} pathway")
    plt.xticks(rotation=0)
    plt.legend(title="Flow type")

    # gridlines and ticks
    plt.grid(axis="y", which="major", linestyle="--", color="gray", alpha=0.5)
    plt.gca().yaxis.set_minor_locator(plt.MultipleLocator(10))
    plt.grid(axis="y", which="minor", linestyle=":", color="gray", alpha=0.3)

    # scale font
    plt.rcParams.update({'font.size': 11})

    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()

# Collect data
df_nd_model = collect_model_flows(inputs, "ND")
df_ffr_model = collect_model_flows(inputs, "FFR")

df_nd_cp30 = collect_cp30_flows(inputs.cp30, "ND")
df_ffr_cp30 = collect_cp30_flows(inputs.cp30, "FFR")

df_nd_all = pd.concat([df_nd_cp30, df_nd_model], ignore_index=True)
df_ffr_all = pd.concat([df_ffr_cp30, df_ffr_model], ignore_index=True)

# Plot
plot_flows(df_nd_all, "ND", output.flows_flex_nd)
plot_flows(df_ffr_all, "FFR", output.caps_flex_ffr)
