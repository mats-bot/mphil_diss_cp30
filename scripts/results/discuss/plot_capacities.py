import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import calliope

inputs = snakemake.input

scenario_labels = {
    "CP30": "CP30",
    "ND": "Base ND",
    "FFR": "Base FFR",
    "S1_ND": "S1 ND",
    "S1_FFR": "S1 FFR",
    "S2_ND": "S2 ND",
    "S2_FFR": "S2 FFR",
    "S3_ND": "S3 ND",
    "S3_FFR": "S3 FFR",
    "S4_ND": "S4 ND",
    "S4_FFR": "S4 FFR",
    "S5_ND": "S5 ND",
    "S5_FFR": "S5 FFR",
}
tech_colors = {
    "onshore_wind": "#33a02c",
    "offshore_wind": "#1f78b4",  
    "solar_pv": "#ffb703",       
}
techs_of_interest = list(tech_colors.keys())


def collect_excel_data(filepath):
    df = pd.read_excel(filepath, sheet_name="CP.04", usecols="M:Q", skiprows=5, nrows=20)
    df = df.rename(columns={
        "Technology": "tech_label",
        "By 2030": "capacity_GW"
    })
    df = df.dropna(subset=["tech_label", "capacity_GW"])

    tech_map = {
    "Solar - Further Flex and Renewables": "solar_pv",
    "Solar - New Dispatch": "solar_pv",
    "Onshore wind - Further Flex and Renewables": "onshore_wind",
    "Onshore wind - New Dispatch": "onshore_wind",
    "Offshore wind - Further Flex and Renewables": "offshore_wind",
    "Offshore wind - New Dispatch": "offshore_wind",
    }

    records = []
    for _, row in df.iterrows():
        if row["tech_label"] not in tech_map:
            continue
        tech = tech_map[row["tech_label"]]
        pathway = "ND" if "New Dispatch" in row["tech_label"] else "FFR"

        records.append({
            "scenario": "CP30",  
            "pathway": pathway,
            "tech": tech,
            "capacity_GW": row["capacity_GW"]
        })

    return pd.DataFrame(records)


def collect_capacity_data(inputs, scenario_labels, techs_of_interest):
    records = []
    for scen in scenario_labels:
        filepath = inputs.get(scen)
        if filepath is None:
            print(f"No file for scenario {scen}, skipping.")
            continue

        print(f"Loading {filepath} for {scen}...")
        model = calliope.read_netcdf(filepath)

        for tech in techs_of_interest:
            try:
                df = model.results.flow_cap.loc[{"techs": tech}]
                cap_val = df.sum().item() 
            except KeyError:
                cap_val = 0.0
            print(f"{scen} | {tech} | {cap_val}")
            records.append({
                "scenario": scen,
                "pathway": "ND" if "ND" in scen else "FFR",
                "tech": tech,
                "capacity_GW": cap_val
            })

    return pd.DataFrame.from_records(records)


def plot_pathway(df, pathway, output_path):
    subset = df[df["pathway"] == pathway]

    plt.figure(figsize=(8, 6))
    scenario_order = [
    "CP30", "ND", "S1_ND", "S2_ND", "S3_ND", "S4_ND", "S5_ND",
    "FFR", "S1_FFR", "S2_FFR", "S3_FFR", "S4_FFR", "S5_FFR"
    ]

    sns.barplot(
        data=subset,
        x="scenario",
        y="capacity_GW",
        hue="tech",
        palette=tech_colors,
        order=[s for s in scenario_order if s in subset["scenario"].unique()]
    )

    plt.xticks(rotation=45, ha="right")
    plt.ylabel("Installed capacity [GW]")
    plt.xlabel("Scenario")
    plt.title(f"Installed capacity comparison – {pathway} pathway")
    plt.legend(title="Technology")

    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()

df = collect_capacity_data(inputs, scenario_labels, techs_of_interest)
excel_path = snakemake.input.cp30  
df_excel = collect_excel_data(excel_path)

df_all = pd.concat([df, df_excel], ignore_index=True)

for pathway, out_path in [("ND", snakemake.output.caps_comp_nd),
                          ("FFR", snakemake.output.caps_comp_ffr)]:
    plot_pathway(df_all, pathway, out_path)