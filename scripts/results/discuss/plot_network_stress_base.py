import geopandas as gpd
import matplotlib.pyplot as plt
import yaml
from shapely.geometry import LineString
import numpy as np
import calliope
import pandas as pd
import matplotlib.cm as cm
import matplotlib.colors as mcolors


def load_yaml(path):
    with open(path) as f:
        return yaml.safe_load(f)

def plot_zones(ax, zones_gdf):
    zones_gdf.plot(ax=ax, color='lightsteelblue', edgecolor='black', legend=False)
    for idx, row in zones_gdf.iterrows():
        ax.annotate(
            row['z1'], 
            xy=(row.geometry.centroid.x, row.geometry.centroid.y + 0.2),
            ha='center', va='center', fontsize=14, fontweight='bold', color='black'
        )

def plot_centroids(ax, centroids_gdf, marker="*", color="black", size=100):
    centroids_gdf.plot(ax=ax, marker=marker, color=color, markersize=size)

def plot_onshore_links(ax, centroids_gdf, onshore_yaml_path, usage_df, cmap, norm):
    data = load_yaml(onshore_yaml_path)

    for link_name, link_info in data["techs"].items():
        start_zone = link_info["link_from"]
        end_zone = link_info["link_to"]

        try:
            start_point = centroids_gdf.loc[centroids_gdf['z1'] == start_zone].geometry.values[0]
            end_point = centroids_gdf.loc[centroids_gdf['z1'] == end_zone].geometry.values[0]
        except IndexError:
            continue

        if not link_info.get("flow_cap_max"):
            color = "purple"
        else:
            util_row = usage_df.loc[usage_df["tech"] == link_name]
            print(link_name, util_row["utilization_pct"].values)

            if not util_row.empty:
                util = float(util_row["utilization_pct"].values[0])
                color = cmap(norm(util))
            else:
                color = "grey"

        line = LineString([start_point, end_point])
        ax.plot(*line.xy, color=color, linewidth=2, alpha=0.9)


curvature_dict = {
    "z4_z7_offshore": (1.5, -1),
    "z4_z6_offshore": (1, 1),  
    "z2_z8_offshore": (0.7, 1),
    "z12_z17_offshore": (1.6, 1)
}

apex_offset_dict = {
    "z4_z6_offshore": (3, -0.1)  
}

def plot_offshore_links(ax, centroids_gdf, offshore_yaml_path, usage_df, cmap, norm,
            curvature_dict=None, apex_offset_dict=None, default_color="darkorange", linewidth=2):
    data = load_yaml(offshore_yaml_path)

    for link_name, link_info in data["techs"].items():
        start_zone = link_info["link_from"]
        end_zone = link_info["link_to"]

        try:
            start_point = centroids_gdf.loc[centroids_gdf['z1'] == start_zone].geometry.values[0]
            end_point = centroids_gdf.loc[centroids_gdf['z1'] == end_zone].geometry.values[0]
        except IndexError:
            continue

        factor = 0.5
        direction = 1
        if curvature_dict and link_name in curvature_dict:
            factor, direction = curvature_dict[link_name]

        mid_x = (start_point.x + end_point.x) / 2
        mid_y = (start_point.y + end_point.y) / 2

        dx = end_point.x - start_point.x
        dy = end_point.y - start_point.y
        perp_x = -dy * factor * direction
        perp_y = dx * factor * direction

        control = (mid_x + perp_x, mid_y + perp_y)

        if apex_offset_dict and link_name in apex_offset_dict:
            offset_x, offset_y = apex_offset_dict[link_name]
            control = (control[0] + offset_x, control[1] + offset_y)

        t = np.linspace(0, 1, 100)
        x = (1-t)**2 * start_point.x + 2*(1-t)*t * control[0] + t**2 * end_point.x
        y = (1-t)**2 * start_point.y + 2*(1-t)*t * control[1] + t**2 * end_point.y

        if not link_info.get("flow_cap_max"):  
            color = "purple"
        else:
            util_row = usage_df.loc[usage_df["tech"] == link_name]
            if not util_row.empty:
                util = util_row["utilization_pct"].values[0]
                color = cmap(norm(util))
            else:
                color = default_color

        ax.plot(x, y, color=color, linewidth=linewidth, alpha=0.9)



def collect_transmission_usage(filepath, techs_with_caps):
    model = calliope.read_netcdf(filepath)
    records = []
    timesteps = len(model.results.timesteps.values)

    for tech, cap in techs_with_caps.items():
        try:
            flow_in = model.results.flow_in.loc[{"techs": tech}].sum().item()
        except KeyError:
            flow_in = 0  

        try:
            flow_cap = model.results.link_flow_cap.loc[{"techs": tech}].item()
        except KeyError:
            flow_cap = 0.0

        capacity = flow_cap * timesteps

        utilization_pct = (flow_in / capacity) * 100 if capacity > 0 else 0

        #print(f"{tech} has %: {utilization_pct}")

        records.append({
            "tech": tech,
            "total_flow_GWh": flow_in,
            "capacity_GWh": capacity,
            "utilization_pct": utilization_pct
        })

    return pd.DataFrame.from_records(records)


def extract_tech_caps(yaml_paths):
    tech_caps = {}
    for path in yaml_paths:
        data = load_yaml(path)
        for tech, info in data["techs"].items():
            if "flow_cap_max" in info and info["flow_cap_max"] is not None:
                tech_caps[tech] = info["flow_cap_max"]
    return tech_caps



zones_gdf = gpd.read_file(snakemake.input.zones)
centroids_gdf = gpd.read_file(snakemake.input.centroids)

techs_with_caps = extract_tech_caps([snakemake.input.onshore_links, snakemake.input.offshore_links])
usage_nd = collect_transmission_usage(snakemake.input.ND, techs_with_caps)
usage_ffr = collect_transmission_usage(snakemake.input.FFR, techs_with_caps)

cmap = cm.RdYlGn_r
norm = mcolors.Normalize(vmin=0, vmax=100)

fig, axes = plt.subplots(1, 2, figsize=(14, 12), sharey=True)

for ax, usage_df, title in zip(axes, [usage_nd, usage_ffr], ["New Dispatch", "Further Flex and Renewables"]):
    plot_zones(ax, zones_gdf)
    plot_centroids(ax, centroids_gdf)
    plot_onshore_links(ax, centroids_gdf, snakemake.input.onshore_links, usage_df, cmap, norm)
    plot_offshore_links(ax, centroids_gdf, snakemake.input.offshore_links, usage_df, cmap, norm,
                        curvature_dict=curvature_dict, apex_offset_dict=apex_offset_dict)
    ax.set_title(title, fontsize=14)
    ax.axis("off")  

fig.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.15, wspace=0.02, hspace=0.02)

sm = cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar_ax = fig.add_axes([0.1, 0.05, 0.8, 0.03])  
cbar = fig.colorbar(sm, cax=cbar_ax, orientation="horizontal")
cbar.set_label("Transmission Utilization (%)", fontsize=11)

plt.savefig(snakemake.output.caps_comp, dpi=300, bbox_inches="tight")