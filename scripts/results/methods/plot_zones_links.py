import geopandas as gpd
import matplotlib.pyplot as plt
import yaml
from shapely.geometry import LineString
from shapely.geometry import Point
import numpy as np
import matplotlib.lines as mlines


country_coords = {
    "bel": (7.0, 51.0),
    "deu": (10.0, 54.0),
    "dnk": (11.5, 56.0),
    "fra": (2.0, 47.0),
    "irl": (-8.0, 53.0),
    "nld": (5.0, 52.0),
    "nor": (10.0, 60.0)
}

country_names = {
    "bel": "Belgium",
    "deu": "Germany",
    "dnk": "Denmark",
    "fra": "France",
    "irl": "Ireland",
    "nld": "Netherlands",
    "nor": "Norway"
}

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

def plot_onshore_links(ax, centroids_gdf, onshore_yaml_path):
    data = load_yaml(onshore_yaml_path)
    for link_name, link_info in data["techs"].items():
        start_zone = link_info["link_from"]
        end_zone = link_info["link_to"]
        
        try:
            start_point = centroids_gdf.loc[centroids_gdf['z1'] == start_zone].geometry.values[0]
            end_point = centroids_gdf.loc[centroids_gdf['z1'] == end_zone].geometry.values[0]
        except IndexError:
            print(f"Warning: centroid not found for {start_zone} or {end_zone}")
            continue
        
        color = 'red' if link_info.get("flow_cap_max") else 'purple'
        
        line = LineString([start_point, end_point])
        ax.plot(*line.xy, color=color, linewidth=2, alpha=0.7)


curvature_dict = {
    "z4_z7_offshore": (1.5, -1),
    "z4_z6_offshore": (1, 1),  
    "z2_z8_offshore": (0.7, 1),
    "z12_z17_offshore": (1.6, 1)
}

apex_offset_dict = {
    "z4_z6_offshore": (3, -0.1)  
}

def plot_offshore_links(ax, centroids_gdf, offshore_yaml_path, curvature_dict=None, apex_offset_dict=None, color="darkorange", linewidth=2):
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

        ax.plot(x, y, color=color, linewidth=linewidth, alpha=0.8)


def plot_interconnectors(ax, centroids_gdf, yaml_file, country_coords, country_names, line_length=2.0):
    with open(yaml_file, 'r') as f:
        data = yaml.safe_load(f)

    for node_id, node_data in data.get("nodes", {}).items():
        try:
            node_point = centroids_gdf.loc[centroids_gdf['z1'] == node_id].geometry.values[0]
        except IndexError:
            continue  
        for tech_name, tech_info in node_data.get("techs", {}).items():
            if not tech_name.startswith("import_"):
                continue

            parts = tech_name.split("_")
            if len(parts) < 3:
                continue
            _, country_code, _ = parts  

            if country_code not in country_coords:
                continue

            country_lon, country_lat = country_coords[country_code]


            dx = country_lon - node_point.x
            dy = country_lat - node_point.y
            length = (dx**2 + dy**2)**0.5
            if length == 0:
                continue

            dx_scaled = dx / length * line_length
            dy_scaled = dy / length * line_length

            ax.arrow(
                node_point.x,
                node_point.y,
                dx_scaled,
                dy_scaled,
                color='green',
                linewidth=1,
                alpha=0.8,
                head_width=0.05, 
                head_length=0.1,    
                length_includes_head=True
            )

            

zones_gdf = gpd.read_file(snakemake.input.zones)
centroids_gdf = gpd.read_file(snakemake.input.centroids)



fig, ax = plt.subplots(figsize=(7, 12))
plot_zones(ax, zones_gdf)
plot_centroids(ax, centroids_gdf)
plot_onshore_links(ax, centroids_gdf, snakemake.input.onshore_links)
plot_offshore_links(ax, centroids_gdf, snakemake.input.offshore_links, curvature_dict=curvature_dict, apex_offset_dict=apex_offset_dict)
plot_interconnectors(ax, centroids_gdf, snakemake.input.interconnectors, country_coords=country_coords, country_names=country_names)


all_geoms = zones_gdf.geometry.unary_union.union(centroids_gdf.geometry.unary_union)
for link_file in [snakemake.input.onshore_links, snakemake.input.offshore_links, snakemake.input.interconnectors]:
    data = load_yaml(link_file)
    for node_id, node_data in data.get("nodes", {}).items():
        for tech_info in node_data.get("techs", {}).values():
            if "link_from" in tech_info and "link_to" in tech_info:
                try:
                    start_pt = centroids_gdf.loc[centroids_gdf['z1'] == tech_info["link_from"]].geometry.values[0]
                    end_pt = centroids_gdf.loc[centroids_gdf['z1'] == tech_info["link_to"]].geometry.values[0]
                    all_geoms = all_geoms.union(LineString([start_pt, end_pt]))
                except IndexError:
                    continue

minx, miny, maxx, maxy = all_geoms.bounds
x_margin_left = (maxx - minx) * 0
x_margin_right = (maxx - minx) * 0
y_margin_top = (maxy - miny) * 0
y_margin_bottom = (maxy - miny) * 0

ax.set_xlim(minx + x_margin_left, maxx - x_margin_right)
ax.set_ylim(miny - y_margin_bottom, maxy + y_margin_top)

star_marker = mlines.Line2D([], [], color='black', marker='*', linestyle='None',
                          markersize=10, label='Zone centroid')
onshore_line = mlines.Line2D([], [], color='red', linewidth=2, label='Onshore link')
onshore_line_uncapped = mlines.Line2D([], [], color='purple', linewidth=2, label='Uncapped onshore link')
offshore_line = mlines.Line2D([], [], color='darkorange', linewidth=2, label='Offshore link')
inter_line = mlines.Line2D([], [], color='green', linewidth=2, label='Interconnector')

ax.legend(handles=[star_marker, onshore_line, onshore_line_uncapped, offshore_line, inter_line],
          loc='upper left', fontsize=13)


ax.set_axis_off()
plt.tight_layout()
plt.savefig(snakemake.output[0], dpi=300)
