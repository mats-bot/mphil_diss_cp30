import pandas as pd
import geopandas as gpd
import yaml

df = pd.read_csv(snakemake.input[0], index_col=0)
capacity_df = pd.read_excel(snakemake.input[1], sheet_name="Onshore - hydropower", header=0)
zones_gdf = gpd.read_file(snakemake.input[2])


hydro_df = capacity_df[
    (capacity_df["Country"] != "NI") &
    (capacity_df["Hydro Type"].str.lower() != "pumped storage")
].copy()

hydro_df["IC (MW)"] = pd.to_numeric(hydro_df["IC (MW)"], errors="coerce")
hydro_df = hydro_df.dropna(subset=["IC (MW)"])


hydro_df["Powerhouse Latitude"] = (
    hydro_df["Powerhouse Latitude"]
    .astype(str)
    .str.split("\n").str[0]
)
hydro_df["Powerhouse Longitude"] = (
    hydro_df["Powerhouse Longitude"]
    .astype(str)
    .str.split("\n").str[0]
)

hydro_df["Powerhouse Latitude"] = pd.to_numeric(hydro_df["Powerhouse Latitude"], errors="coerce")
hydro_df["Powerhouse Longitude"] = pd.to_numeric(hydro_df["Powerhouse Longitude"], errors="coerce")
hydro_df = hydro_df.dropna(subset=["Powerhouse Latitude", "Powerhouse Longitude"])

hydro_gdf = gpd.GeoDataFrame(
    hydro_df,
    geometry=gpd.points_from_xy(hydro_df["Powerhouse Longitude"], hydro_df["Powerhouse Latitude"]),
    crs="EPSG:4326"
)

joined_gdf = gpd.sjoin(hydro_gdf, zones_gdf[["z1", "geometry"]], how="left", predicate="within")

unassigned_mask = joined_gdf['z1'].isna()
unassigned_gdf = hydro_gdf[unassigned_mask].copy()

if not unassigned_gdf.empty:
    for idx, point in unassigned_gdf.geometry.items():
        nearest_zone = None
        min_distance = float('inf')
        for zone_idx, zone in zones_gdf.iterrows():
            distance = point.distance(zone.geometry)
            if distance < min_distance:
                min_distance = distance
                nearest_zone = zone['z1']
        joined_gdf.loc[idx, 'z1'] = nearest_zone

hydro_df['zone'] = joined_gdf['z1']
zone_capacity = hydro_df.groupby("zone")["IC (MW)"].sum()

output_df = pd.DataFrame({
    zone: [cap, cap] for zone, cap in zone_capacity.items()
}, index=["flow_cap_min", "flow_cap_max"])

output_df.to_csv(snakemake.output[1])

techs = ["Hydro"]

techs_yaml = {}

for tech in techs:
    tech_name = tech.lower()

    # flow_cap_data = tech_zone_capacity.get(tech_name, [0.0] * len(zones))
    # flow_cap_max = {"data": flow_cap_data, "dims": ["carriers"], "index": zones}
    techs_yaml[tech_name] = {
        "category": "renewable",
        "cp30_category": "renewable",
        "base_tech": "supply",
        "name": tech,
        "carrier_out": "electricity",
        "flow_out_eff": float(df.loc["efficiency", tech]),  # unitless (fraction)
        "lifetime": int(df.loc["lifetime", tech]),  # years
        "capacity_factor_max": 0.45,  # from Electricity Generation costs 2023 to avoid full load
        "cost_om_annual": {
            "data": float(df.loc["om_annual", tech]),
            "index": "monetary",
            "dims": ["costs"],
        },
        "cost_flow_out": {
            "data": float(df.loc["om_prod", tech]),
            "index": "monetary",
            "dims": ["costs"],
        },
        "cost_flow_in": {
            "data": float(df.loc["fuel_cost", tech]),
            "index": "monetary",
            "dims": ["costs"],
        },
        # Fix hydro capacity to today's values
        # "flow_cap_min": flow_cap_max,
        # "flow_cap_max": flow_cap_max,  
    }

data_tables_yaml = {
    "hydro_capacities": {
        "data": "data/processed/techs/hydro_capacities.csv",
        "rows": "parameters",
        "columns": "nodes",
        "add_dims": {
            "techs": tech_name
        }
    }
}

# Write output YAML
with open(snakemake.output[0], "w") as f:
    yaml.dump({"techs": techs_yaml, "data_tables": data_tables_yaml}, f, sort_keys=False)
