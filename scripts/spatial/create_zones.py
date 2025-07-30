import geopandas as gpd
import yaml

centroids = gpd.read_file(snakemake.input[0])

locations_dict = {
    row["z1"]: {"coordinates": f"{float(row.geometry.y)}, {float(row.geometry.x)}"}
    for _, row in centroids.iterrows()
}

sorted_zones = sorted(locations_dict.items(), key=lambda x: int(x[0][1:]))

locations = {"nodes": dict(sorted_zones)}

with open(snakemake.output[0], "w") as f:
    yaml.dump(locations, f, sort_keys=False, default_flow_style=False, indent=4)
