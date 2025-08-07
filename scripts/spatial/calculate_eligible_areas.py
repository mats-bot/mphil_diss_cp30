import geopandas as gpd
import rasterio
from rasterstats import zonal_stats

# See https://github.com/calliope-project/solar-and-wind-potentials/blob/main/lib/renewablepotentialslib/eligibility.py
tech_types = {
    "not_eligible": 0,
    "rooftop_pv": 250,
    "onshore_wind_and_pv": 180,  # PV and wind compete for this space
    "onshore_wind": 110,
    "offshore_wind": 40,
}
zones_gdf = gpd.read_file(snakemake.input.tzones)

# Each pixel in this map is assigned to a specific technology type that can occupy that pixel
with rasterio.open(snakemake.input.eligible_land) as src:
    land = src.read(1)

# Each pixel in this map gives the area in km² that can be occupied by the technology type that can occupy that pixel
with rasterio.open(snakemake.input.eligible_area) as src:
    area = src.read(1)
    meta = src.meta
    affine = src.transform

    zones_gdf_transformed = zones_gdf.to_crs(meta["crs"])

    for tech_type, land_val in tech_types.items():
        # Create a mask for the eligible land for each technology type
        mask = land == land_val
        eligible_area = area * mask

        # Calculate zonal statistics for the eligible area
        zonal_stats_result = zonal_stats(
            zones_gdf_transformed, eligible_area, affine=affine, stats="sum", nodata=0
        )
        zones_gdf_transformed[tech_type] = [
            result["sum"] for result in zonal_stats_result
        ]

zones_gdf_transformed.drop("geometry", axis=1).to_csv(snakemake.output[0], index=False)
