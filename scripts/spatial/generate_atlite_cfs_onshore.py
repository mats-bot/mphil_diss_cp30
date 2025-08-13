import atlite
import geopandas as gpd
import pandas as pd

zones = (
    gpd.read_file(snakemake.input.zones)
    .set_index("z1")
    .sort_index(key=lambda x: x.str.removeprefix("z").astype(int))
)


solar_dfs = []
wind_dfs = []

for cutout_path in snakemake.input.cutouts:
    cutout = atlite.Cutout(path=cutout_path)

    solar_cf = cutout.pv(
        panel="CSi", orientation="latitude_optimal", shapes=zones, per_unit=True
    ).to_pandas()

    solar_dfs.append(solar_cf)

    wind_cf = cutout.wind(
        turbine="NREL_ReferenceTurbine_2020ATB_7MW",
        shapes=zones,
        add_cutout_windspeed=True,
        per_unit=True,
    ).to_pandas()

    wind_dfs.append(wind_cf)


pd.concat(solar_dfs).to_csv(snakemake.output["solar_cf"])
pd.concat(wind_dfs).to_csv(snakemake.output["onshore_cf"])
