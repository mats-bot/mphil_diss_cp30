import atlite
import geopandas as gpd
import pandas as pd

zones = gpd.read_file(snakemake.input.zones).set_index("Site Name")


offshore_dfs = []

for cutout_path in snakemake.input.cutouts:
    cutout = atlite.Cutout(path=cutout_path)

    cf_offshore = cutout.wind(
        turbine="NREL_ReferenceTurbine_2020ATB_15MW_offshore",
        shapes=zones,
        add_cutout_windspeed=True,
        per_unit=True,
    )

    df_month = cf_offshore.to_pandas()
    offshore_dfs.append(df_month)

aggregate_df = pd.concat(offshore_dfs)

aggregate_df.to_csv(snakemake.output["offshore_cf"])
