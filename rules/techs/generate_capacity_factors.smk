rule generate_CFs:
    input:
        "uploaded_data/tzones.gpkg",
        "data/processed/techs/offshore_wind_projects.csv"
    output:
        onshore_cutouts = expand("data/intermediates/techs/onshore_cutouts/{weather_year}/cutout_{month}.nc", weather_year="{weather_year}", month=[f"{m:02d}" for m in range(1,13)]),
        solar_cf = "data/processed/spatial/solar_cf_{weather_year}.csv",
        onshore_cf = "data/processed/spatial/onshore_cf_{weather_year}.csv",
        offshore_cutouts = expand("data/intermediates/techs/offshore_cutouts/{weather_year}/offshore_cutout_{month}.nc", weather_year="{weather_year}", month=[f"{m:02d}" for m in range(1,13)]),
        offshore_cf = "data/processed/spatial/offshore_cfs_{weather_year}.csv"
    conda:
        "../../envs/atlite_data.yaml"
    script:
        "../../scripts/techs/generate_atlite_cfs.py"