rule generate_CFs:
    input:
        "uploaded_data/tzones.gpkg",
        "data/processed/techs/offshore_wind_projects.csv"

    output:
        "data/intermediates/spatial/atlite_cutout1_{year}.nc",
        "data/intermediates/spatial/atlite_cutout2_{year}.nc",
        "data/processed/spatial/solar_cf_{year}.csv",
        "data/processed/spatial/onshore_cf_{year}.csv",
        "data/processed/spatial/offshore_cfs_{year}.csv"
    conda:
        "../../envs/atlite_data.yaml"
    script:
        "../../scripts/techs/generate_atlite_cfs.py"