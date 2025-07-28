rule generate_CFs:
    input:
        "uploaded_data/tzones.gpkg",
        "data/processed/techs/offshore_wind_projects.csv"
    output:
        # "data/intermediates/spatial/atlite_cutout1a_{year}.nc",
        # "data/intermediates/spatial/atlite_cutout1b_{year}.nc",
        # "data/intermediates/spatial/atlite_cutout1c_{year}.nc",
        # "data/intermediates/spatial/atlite_cutout1d_{year}.nc",
        protected(directory("data/intermediates/techs/onshore_cutouts/{year}")),
        "data/processed/spatial/solar_cf_{year}.csv",
        "data/processed/spatial/onshore_cf_{year}.csv",
        protected(directory("data/intermediates/techs/offshore_cutouts/{year}")),
        "data/processed/spatial/offshore_cfs_{year}.csv"
    conda:
        "../../envs/atlite_data.yaml"
    script:
        "../../scripts/techs/generate_atlite_cfs.py"