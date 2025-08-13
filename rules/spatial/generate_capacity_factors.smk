MONTHS = [f"{m:02d}" for m in range(1,13)]
rule generate_onshore_cutouts:
    input:
        "uploaded_data/tzones.gpkg",
    output:
        protected("data/intermediates/spatial/onshore_cutouts/{weather_year}/cutout_{month}.nc"),
    params:
        features = ["onshore_wind", "pv"]
    conda:
        "../../envs/atlite_data.yaml"
    script:
        "../../scripts/spatial/generate_atlite_cutouts.py"

rule generate_offshore_region:
    message: "Create a geospatial region around planned offshore wind sites"
    input:
        "data/processed/techs/offshore_wind_projects.csv"
    output:
        "data/processed/spatial/offshore_wind_projects.gpkg"
    conda:
        "../../envs/gpkg_data.yaml"
    script:
        "../../scripts/spatial/offshore_regions.py"

rule generate_offshore_cutouts:
    input:
        rules.generate_offshore_region.output[0]
    output:
        protected("data/intermediates/spatial/offshore_cutouts/{weather_year}/cutout_{month}.nc"),
    params:
        features = ["offshore_wind"]
    conda:
        "../../envs/atlite_data.yaml"
    script:
        "../../scripts/spatial/generate_atlite_cutouts.py"


rule generate_onshore_CFs:
    input:
        zones = "uploaded_data/tzones.gpkg",
        cutouts = expand(
            "data/intermediates/spatial/onshore_cutouts/{{weather_year}}/cutout_{month}.nc",
            month=MONTHS
        )
    output:
        solar_cf = "data/intermediates/spatial/solar_cf_{weather_year}_raw.csv",
        onshore_cf = "data/intermediates/spatial/onshore_cf_{weather_year}_raw.csv",
    conda:
        "../../envs/atlite_data.yaml"
    script:
        "../../scripts/spatial/generate_atlite_cfs_onshore.py"

rule generate_offshore_CFs:
    input:
        zones = rules.generate_offshore_region.output[0],
        cutouts = expand(
            "data/intermediates/spatial/offshore_cutouts/{{weather_year}}/cutout_{month}.nc",
            month=MONTHS
        )
    output:
        offshore_cf = "data/intermediates/spatial/offshore_cf_{weather_year}_raw.csv"
    conda:
        "../../envs/atlite_data.yaml"
    script:
        "../../scripts/spatial/generate_atlite_cfs_offshore.py"


rule shift_cfs_year:
    input:
        "data/intermediates/spatial/{tech}_cf_{weather_year}_raw.csv",
    output:
        "data/processed/spatial/{tech}_cf_{weather_year}.csv",
    conda:
        "../../envs/data_processing.yaml"
    script:
        "../../scripts/spatial/shift_year.py"
