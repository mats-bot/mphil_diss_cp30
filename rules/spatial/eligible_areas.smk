rule calculate_eligible_areas:
    message: "Calculate eligible areas for each renewable technology type"
    input:
        tzones="uploaded_data/tzones.gpkg",
        eligible_land="data/raw/technically-eligible-land.tif",
        eligible_area="data/raw/technically-eligible-area-km2.tif"
    conda: "../../envs/geo.yaml"
    output: "data/interim/eligible_areas_km2.csv"
    script: "../../scripts/spatial/calculate_eligible_areas.py"
