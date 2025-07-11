rule generate_zones:
    input:
        "uploaded_data/tzones_centroids.gpkg"
    output:
        "spatial/zones.yaml"
    conda:
        "../../envs/gpkg_data.yaml"
    script:
        "../../scripts/spatial/create_zones.py"
    