rule compare_connection_queues:
    input:
        "data/processed/techs/renewables_2030_queue.csv",
        "data/raw/techs/CP30_workbook.xlsx",
        "data/processed/techs/offshore_wind_projects_aggregated.csv"
    output:
        "results/figures/methods/total_2030_capacities.png"
    conda:
        "../../envs/plotting.yaml"
    script:
        "../../scripts/results/methods/plot_2030_queue.py"

rule plot_zones_links:
    input:
        zones = "uploaded_data/tzones.gpkg",
        centroids = "uploaded_data/tzones_centroids.gpkg",
        onshore_links = "spatial/onshore_transmission.yaml",
        offshore_links = "spatial/offshore_transmission.yaml",
        interconnectors = "techs/import_export.yaml"
    output:
        "results/figures/methods/zones_links.png"
    conda:
        "../../envs/plotting.yaml"
    script:
        "../../scripts/results/methods/plot_zones_links.py"
