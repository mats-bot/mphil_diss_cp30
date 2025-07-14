rule download_GTCs:
    output: 
        "data/raw/FES_GTCs.xlsx"
    conda:
        "../../envs/data_processing.yaml"
    shell:
        """
        curl --ssl-no-revoke -sSL https://www.neso.energy/document/352111/download -o {output}
        """

rule clean_onshore_transmission_caps:
    input: 
        "data/raw/FES_GTCs.xlsx",
        "uploaded_data/GTC_boundaries.xlsx"
    output:
        "data/intermediates/spatial/GTCs_linked.csv"
    conda:
        "../../envs/data_processing.yaml"
    script:
        "../../scripts/spatial/clean_onshore_links.py"


rule calculate_onshore_transmission_losses:
    input:
        "data/intermediates/spatial/GTCs_linked.csv",
        "uploaded_data/tzones_centroids.gpkg"
    output:
        "data/processed/spatial/onshore_transmission.csv"
    conda: 
        "../../envs/gpkg_data.yaml"
    script:
        "../../scripts/spatial/calculate_onshore_losses.py"


rule generate_onshore_transmission_yaml:
    input:
        "data/processed/spatial/onshore_transmission.csv"
    output: 
        "spatial/onshore_transmission.yaml"
    conda:
        "../../envs/data_processing.yaml"
    script:
        "../../scripts/spatial/create_onshore_links.py"

