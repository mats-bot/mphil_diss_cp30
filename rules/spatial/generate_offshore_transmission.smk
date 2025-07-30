rule calculate_onshore_transmission_losses:
    input:
        "uploaded_data/offshore_links.xlsx"
    output:
        "spatial/offshore_transmission.yaml"
    conda:
        "../../envs/data_processing.yaml"
    script:
        "../../scripts/spatial/create_offshore_links.py"