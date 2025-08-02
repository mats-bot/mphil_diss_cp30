rule create_offshore_links_yaml:
    input:
        "uploaded_data/offshore_links.xlsx"
    output:
        "spatial/offshore_transmission.yaml"
    conda:
        "../../envs/data_processing.yaml"
    script:
        "../../scripts/spatial/create_offshore_links.py"