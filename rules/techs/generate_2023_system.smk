rule download_2024_operational:
    output: 
        "data/raw/techs/DUKES_power_2024.xlsx"
    conda:
        "../../envs/data_processing.yaml"
    shell:
        """
        curl --ssl-no-revoke -sSL https://assets.publishing.service.gov.uk/media/66a7daa4fc8e12ac3edb068e/DUKES_5.11.xlsx -o {output}
        """


rule download_CP30_technologies:
    output:
        "data/raw/techs/CP30_workbook.xlsx"
    conda:
        "../../envs/data_processing.yaml"
    shell:
        """
        curl --ssl-no-revoke -sSL https://www.neso.energy/document/346781/download -o {output}
        """

rule download_REPD_queue:
    output:
        "data/raw/techs/REPD_queue.xlsx"
    conda:
        "../../envs/data_processing.yaml"
    shell:
        """
        curl --ssl-no-revoke -sSL https://assets.publishing.service.gov.uk/media/6841899cd5f6bf75703e878b/repd-q1-apr-2025.xlsx -o {output}
        """

rule clean_REPD_queue:
    input:
        "data/raw/techs/REPD_queue.xlsx",
        "uploaded_data/GSP_unmatched.xlsx",
        "uploaded_data/tzones.gpkg"
    output:
        "data/intermediates/techs/offshore_wind_queue.csv",
        "data/intermediates/techs/renewables_2023.csv"
    conda:
        "../../envs/gpkg_data.yaml"
    script:
        "../../scripts/techs/clean_REPD.py"

rule generate_CP30_techs:
    input:
        "data/raw/techs/CP30_workbook.xlsx"
    output:
        "data/intermediates/techs/CP30_techs.csv"
    conda:
        "../../envs/data_processing.yaml"
    script:
        "../../scripts/techs/generate_CP30_techs.py"


rule clean_DUKES: 
    input:
        "data/raw/techs/DUKES_power_2024.xlsx",
        "uploaded_data/GSP_unmatched.xlsx",
        "uploaded_data/tzones.gpkg"
    output: 
        "data/intermediates/techs/fossil_nuclear_2023.csv"
    conda:
        "../../envs/gpkg_data.yaml"
    script:
        "../../scripts/techs/clean_2024_DUKES.py"


rule clean_offshore_wind:
    input: 
        "data/intermediates/techs/offshore_wind_queue.csv",
        "uploaded_data/tzones.gpkg"
    output:
        "data/processed/techs/offshore_wind_projects.csv"
    conda:
        "../../envs/gpkg_data.yaml"
    script:
        "../../scripts/techs/clean_offshore_wind.py"



