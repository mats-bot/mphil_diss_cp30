rule download_demand_data:
    output:
        "data/raw/demand/demand.xlsx"
    conda:
        "../../envs/data_processing.yaml"
    shell:
        """
        curl --ssl-no-revoke -sSL https://www.neso.energy/document/321056/download -o {output}
        """


rule extract_demand_GSPs:
    input:
        "data/raw/demand/demand.xlsx"
    output:
        "data/intermediates/demand/GSP_demand_peaks.csv",
    conda:
        "../../envs/data_processing.yaml"
    script:
        "../../scripts/demand/extract_GSP_demand_peaks.py"



rule extract_GSP_coordinates:
    input:
        "data/raw/demand/demand.xlsx"
    output:
        "data/intermediates/demand/GSP_coords.csv"
    conda:
        "../../envs/data_processing.yaml"
    script:
        "../../scripts/demand/extract_GSP_coords.py"



rule combine_demand_files:
    input:
        demand = "data/intermediates/demand/GSP_demand_peaks.csv",
        coords = "data/intermediates/demand/GSP_coords.csv",
        unmatched_GSP_coords = "uploaded_data/GSP_unmatched.xlsx"
    output:
        "data/intermediates/demand/GSP_demand_coords.csv"
    conda:
        "../../envs/data_processing.yaml"
    script:
        "../../scripts/demand/combine_demand_coords.py"


