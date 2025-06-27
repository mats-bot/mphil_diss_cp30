rule download_demand_data:
    output:
        "data/raw/demand.xlsx"
    conda: 
        "../envs/data_processing.yaml"
    shell:
        """
        curl --ssl-no-revoke -sSL https://www.neso.energy/document/321056/download -o {output}
        """
        # curl --ssl-no-revoke -L https://www.neso.energy/document/321056/download -o {output}


rule extract_demand_GSPs:
    input:
        "data/raw/demand.xlsx"
    output:
        "data/intermediates/GSP_demand_peaks.csv"
    conda:
        "../envs/data_processing.yaml"
    script:
        "../scripts/extract_GSP_demand_peaks.py"
        


rule extract_GSP_coordinates:
    input:
        "data/raw/demand.xlsx"
    output:
        "data/intermediates/GSP_coords.csv"
    conda:
        "../envs/data_processing.yaml"
    script:
        "../scripts/extract_GSP_coords.py"



rule combine_demand_files:
    input:
        demand = "data/intermediates/GSP_demand_peaks.csv",
        coords = "data/intermediates/GSP_coords.csv"

    output:
        "data/intermediates/GSP_demand_coords.csv"
    conda:
        "../envs/data_processing.yaml"
    script:
        "../scripts/combine_demand_coords.py"


