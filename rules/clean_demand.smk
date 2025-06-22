rule download_demand_data:
    output:
        "data/raw/demand.xlsx"
    conda: 
        "../envs/data_processing.yaml"
    shell:
        """
        curl -L https://www.neso.energy/document/321056/download -o {output}
        """


rule clean_demand:
    input:
        "data/raw/demand.xlsx"
    output:
        "data/processed/cleaned_demand.csv"
    conda:
        "../envs/data_processing.yaml"
    script:
        "../scripts/clean_demand.py"