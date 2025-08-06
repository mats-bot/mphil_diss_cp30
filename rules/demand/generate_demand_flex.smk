rule download_demand_flex_data:
    output:
        "data/raw/demand/FES24_workbook.xlsx"
    conda:
        "../../envs/data_processing.yaml"
    shell:
        """
        curl --ssl-no-revoke -sSL https://www.neso.energy/document/321051/download -o {output}
        """

rule normalize_flex_data:
    input:
        "data/raw/demand/FES24_workbook.xlsx"
    output:
        "data/intermediates/demand/ND_flexibility_caps.csv",
        "data/intermediates/demand/FFR_flexibility_caps.csv",
    conda:
        "../../envs/data_processing.yaml"
    script:
        "../../scripts/demand/normalize_flex.py"