rule download_GTCs:
    output: 
        "data/raw/FES_GTCs.xlsx"
    conda:
        "../../envs/data_processing.yaml"
    shell:
        """
        curl --ssl-no-revoke -sSL https://www.neso.energy/document/352111/download -o {output}
        """
