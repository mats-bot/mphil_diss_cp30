rule download_pypsa_techdata:
    output:
        "data/raw/techs/pypsa_2030_costs.csv"
    run:
        import os
        import requests

        url = "https://raw.githubusercontent.com/PyPSA/technology-data/cbf81d1429ded0d3eb18eb24a1d1bf3498b681cf/outputs/costs_2030.csv"
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)

        r = requests.get(url)
        r.raise_for_status()

        with open(output[0], "wb") as f:
            f.write(r.content)


rule download_DESNZ_costs:
    output:
        "data/raw/techs/DESNZ_costs.xlsx"
    conda:
        "../../envs/data_processing.yaml"
    shell:
        """
        curl --ssl-no-revoke -sSL https://assets.publishing.service.gov.uk/media/6555cb6d046ed4000d8b99bb/annex-a-additional-estimates-and-key-assumptions.xlsx -o {output}
        """

rule clean_DESNZ_costs:
    input:
        "data/raw/techs/DESNZ_costs.xlsx"
    output:
        "data/intermediates/techs/DESNZ_costs.csv"
    conda:
        "../../envs/data_processing.yaml"
    script:
        "../../scripts/techs/clean_DESNZ_costs.py"


rule add_missing_costs:
    input:
        "data/intermediates/techs/DESNZ_costs.csv"
    output:
        "data/processed/techs/storage_costs.csv",
        "data/processed/techs/generation_costs.csv"
    conda:
        "../../envs/data_processing.yaml"
    script:
        "../../scripts/techs/add_missing_costs.py"