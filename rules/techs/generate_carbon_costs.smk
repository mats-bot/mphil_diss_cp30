rule download_pypsa_emissions:
    output:
        "data/raw/techs/pypsa_GB_emissions.csv"
    conda:
        "../../envs/data_processing.yaml"
    shell:
        """
        curl --ssl-no-revoke -sSL https://github.com/andrewlyden/PyPSA-GB/raw/refs/heads/master/data/emissions_intensity_by_types.csv -o {output}
        """

rule generate_emissions_files:
    input: 
        "data/raw/techs/pypsa_GB_emissions.csv"
    output:
        "sensitivities/S3/emissions.yaml"
    conda:
        "../../envs/data_processing.yaml"
    script:
        "../../scripts/techs/process_emissions.py"
