rule generate_CFs:
    input:
        "spatial"
    output:
        "data/processed/techs/storage_costs.csv", ####

    conda:
        "../../envs/atlite_data.yaml"
    script:
        "../../scripts/techs/generate_atlite_cfs.py"