rule generate_minimum_capacities:
    input: 
        "data/processed/techs/fossil_nuclear_2023.csv",
        "data/processed/techs/renewables_2023.csv",
    output: 
        "data/processed/spatial/min_zonal_caps.csv",
        "data/processed/spatial/max_zonal_caps.csv",
        "spatial/capacities_2023.yaml"
    conda: 
        "../../envs/data_processing.yaml"
    script:
        "../../scripts/techs/generate_2023_caps.py"

rule generate_renewables_max_zonal_capacities:
    input:
        "data/processed/techs/renewables_2030_queue.csv"
    output:
        "spatial/capacities_2030_REPD.yaml",
        "data/processed/techs/unconstrained_techs_2030.csv",
        "data/processed/techs/CHP_zero_caps.csv"
    conda: 
        "../../envs/data_processing.yaml"
    script:
        "../../scripts/techs/generate_2030_caps_REPD.py"

rule generate_GB_max_capacities_FES:
    input:
        "data/raw/demand/FES24_workbook.xlsx"
    output:
        "spatial/capacities_2030_FFR.yaml",
        "spatial/capacities_2030_ND.yaml",
    conda:
        "../../envs/data_processing.yaml"
    script:
        "../../scripts/techs/generate_2030_caps_FES.py"



