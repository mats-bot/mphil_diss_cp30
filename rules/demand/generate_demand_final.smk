rule prep_demand_data:
    input:
        "data/intermediates/demand/tzone_demand_aggregated.csv",
        "data/intermediates/demand/S2_tzone_demand_aggregated.csv",
    output:
        directory("data/intermediates/demand/demand_by_type"),
        directory("data/intermediates/demand/S2_demand_by_type")

    conda:
        "../../envs/data_processing.yaml"
    script:
        "../../scripts/demand/process_final_demand.py"

rule define_flexible_demand_share:
    input:
        "data/intermediates/demand/demand_by_type/",
        "data/intermediates/demand/{scenario}_flexibility_caps.csv",
        "data/intermediates/demand/S2_demand_by_type/"
    output:
        directory("data/processed/demand/split_by_flex/{scenario}"),
        directory("data/processed/demand/S2_split_by_flex/{scenario}")
    params:
        scenario=lambda wildcards: wildcards.scenario
    conda:
        "../../envs/data_processing.yaml"
    script:
        "../../scripts/demand/add_demand_flex.py"

rule generate_demand_yaml:
    input:
        "data/processed/demand/split_by_flex/{scenario}",
        "data/processed/demand/S2_split_by_flex/{scenario}"
    params:
        scenario=lambda wildcards: wildcards.scenario
    output:
        "demand/demand_{scenario}.yaml",
        "sensitivities/S2/demand_S2_{scenario}.yaml",
        "sensitivities/S5/demand_S5_{scenario}.yaml"

    conda:
        "../../envs/data_processing.yaml"
    script:
        "../../scripts/demand/generate_demand_yaml.py"

