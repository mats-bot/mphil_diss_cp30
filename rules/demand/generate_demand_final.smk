rule prep_demand_data:
    input:
        "data/intermediates/demand/tzone_demand_aggregated.csv"
    output:
        directory("data/intermediates/demand/demand_by_type")
    conda:
        "../../envs/data_processing.yaml"
    script:
        "../../scripts/demand/process_final_demand.py"

rule define_flexible_demand_share:
    input:
        "data/intermediates/demand/demand_by_type/",
        "data/intermediates/demand/{scenario}_flexibility_caps.csv"
    output:
        directory("data/processed/demand/split_by_flex/{scenario}")
    params:
        scenario=lambda wildcards: wildcards.scenario
    conda:
        "../../envs/data_processing.yaml"
    script:
        "../../scripts/demand/add_demand_flex.py"

rule generate_demand_yaml:
    input:
        "data/processed/demand/split_by_flex/{scenario}"
    params:
        scenario=lambda wildcards: wildcards.scenario
    output:
        "demand_{scenario}.yaml"
    conda:
        "../../envs/data_processing.yaml"
    script:
        "../../scripts/demand/generate_demand_yaml.py"

