rule prep_demand_data:
    input:
        "data/intermediates/demand/tzone_demand_aggregated.csv"
    output:
        directory("data/processed/demand/split_by_zone_and_type")
    conda:
        "../../envs/data_processing.yaml"
    script:
        "../../scripts/demand/process_final_demand.py"

rule generate_demand_yaml:
    input:
        "data/processed/demand/split_by_zone_and_type/"
    output:
        "demand.yaml"
    conda:
        "../../envs/data_processing.yaml"
    script:
        "../../scripts/demand/generate_demand_yaml.py"
