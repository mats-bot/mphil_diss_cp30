configfile: "config/default.yaml"

include: "rules/clean_demand.smk"

rule all:
    input:
        config["demand_data"]["expanded"]


