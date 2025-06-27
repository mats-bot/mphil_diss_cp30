configfile: "config/default.yaml"

include: "rules/clean_demand.smk"
include: "rules/estimate_demand_timeseries.smk" 

rule all:
    input:
        config["demand_data"]["expanded"]


