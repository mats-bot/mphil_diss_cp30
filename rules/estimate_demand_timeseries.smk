rule download_2013_demand_timeseries:
    output:
        "data/raw/2013_demand_bihourly.csv"
    conda:
        "../envs/data_processing.yaml"
    params:
        year = config["demand_data"]["year"]
    script:
        "../scripts/download_2013_demand.py"


rule clean_2013_demand_timeseries:
    input:
        "data/raw/2013_demand_bihourly.csv"
    output:
        "data/intermediates/2013_demand_cleaned.csv"
    conda:
        "../envs/data_processing.yaml"
    script:
        "../scripts/clean_2013_demand.py"


rule estimate_demand_timeseries:
    input:
        HISdemand = "data/intermediates/2013_demand_cleaned.csv",
        FESdemand = "data/intermediates/GSP_demand_coords.csv"
    output:
        "data/intermediates/GSP_timeseries.csv"
    conda:
        "../envs/data_processing.yaml"
    script:
        "../scripts/estimate_demand_hourly.py"
