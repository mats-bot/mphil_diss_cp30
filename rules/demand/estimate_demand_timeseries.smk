rule download_2013_demand_timeseries:
    output:
        "data/raw/demand/2013_demand_bihourly.csv",
        "data/raw/demand/2017_demand_bihourly.csv"
    conda:
        "../../envs/data_processing.yaml"
    params:
        year = config["demand_data"]["year"],
        S2_year = 2017
    script:
        "../../scripts/demand/download_2013_demand.py"


rule clean_2013_demand_timeseries:
    input:
        "data/raw/demand/2013_demand_bihourly.csv",
        "data/raw/demand/2017_demand_bihourly.csv"
    output:
        "data/intermediates/demand/2013_demand_cleaned.csv",
        "data/intermediates/demand/2017_demand_cleaned.csv",
    conda:
        "../../envs/data_processing.yaml"
    script:
        "../../scripts/demand/clean_2013_demand.py"


rule estimate_demand_timeseries:
    input:
        HISdemand = "data/intermediates/demand/2013_demand_cleaned.csv",
        FESdemand = "data/intermediates/demand/GSP_demand_coords.csv",
        S2_HISdemand = "data/intermediates/demand/2017_demand_cleaned.csv"
    output:
        "data/intermediates/demand/GSP_timeseries.csv",
        "data/intermediates/demand/S2_GSP_timeseries.csv",
    conda:
        "../../envs/data_processing.yaml"
    script:
        "../../scripts/demand/estimate_demand_hourly.py"
