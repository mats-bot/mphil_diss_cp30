rule map_GSPs_to_tzones:
    input:
        "data/intermediates/GSP_timeseries.csv",
        "tzones.gpkg"
    output:
        "data/intermediates/tzone_demand.csv"
    conda:
        "../envs/data_processing.yaml"
    script:
        "../scripts/demand_to_tzones.py"


