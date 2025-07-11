rule map_GSPs_to_tzones:
    input:
        "data/intermediates/GSP_timeseries.csv",
        "uploaded_data/tzones.gpkg"
    output:
        "data/intermediates/tzone_demand.csv"
    conda:
        "../envs/gpkg_data.yaml"
    script:
        "../scripts/demand_to_tzones.py"


rule aggregate_demand_tzones:
    input: 
        "data/intermediates/tzone_demand.csv"
    output:
        "data/intermediates/tzone_demand_aggregated.csv"
    conda:
        "../envs/data_processing.yaml"
    script:
        "../scripts/aggregate_tzone_demand.py"