rule map_GSPs_to_tzones:
    input:
        "data/intermediates/demand/GSP_timeseries.csv",
        "uploaded_data/tzones.gpkg",
        "data/intermediates/demand/S2_GSP_timeseries.csv"
    output:
        "data/intermediates/demand/tzone_demand.csv",
        "data/intermediates/demand/S2_tzone_demand.csv"

    conda:
        "../../envs/gpkg_data.yaml"
    script:
        "../../scripts/demand/demand_to_tzones.py"


rule aggregate_demand_tzones:
    input: 
        "data/intermediates/demand/tzone_demand.csv",
        "data/intermediates/demand/S2_tzone_demand.csv",
    output:
        "data/intermediates/demand/tzone_demand_aggregated.csv",
        "data/intermediates/demand/S2_tzone_demand_aggregated.csv",
    conda:
        "../../envs/data_processing.yaml"
    script:
        "../../scripts/demand/aggregate_tzone_demand.py"