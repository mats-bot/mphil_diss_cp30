rule compare_connection_queues:
    input:
        "data/processed/spatial/renewables_2030_queue.csv",
        "data/raw/techs/CP30_workbook.xlsx"
    output:
        "results/figures/methods/total_2030_capacities.png"
    conda:
        "../../envs/plotting.yaml"
    script:
        "../../scripts/results/methods/plot_2030_queue.py"
