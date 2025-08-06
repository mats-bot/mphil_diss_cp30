configfile: "config/default.yaml"

# Generate data files and definitions for model
include: "rules/demand/clean_demand.smk"
include: "rules/demand/estimate_demand_timeseries.smk"
include: "rules/demand/aggregate_demand_tzone.smk"
include: "rules/demand/generate_demand_final.smk"
include: "rules/spatial/create_zones.smk"
include: "rules/spatial/generate_onshore_transmission.smk"
include: "rules/techs/generate_2023_system.smk"
include: "rules/techs/generate_monetary_costs.smk"
include: "rules/techs/generate_capacity_factors.smk"
include: "rules/techs/generate_tech_files.smk"
include: "rules/spatial/generate_offshore_transmission.smk"
include: "rules/demand/generate_demand_flex.smk"


include: ""


rule all:
    input:
        "results/model_results.nc"

rule run_calliope:
    input:
        "model.yml",
        "demand_ND.yaml", ## Sens
        "cp30_constraint.yaml",
        "techs/transmission.yaml",
        "spatial/onshore_transmission.yaml",
        "spatial/offshore_transmission.yaml",
        "techs/CCS.yaml",
        "techs/fossil_fuels.yaml",
        "techs/low_carbon_thermal.yaml",
        "techs/offshore_wind.yaml",
        "techs/other_renewables.yaml",
        "techs/solar_onshore_wind.yaml",
        "techs/storage.yaml",
        "techs/import_export.yaml",
        "spatial/nodes_techs.yaml",
#        "spatial/capacities_2023.yaml",
        "techs/thermal_power_constraints.yaml"
    output:
        "results/model_results.nc"
    conda:
        "environment.yml"
    params:
        response_hrs=4, # Max flexibility window
        resolution_hrs=1,
    script:
        "scripts/run_calliope.py"
    # shell:
    #     """
    #     calliope run model.yml --save_netcdf={output}
    #     """

rule serve_calligraph:
    input:
        "results/model_results.nc"
    conda:
        "environment.yml"
    shell:
        """
        calligraph {input}
        """
        

rule dag_dot:
    output: temp("data/interim/dag.dot")
    shell: "snakemake --rulegraph > {output}"


rule rulegraph:
    message: "Plot dependency graph of the workflow."
    input: rules.dag_dot.output[0]
    # Output is deliberately omitted so rule is executed each time.
    conda: "envs/dag.yaml"
    shell: "dot -Tpdf {input} -o rulegraph.pdf"




