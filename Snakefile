configfile: "config/default.yaml"

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
include: "rules/prep_calliope.smk"




rule all:
    input:
        "results/model.nc"
    default_target: True


rule run_calliope:
    input:
        model="full_model.yaml"
    output:
        "results/model.nc"
    shell:
        """
        calliope run {input.model}
        mv model.nc {output}
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




