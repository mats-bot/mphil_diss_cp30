configfile: "config/default.yaml"

include: "rules/clean_demand.smk"
include: "rules/estimate_demand_timeseries.smk"

rule all:
    input:
        "data/intermediates/GSP_timeseries.csv"
    default_target: true

rule dag_dot:
    output: temp("data/interim/dag.dot")
    shell: "snakemake --rulegraph > {output}"


rule rulegraph:
    message: "Plot dependency graph of the workflow."
    input: rules.dag_dot.output[0]
    # Output is deliberately omitted so rule is executed each time.
    conda: "envs/dag.yaml"
    shell: "dot -Tpdf {input} -o rulegraph.pdf"
