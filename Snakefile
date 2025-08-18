configfile: "config/default.yaml"

# Generate data files and definitions for model
include: "rules/demand/clean_demand.smk"
include: "rules/demand/estimate_demand_timeseries.smk"
include: "rules/demand/aggregate_demand_tzone.smk"
include: "rules/demand/generate_demand_final.smk"
include: "rules/spatial/create_zones.smk"
include: "rules/spatial/generate_onshore_transmission.smk"
include: "rules/spatial/eligible_areas.smk"
include: "rules/techs/generate_2023_system.smk"
include: "rules/techs/generate_monetary_costs.smk"
include: "rules/spatial/generate_capacity_factors.smk"
include: "rules/techs/generate_tech_files.smk"
include: "rules/spatial/generate_offshore_transmission.smk"
include: "rules/demand/generate_demand_flex.smk"
include: "rules/techs/generate_capacities.smk"
include: "rules/techs/generate_carbon_costs.smk"

# Interpret model results and data used to generate these
include: "rules/results/generate_method_plots.smk"


rule prepare_inputs:
    input:
        "demand/demand_ND.yaml", ## Sens
        "demand/demand_FFR.yaml", ## Sens
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
        "spatial/capacities_2023.yaml",
        "techs/thermal_power_constraints.yaml",
        "spatial/capacities_2030_REPD.yaml",
        "spatial/capacities_2030_ND.yaml", ## Sens
        "spatial/capacities_2030_FFR.yaml", ## Sens
    output:
        temp(touch("model_inputs.done"))

rule run_calliope:
    input:
        model_yaml = "model_B{run_number}.yml",
        other_inputs = rules.prepare_inputs.output[0]
    output:
        "results/model_results_B{run_number}.nc"
    conda:
        "envs/calliope.yaml"
    params:
        response_hrs=4, # Max flexibility window
        resolution_hrs=1,
        data_scaling=config["data_scaling"]
    script:
        "scripts/run_calliope.py"

rule run_scenarios:
    input:
        model_yaml = "model_{run_number}_{sens}.yml",
        other_inputs = rules.prepare_inputs.output[0]
    output:
        "results/model_results_{run_number}_{sens}.nc"
    conda:
        "envs/calliope.yaml"
    params:
        response_hrs=4, # Max flexibility window
        resolution_hrs=1,
        data_scaling=config["data_scaling"]
    script:
        "scripts/run_calliope.py"

rule serve_calligraph:
    input:
        "results/model_results_B2.nc"
    conda:
        "envs/calliope.yaml"
    shell:
        """
        calligraph {input}
        """

        

# snakemake --cores 1 results/model_results_S1_FFR.nc results/model_results_S2_FFR.nc results/model_results_S3_FFR.nc results/model_results_S4_FFR.nc results/model_results_S5_FFR.nc
# snakemake --cores 1 results/model_results_B2.nc results/model_results_S1_FFR.nc results/model_results_S1_ND.nc results/model_results_S2_FFR.nc results/model_results_S2_ND.nc results/model_results_S3_FFR.nc results/model_results_S3_ND.nc results/model_results_S4_FFR.nc results/model_results_S4_ND.nc results/model_results_S5_FFR.nc results/model_results_S5_ND.nc results/model_results_S6_ND.nc

# snakemake --rerun-incomplete --cores 1 results/model_results_S5_ND.nc results/model_results_S5_FFR.nc


rule dag_dot:
    output: temp("data/interim/dag.dot")
    shell: "snakemake all --rulegraph > {output}"


rule rulegraph:
    message: "Plot dependency graph of the workflow."
    input: rules.dag_dot.output[0]
    # Output is deliberately omitted so rule is executed each time.
    conda: "envs/dag.yaml"
    shell: "dot -Tpdf {input} -o rulegraph.pdf"



# snakemake --touch data/intermediates/spatial/offshore_cutouts/2017/cutout_01.nc data/intermediates/spatial/offshore_cutouts/2017/cutout_02.nc data/intermediates/spatial/offshore_cutouts/2017/cutout_04.nc data/intermediates/spatial/offshore_cutouts/2017/cutout_03.nc data/intermediates/spatial/offshore_cutouts/2017/cutout_05.nc data/intermediates/spatial/offshore_cutouts/2017/cutout_06.nc data/intermediates/spatial/offshore_cutouts/2017/cutout_07.nc data/intermediates/spatial/offshore_cutouts/2017/cutout_08.nc data/intermediates/spatial/offshore_cutouts/2017/cutout_09.nc data/intermediates/spatial/offshore_cutouts/2017/cutout_10.nc data/intermediates/spatial/offshore_cutouts/2017/cutout_11.nc data/intermediates/spatial/offshore_cutouts/2017/cutout_12.nc