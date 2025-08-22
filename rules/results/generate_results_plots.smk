rule plot_flex_flows:
    input:
        cp30 = "data/raw/techs/CP30_workbook.xlsx",
        B1 = "results/model_results_B1.nc",
        B2 = "results/model_results_B2.nc",
        S1_ND = "results/model_results_S1_ND.nc",
        S1_FFR = "results/model_results_S1_FFR.nc",
        S2_FFR = "results/model_results_S2_FFR.nc",
        S2_ND = "results/model_results_S2_ND.nc",
        S3_FFR = "results/model_results_S3_FFR.nc",
        S3_ND = "results/model_results_S3_ND.nc",
        S4_FFR = "results/model_results_S4_FFR.nc",
        S4_ND = "results/model_results_S4_ND.nc",
        S5_FFR = "results/model_results_S5_FFR.nc",
        S5_ND = "results/model_results_S5_ND.nc",
    output:
        flows_flex_nd="results/figures/discuss/flex_flows_ND.png",
        flows_flex_ffr="results/figures/discuss/flex_flows_FFR.png"
    conda:
        "../../envs/plotting.yaml"
    script:
        "../../scripts/results/discuss/plot_flex_flows.py"


rule compare_capacities:
    input:
        cp30 = "data/raw/techs/CP30_workbook.xlsx",
        ND = "results/model_results_B1.nc",
        FFR = "results/model_results_B2.nc",
    output:
        ND="results/figures/discuss/caps_ND.png",
        FFR="results/figures/discuss/caps_FFR.png",
    conda:
        "../../envs/plotting.yaml"
    script:
        "../../scripts/results/discuss/plot_compare_capacities.py"

rule compare_generation:
    input:
        cp30 = "data/raw/techs/CP30_workbook.xlsx",
        ND = "results/model_results_B1.nc",
        FFR = "results/model_results_B2.nc",
    output:
        ND="results/figures/discuss/flow_ND.png",
        FFR="results/figures/discuss/flow_FFR.png",
    conda:
        "../../envs/plotting.yaml"
    script:
        "../../scripts/results/discuss/plot_compare_flows.py"

rule plot_network_stress_base:
    input:
        zones = "uploaded_data/tzones.gpkg",
        centroids = "uploaded_data/tzones_centroids.gpkg",
        onshore_links = "spatial/onshore_transmission.yaml",
        offshore_links = "spatial/offshore_transmission.yaml",
        ND = "results/model_results_B1.nc",
        FFR = "results/model_results_B2.nc",
    output:
        caps_comp="results/figures/discuss/network_stress_base.png",
    conda:
        "../../envs/plotting.yaml"
    script:
        "../../scripts/results/discuss/plot_network_stress_base.py"


rule sens_base_comp:
    input:
        B1 = "results/model_results_B1.nc",
        B2 = "results/model_results_B2.nc",
        S1_ND = "results/model_results_S1_ND.nc",
        S1_FFR = "results/model_results_S1_FFR.nc",
        S2_FFR = "results/model_results_S2_FFR.nc",
        S2_ND = "results/model_results_S2_ND.nc",
        S3_FFR = "results/model_results_S3_FFR.nc",
        S3_ND = "results/model_results_S3_ND.nc",
        S4_FFR = "results/model_results_S4_FFR.nc",
        S4_ND = "results/model_results_S4_ND.nc",
        S5_FFR = "results/model_results_S5_FFR.nc",
        S5_ND = "results/model_results_S5_ND.nc",
        emissions = "sensitivities/S3/emissions.yaml",
    output:
        ND = "results/figures/discuss/ND_spooder.png",
        FFR = "results/figures/discuss/FFR_spooder.png"
    conda:
        "../../envs/plotting.yaml"
    script:
        "../../scripts/results/discuss/plot_spooders.py"









rule compare_scenario_capacities:
    input:
        cp30 = "data/raw/techs/CP30_workbook.xlsx",
        S1_ND = "results/model_results_S1_ND.nc",
        S1_FFR = "results/model_results_S1_FFR.nc",
        S2_ND = "results/model_results_S2_ND.nc",
        S2_FFR = "results/model_results_S2_FFR.nc",
        S3_ND = "results/model_results_S3_ND.nc",
        S3_FFR = "results/model_results_S3_FFR.nc",
        S4_ND = "results/model_results_S4_ND.nc",
        S4_FFR = "results/model_results_S4_FFR.nc",
        S5_ND = "results/model_results_S5_ND.nc",
        S5_FFR = "results/model_results_S5_FFR.nc",
    output:
        S1_ND = "results/figures/discuss/scenarios/S1_ND_capacity.png",
        S1_FFR = "results/figures/discuss/scenarios/S1_FFR_capacity.png",
        S2_ND = "results/figures/discuss/scenarios/S2_ND_capacity.png",
        S2_FFR = "results/figures/discuss/scenarios/S2_FFR_capacity.png",
        S3_ND = "results/figures/discuss/scenarios/S3_ND_capacity.png",
        S3_FFR = "results/figures/discuss/scenarios/S3_FFR_capacity.png",
        S4_ND = "results/figures/discuss/scenarios/S4_ND_capacity.png",
        S4_FFR = "results/figures/discuss/scenarios/S4_FFR_capacity.png",
        S5_ND = "results/figures/discuss/scenarios/S5_ND_capacity.png",
        S5_FFR = "results/figures/discuss/scenarios/S5_FFR_capacity.png",
    conda:
        "../../envs/plotting.yaml",
    script:
        "../../scripts/results/discuss/plot_scen_caps.py"

rule compare_scenario_flows:
    input:
        cp30 = "data/raw/techs/CP30_workbook.xlsx",
        S1_ND = "results/model_results_S1_ND.nc",
        S1_FFR = "results/model_results_S1_FFR.nc",
        S2_ND = "results/model_results_S2_ND.nc",
        S2_FFR = "results/model_results_S2_FFR.nc",
        S3_ND = "results/model_results_S3_ND.nc",
        S3_FFR = "results/model_results_S3_FFR.nc",
        S4_ND = "results/model_results_S4_ND.nc",
        S4_FFR = "results/model_results_S4_FFR.nc",
        S5_ND = "results/model_results_S5_ND.nc",
        S5_FFR = "results/model_results_S5_FFR.nc",
    output:
        capacity_png = expand("results/figures/discuss/scenarios/{scenario}_capacity.png",
                              scenario=["S1_ND","S1_FFR","S2_ND","S2_FFR","S3_ND","S3_FFR","S4_ND","S4_FFR","S5_ND","S5_FFR"]),
    conda:
        "../../envs/plotting.yaml",
    script:
        "../../scripts/results/discuss/plot_scen_caps.py"