rule compare_installed_caps:
    input:
        cp30 = "data/raw/techs/CP30_workbook.xlsx",
        ND = "results/model_results_B1.nc",
        FFR = "results/model_results_B2.nc",
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
        caps_comp_nd="results/figures/discuss/capacities_ND.png",
        caps_comp_ffr="results/figures/discuss/capacities_FFR.png"
    conda:
        "../../envs/plotting.yaml"
    script:
        "../../scripts/results/discuss/plot_capacities.py"


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

