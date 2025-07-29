# Combine all yaml files

tech_files = [
    "CCS.yaml",
    "fossil_fuels.yaml",
    "low_carbon_thermal.yaml",
    "offshore_wind.yaml",
    "other_renewables.yaml",
    "solar_onshore_wind.yaml",
    "storage.yaml",
    "transmission.yaml"
]

spatial_files = [
    "fossil_nuclear_2023_caps.yaml",
    "onshore_transmission.yaml",
    "renewables_2023_caps.yaml",
    "zones.yaml",
]


rule combine_techs:
    input:
        expand("techs/{file}", file=tech_files)
    output:
        "techs.yaml"
    conda:
        "../envs/data_processing.yaml"
    script:
        "../scripts/merge_techs.py" 

rule combine_spatial:
    input:
        expand("spatial/{file}", file=spatial_files)
    output:
        "spatial.yaml"
    conda:
        "../envs/data_processing.yaml"
    script:
        "../scripts/merge_spatial.py"


rule create_model:
    input:
        techs="techs.yaml",
        spatial="spatial.yaml",
        demand="demand.yaml"
    output:
        "full_model.yaml"
    script:
        "../scripts/create_model_yaml.py"