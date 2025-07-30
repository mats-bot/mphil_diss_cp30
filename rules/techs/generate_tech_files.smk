weather_year = config["weather_year"]

# redundant in current version
rule generate_fossil_nuclear_capacities:
    input: 
        "data/processed/techs/fossil_nuclear_2023.csv"
    output: 
        "spatial/fossil_nuclear_2023_caps.yaml"
    conda: 
        "../../envs/data_processing.yaml"
    script:
        "../../scripts/techs/generate_fossil_nuclear_caps.py"

# redundant in current version
rule generate_other_renewable_capacities:
    input: 
        "data/processed/techs/renewables_2023.csv"
    output: 
        "spatial/renewables_2023_caps.yaml"
    conda:
        "../../envs/data_processing.yaml"
    script:
        "../../scripts/techs/generate_renewables_caps.py"

        
rule generate_fossil_fuel_definitions:
    input:
        "data/processed/techs/generation_costs.csv",
        "data/processed/techs/fossil_nuclear_2023.csv"
    output:
        "techs/fossil_fuels.yaml"
    conda: 
        "../../envs/data_processing.yaml"
    script:
        "../../scripts/techs/generate_fossil_defs.py"

rule generate_other_thermal_definitions:
    input: 
        "data/processed/techs/generation_costs.csv",
        "data/processed/techs/renewables_2023.csv",
        "data/processed/techs/fossil_nuclear_2023.csv"
    output:
        "techs/low_carbon_thermal.yaml"
    conda: 
        "../../envs/data_processing.yaml"
    script:
        "../../scripts/techs/generate_thermal_defs.py"

rule generate_ccs_techs_definitions:
    input: 
        "data/processed/techs/generation_costs.csv"
    output:
        "techs/CCS.yaml"
    conda: 
        "../../envs/data_processing.yaml"
    script:
        "../../scripts/techs/generate_CCS_defs.py"

rule generate_solar_onshore_wind_definitions:
    input: 
        "data/processed/techs/generation_costs.csv",
        "data/processed/spatial/onshore_cf_2030.csv",
        "data/processed/spatial/solar_cf_2030.csv",
        "data/processed/techs/renewables_2023.csv"
    output:
        "techs/solar_onshore_wind.yaml"
    conda: 
        "../../envs/data_processing.yaml"
    script:
        "../../scripts/techs/generate_land_VRE_defs.py"

rule generate_offshore_wind_project_defs:
    input: 
        "data/processed/techs/generation_costs.csv",
        "data/processed/techs/offshore_wind_projects.csv",
        "data/processed/spatial/offshore_cfs_2030.csv"
    output:
        "techs/offshore_wind.yaml",
        "data/processed/techs/offshore_wind_projects_aggregated.csv",
        "data/processed/spatial/offshore_cfs_2030_aggregated.csv"
    conda: 
        "../../envs/data_processing.yaml"
    script:
        "../../scripts/techs/generate_offshore_wind_defs.py"

rule generate_other_renewables_definitions:
    input: 
        "data/processed/techs/generation_costs.csv",
        "data/processed/techs/renewables_2023.csv"
    output:
        "techs/other_renewables.yaml"
    conda: 
        "../../envs/data_processing.yaml"
    script:
        "../../scripts/techs/generate_other_ren_defs.py"
        
rule generate_storage_tech_definitions:
    input:
        "data/processed/techs/storage_costs.csv",
        "data/processed/techs/renewables_2023.csv"
    output:
        "techs/storage.yaml"
    conda: 
        "../../envs/data_processing.yaml"
    script:
        "../../scripts/techs/generate_storage_defs.py"

rule generate_transmission_definitions:
    output:
        "techs/transmission.yaml"
    conda:
        "../../envs/data_processing.yaml"
    script:
        "../../scripts/techs/generate_transmission_defs.py"

rule generate_nodes_techs:
    input:
        tech_files=expand("techs/{tech}.yaml", tech=["CCS", "fossil_fuels", "low_carbon_thermal", 
        "offshore_wind", "other_renewables", "solar_onshore_wind", "storage"]),  
        nodes_file="spatial/zones.yaml",
        offshore_df="data/processed/techs/offshore_wind_projects.csv"
    output:
        nodes_yaml="spatial/nodes_techs.yaml"
    conda:
        "../../envs/data_processing.yaml"
    script:
        "../../scripts/spatial/generate_nodes_techs.py"


