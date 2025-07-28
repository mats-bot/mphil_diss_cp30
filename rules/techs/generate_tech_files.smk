weather_year = config["weather_year"]


rule generate_fossil_nuclear_capacities:
    input: 
        "data/processed/techs/fossil_nuclear_2023.csv"
    output: 
        "spatial/fossil_nuclear_2023_caps.yaml"
    conda: 
        "../../envs/data_processing.yaml"
    script:
        "../../scripts/techs/generate_fossil_nuclear_caps.py"

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
        "data/processed/techs/generation_costs.csv"
    output:
        "techs/fossil_fuels.yaml"
    conda: 
        "../../envs/data_processing.yaml"
    script:
        "../../scripts/techs/generate_fossil_defs.py"

rule generate_other_thermal_definitions:
    input: 
        "data/processed/techs/generation_costs.csv"
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
        f"data/processed/spatial/onshore_cf_{weather_year}.csv",
        f"data/processed/spatial/solar_cf_{weather_year}.csv"
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
        f"data/processed/spatial/offshore_cfs_{weather_year}.csv"
    output:
        "techs/offshore_wind.yaml"
    conda: 
        "../../envs/data_processing.yaml"
    script:
        "../../scripts/techs/generate_offshore_wind_defs.py"

rule generate_other_renewables_definitions:
    input: 
        "data/processed/techs/generation_costs.csv"
    output:
        "techs/other_renewables.yaml"
    conda: 
        "../../envs/data_processing.yaml"
    script:
        "../../scripts/techs/generate_other_ren_defs.py"
        
rule generate_storage_tech_definitions:
    input:
        "data/processed/techs/storage_costs.csv"
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



