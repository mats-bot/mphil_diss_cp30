import pandas as pd
import yaml

df1 = pd.read_csv(snakemake.input[0])
df2 = pd.read_csv(snakemake.input[1])

df1.columns.values[2] = 'capacity'
df2.columns.values[2] = 'capacity'

df = pd.concat([df1, df2], ignore_index=True)

df["base_tech"] = df["CP30 technology"].str.lower()

# for fossil fuels _existing / _new
df["tech"] = df["base_tech"].replace({
    "diesel": "diesel_existing",
#    "coal": "coal_existing",
    "gas_ccgt": "gas_ccgt_existing",
    "gas_ocgt": "gas_ocgt_existing"
})


# Excluding techs defined elsewhere or "new" techs or coal (removed from model)
exclude_techs = ["offshore_wind", "hydro", "gas_ccgt_new", "gas_ccgt_chp_new", "gas_ocgt_new", "other_renewables", "coal"]
df = df[~df["tech"].isin([t.lower() for t in exclude_techs])]

# Techs where no new investment allowed - set maximum but no minimum
no_investment_techs = ["gas_ccgt_existing", "gas_ocgt_existing", "diesel_existing", "nuclear", "pumped_hydro", "biomass"]



pivot = df.pivot_table(index="tech", columns="zone", values="capacity", aggfunc="sum", fill_value=0)
pivot_min = pivot_min = pivot.loc[~pivot.index.isin(no_investment_techs)]
pivot_max = pivot.loc[pivot.index.isin(no_investment_techs)]

pivot_min.to_csv(snakemake.output[0])
pivot_max.to_csv(snakemake.output[1])

flow_cap_min_table = {
    "data": "data/processed/spatial/min_zonal_caps.csv",
    "rows": "techs",
    "columns": "nodes",
    "add_dims": {"parameters": "flow_cap_min"},
    "rename_dims": {
        "tech": "techs"
    }
}

flow_cap_max_table = {
    "data": "data/processed/spatial/max_zonal_caps.csv",
    "rows": "techs",
    "columns": "nodes",
    "add_dims": {"parameters": "flow_cap_max"},
    "rename_dims": {
        "tech": "techs"
    }
}

output_data = {
    "data_tables": {
        "flow_cap_min": flow_cap_min_table,
        "flow_cap_max": flow_cap_max_table,
    }
}

with open(snakemake.output[2], "w") as f:
    yaml.dump(output_data, f, sort_keys=False)