import yaml

combined_techs = {}

for infile_path in snakemake.input:
    with open(infile_path) as infile:
        data = yaml.safe_load(infile)
        if not data:
            continue
        if "techs" in data:
            data = data["techs"]
        for tech_name, tech_info in data.items():
            if tech_name in combined_techs:
                raise ValueError(f"Duplicate tech key found: {tech_name}")
            combined_techs[tech_name] = tech_info

with open(snakemake.output[0], "w") as outfile:
    yaml.dump(combined_techs, outfile, sort_keys=False)