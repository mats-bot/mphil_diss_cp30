import yaml
from collections import defaultdict

merged = defaultdict(dict)

for filename in snakemake.input:
    with open(filename) as f:
        data = yaml.safe_load(f)
        for top_key, content in data.items():
            if isinstance(content, dict):
                # merge dict content per top_key
                merged[top_key].update(content)
            else:
                # for non-dict values, just set or extend as needed
                merged[top_key] = content

with open(snakemake.output[0], "w") as f:
    yaml.dump(dict(merged), f, sort_keys=False)
