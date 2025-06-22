import sys
import pandas as pd

input_file = sys.argv[1]
output_file = sys.argv[2]

df = pd.read_excel(input_file, sheet_name="Active")

inter1 = df.to_csv(output_file, index=False)


print(f"Read {len(df)} rows from '{"Active"}' in {input_file}")
print(f"Saved cleaned data to {output_file}")