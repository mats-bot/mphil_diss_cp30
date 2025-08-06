import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv(snakemake.input[0], index_col=0)

cp30 = pd. read_excel(snakemake.input[1])

total_caps_mw = df.sum(axis=1)
total_caps_gw = total_caps_mw / 1000 # convert to GW


plt.figure(figsize=(10,6))
total_caps_gw.plot(kind='bar', color='skyblue', edgecolor='black')

plt.xlabel("Technologies")
plt.ylabel("Capacity (GW)")
plt.title("REPD connections queue (May 2025)")
plt.xticks(rotation=45, ha='right')
plt.tight_layout()

plt.savefig(snakemake.output[0], dpi=300)
plt.show()