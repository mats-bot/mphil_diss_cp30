import pandas as pd
import matplotlib.pyplot as plt

rens_df = pd.read_csv(snakemake.input[0], index_col=0)
cp30_df = pd.read_excel(snakemake.input[1], sheet_name="CP.10")
offshore_df = pd.read_csv(snakemake.input[2], index_col=0)

tech_to_group = {
    'solar_pv': 'Solar',
    'onshore_wind': 'Onshore Wind',
    'battery': 'Storage',
    'pumped_hydro': 'Storage',
    'laes': 'Storage',
    'caes': 'Storage',
}
rens_renamed = rens_df.rename(index=tech_to_group)
rens_filtered = rens_renamed[rens_renamed.index.isin(set(tech_to_group.values()))]
rens_grouped = rens_filtered.groupby(level=0).sum()
rens_grouped_gw = rens_grouped.sum(axis=1) / 1000

offshore_wind_gw = offshore_df.iloc[1, :].sum() / 1000
rens_grouped_gw['Offshore Wind'] = offshore_wind_gw

cp30_techs = cp30_df.iloc[6:10, 12].tolist()
cp30_queue = cp30_df.iloc[6:10, 16].astype(float).values
cp30_2023 = cp30_df.iloc[6:10, 15].astype(float).values
cp30_ffr = cp30_df.iloc[6:10, 18].astype(float).values
cp30_nd = cp30_df.iloc[6:10, 19].astype(float).values

cp30_plot_df = pd.DataFrame({
    'Tech': cp30_techs,
    'Queue': cp30_queue,
    '2023 Capacity': cp30_2023,
    'FFR Capacity': cp30_ffr,
    'ND Capacity': cp30_nd
}).set_index('Tech')

rens_aligned = rens_grouped_gw.reindex(cp30_plot_df.index).fillna(0)


fig, ax = plt.subplots(figsize=(12, 7))

width = 0.15
x = range(len(cp30_plot_df))

ax.bar([p - 2*width for p in x], cp30_plot_df['2023 Capacity'], width=width, label='2023 Capacity (CP30)')
ax.bar([p - width for p in x], cp30_plot_df['FFR Capacity'], width=width, label='FFR Capacity')
ax.bar(x, cp30_plot_df['ND Capacity'], width=width, label='ND Capacity')
ax.bar([p + width for p in x], cp30_plot_df['Queue'], width=width, label='CP30 Queue')
ax.bar([p + 2*width for p in x], rens_aligned, width=width, label='REPD Queue')

ax.set_xticks(x)
ax.set_xticklabels(cp30_plot_df.index, rotation=45, ha='right')
ax.set_ylabel('Capacity (GW)')
ax.set_title('Capacity Comparison by Technology')
ax.legend()
plt.tight_layout()
plt.savefig(snakemake.output[0], dpi=300)
plt.show()