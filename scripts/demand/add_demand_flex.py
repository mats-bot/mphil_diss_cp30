import pandas as pd
import os
import shutil

def split_demand_by_flexibility_combined(demand_csv, flex_caps_df, output_folder):

    df = pd.read_csv(demand_csv)
    zones = df.columns.drop('timesteps')

    demand_type = os.path.basename(demand_csv).split('_')[1].split('.')[0]

    matching_rows = flex_caps_df[flex_caps_df['demand_type'] == demand_type]

    if not matching_rows.empty:
        total_flex_gw = matching_rows['flex_gw'].sum()

        gb_series = df[zones].sum(axis=1)
        gb_peak = gb_series.max()

        flex_share = min(total_flex_gw / gb_peak, 1.0)

        flex_df = df.copy()
        inflex_df = df.copy()

        for zone in zones:
            flex_df[zone] = df[zone] * flex_share
            inflex_df[zone] = df[zone] * (1 - flex_share)

        os.makedirs(output_folder, exist_ok=True)
        base_name = os.path.basename(demand_csv).replace('.csv', '')
        flex_df.to_csv(os.path.join(output_folder, f"{base_name}_flex.csv"), index=False)
        inflex_df.to_csv(os.path.join(output_folder, f"{base_name}_inflex.csv"), index=False)
    else:
        os.makedirs(output_folder, exist_ok=True)
        shutil.copy(demand_csv, output_folder)

def process_all_demand_files_combined(input_folder, flex_caps_df, output_folder):
    for fname in os.listdir(input_folder):
        if fname.endswith('.csv'):
            full_path = os.path.join(input_folder, fname)
            split_demand_by_flexibility_combined(full_path, flex_caps_df, output_folder)


flex_types = ['C', 'I', 'H', 'D', 'R', 'E']  

flex_caps_df = pd.read_csv(snakemake.input[1])
process_all_demand_files_combined(
    input_folder=snakemake.input[0],
    flex_caps_df=flex_caps_df,    
    output_folder=snakemake.output[0],
)           


flex_caps_df = pd.read_csv(snakemake.input[1])
process_all_demand_files_combined(
    input_folder=snakemake.input[2],
    flex_caps_df=flex_caps_df,    
    output_folder=snakemake.output[1],
)           