import pandas as pd
from datetime import datetime, timedelta
import os

# Hard-coded 2030
def convert_demand_csv(input_csv, output_folder, target_year=30):
    df = pd.read_csv(input_csv)

    df_year = df[df['year'] == target_year]

    hour_cols = sorted([col for col in df.columns if col.startswith('h_')],
                       key=lambda x: int(x.split('_')[1]))

    full_year = target_year + 2000
    os.makedirs(output_folder, exist_ok=True)

    grouped = df_year.groupby('type')

    for dtype, group in grouped:
        start_datetime = datetime(full_year, 1, 1, 0, 0)
        datetimes = [start_datetime + timedelta(hours=i) for i in range(len(hour_cols))]

        zone_values = {}
        for _, row in group.iterrows():
            zone = row['zone']
            values = row[hour_cols].values.astype(float)
            zone_values[zone] = values

        df_out = pd.DataFrame(zone_values)
        df_out.insert(0, 'timesteps', [dt.isoformat() for dt in datetimes])

        filename = os.path.join(output_folder, f"demand_{dtype}.csv")
        df_out.to_csv(filename, index=False)

convert_demand_csv(
    input_csv=snakemake.input[0],
    output_folder=snakemake.output[0]
)
