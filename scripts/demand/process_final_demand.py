import pandas as pd
from datetime import datetime, timedelta
import os

# Hard-coded 2030
def convert_demand_csv(input_csv, output_folder, target_year=30):
    df = pd.read_csv(input_csv)

    df_year = df[df['year'] == target_year]

    hour_cols = sorted([col for col in df.columns if col.startswith('h_')],
                       key=lambda x: int(x.split('_')[1]))

    grouped = df_year.groupby(['zone', 'type'])

    os.makedirs(output_folder, exist_ok=True)

    for (zone, dtype), group in grouped:
        row = group.iloc[0]

        start_datetime = datetime(target_year, 1, 1, 0, 0)
        datetimes = [start_datetime + timedelta(hours=i) for i in range(len(hour_cols))]
        values = row[hour_cols].values.astype(float)

        df_out = pd.DataFrame({
            'datetime': [dt.isoformat() for dt in datetimes],
            'value': values
        })

        filename = os.path.join(output_folder, f"demand_{zone}_{dtype}.csv")
        df_out.to_csv(filename, index=False)
        print(f"Saved {filename}")


convert_demand_csv(
    input_csv=snakemake.input[0],
    output_folder=snakemake.output[0]
)
