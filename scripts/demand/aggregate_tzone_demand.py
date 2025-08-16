import pandas as pd
import numpy as np

def aggregate_demand(demand_csv, output_csv):
    demand_df = pd.read_csv(demand_csv)
    demand_df = demand_df.drop(columns=["GSP", "Latitude", "Longitude"])

    demand_df['Demand'] = demand_df['Demand'].str.strip('[]').str.split().apply(
        lambda x: np.array(x, dtype=float)
    )

    assert all(demand_df['Demand'].apply(len) == 8760), "Found arrays with incorrect length"

    # Aggregates data in terms of year, type of demand (e.g. commercial) and year
    aggregated_df = (
        demand_df.groupby(['zone', 'year', 'type'])
        ['Demand'].apply(lambda x: np.sum(np.vstack(x), axis=0))
        .reset_index()
    )

    # Appply the 4 Twh demand yearly reduction for I&C demand considered by NESO. Reduction applied 
    # uniformly (ignores peaks)
    for year, year_group in aggregated_df.groupby('year'):
        year_ic_mask = year_group['type'].isin(['I', 'C'])
        total_ic_demand = year_group.loc[year_ic_mask, 'Demand'].apply(sum).sum() / 1e6  #convert to TWh

        reduction_factor = 1 - (4 / total_ic_demand)

        ic_mask = (aggregated_df['year'] == year) & aggregated_df['type'].isin(['I', 'C'])
        aggregated_df.loc[ic_mask, 'Demand'] *= reduction_factor

        print(f"Year {year}: Reduced I&C demand by 4TWh (factor: {reduction_factor:.3f})")

    hourly_cols = [f"h_{i}" for i in range(8760)]  
    demand_cols = pd.DataFrame(aggregated_df['Demand'].tolist(), columns=hourly_cols)
    aggregated_df = pd.concat([aggregated_df.drop('Demand', axis=1), demand_cols], axis=1)

    aggregated_df.to_csv(output_csv, index=False)


aggregate_demand(snakemake.input[0], snakemake.output[0])

# S2
aggregate_demand(snakemake.input[1], snakemake.output[1])
