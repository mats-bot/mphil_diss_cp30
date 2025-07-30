import pandas as pd

def shift_year_in_csv(input_csv, output_csv, old_year=2013, new_year=2030):
    df = pd.read_csv(input_csv, index_col=0, parse_dates=True)
    new_index = df.index.map(lambda ts: ts.replace(year=new_year))
    df.index = new_index
    df.to_csv(output_csv)

old_year = 2013
new_year = 2030

shift_year_in_csv(snakemake.input[0], snakemake.output[0], old_year, new_year)
shift_year_in_csv(snakemake.input[1], snakemake.output[1], old_year, new_year)
shift_year_in_csv(snakemake.input[2], snakemake.output[2], old_year, new_year)