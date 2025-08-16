import kagglehub
from kagglehub import KaggleDatasetAdapter

FILE_PATH = f"historic_demand_year_{snakemake.params.year}.csv"

# from code in https://www.kaggle.com/datasets/albertovidalrod/electricity-consumption-uk-20092022/data?select=historic_demand_year_2013.csv
df = kagglehub.load_dataset(
    KaggleDatasetAdapter.PANDAS,
    "albertovidalrod/electricity-consumption-uk-20092022",
    FILE_PATH,
)


df.to_csv(snakemake.output[0], index=False)



S2_FILE_PATH = f"historic_demand_year_2017.csv"
df = kagglehub.load_dataset(
    KaggleDatasetAdapter.PANDAS,
    "albertovidalrod/electricity-consumption-uk-20092022",
    S2_FILE_PATH,
)
df.to_csv(snakemake.output[1], index=False)