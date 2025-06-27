import kagglehub
from kagglehub import KaggleDatasetAdapter
import pandas as pd

FILE_PATH = "historic_demand_year_2013.csv"

# from code in https://www.kaggle.com/datasets/albertovidalrod/electricity-consumption-uk-20092022/data?select=historic_demand_year_2013.csv
df = kagglehub.load_dataset(
    KaggleDatasetAdapter.PANDAS,
    "albertovidalrod/electricity-consumption-uk-20092022",
    FILE_PATH,
)

df.to_csv("data/raw/2013_demand_bihourly.csv", index=False)