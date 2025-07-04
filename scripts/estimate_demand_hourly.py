import numpy as np
import pandas as pd


def create_timeseries(HIS_norm, peak_hour, AM_min_hour, PM_min_hour, P, A, M):
    """Returns an hourly timeseries for yearly data entry from FES regional breakdown based on peaks/minimums.

    Parameters:
        HIS_norm: pd.series of 2013 national demand data, indexed by hour of year (0-8759)
        peak_hour, AM_min_hour. PM_min_hour: scalars, 2013 data hours (1-24)
        P, A, M: scalars, peak and mininimum demand values from FES data (MW)

    Returns:
        pd.series of hourly demand indexed by GSP point, demand type, year (MWh)
    """
    # Correction factors at anchor points
    k_peak = P / HIS_norm.loc[peak_hour]
    k_am = A / HIS_norm.loc[AM_min_hour]
    k_pm = M / HIS_norm.loc[PM_min_hour]

    # initialize array
    k = np.zeros(len(HIS_norm))

    # sort anchor order
    anchors = [(AM_min_hour, k_am), (peak_hour, k_peak), (PM_min_hour, k_pm)]
    anchors_sorted = sorted(anchors, key=lambda x: x[0])

    # Interpolation
    for i in range(3):
        h_start, k_start = anchors_sorted[i]
        h_end, k_end = anchors_sorted[(i + 1) % 3]  # wraps around

        if h_end > h_start:
            hours = np.arange(h_start, h_end)
        else:
            # wraps around year boundary
            hours = np.concatenate(
                [np.arange(h_start, len(HIS_norm)), np.arange(0, h_end)]
            )

        interp_vals = np.linspace(k_start, k_end, len(hours), endpoint=False)
        k[hours] = interp_vals

    # Assign correction values
    for h, k_val in anchors:
        k[h] = k_val

    # Scale timeseries
    scaled = HIS_norm.values * k

    return pd.Series(scaled, index=HIS_norm.index)


# file inputs
HIS = pd.read_csv(snakemake.input["HISdemand"])
FES = pd.read_csv(snakemake.input["FESdemand"])

# convert data to hours in each year (0-8759)
HIS["SETTLEMENT_DATE"] = pd.to_datetime(HIS["SETTLEMENT_DATE"])
HIS["Hour_yearly"] = (HIS["SETTLEMENT_DATE"].dt.dayofyear - 1) * 24 + (HIS["HOUR"] - 1)
HIS = HIS.sort_values("Hour_yearly")

# normalize historical demand
HIS_peak = HIS["ND"].max()
HIS["ND_norm"] = HIS["ND"] / HIS_peak

# find peak hours in 2013 data
peak_hour = HIS.loc[HIS["ND_norm"].idxmax(), "HOUR"]
AM_HIS = HIS[HIS["HOUR"].between(1, 12)]
AM_min_hour = AM_HIS.loc[AM_HIS["ND_norm"].idxmin(), "HOUR"]
PM_HIS = HIS[HIS["HOUR"].between(13, 24)]
PM_min_hour = PM_HIS.loc[PM_HIS["ND_norm"].idxmin(), "HOUR"]

# set index for scaling function
HIS_norm_series = HIS.set_index("Hour_yearly")["ND_norm"]

results = []

# Iterates over every row in FES data
for idx, row in FES.iterrows():
    scaled_ts = create_timeseries(
        HIS_norm_series,
        peak_hour=int(peak_hour),
        AM_min_hour=int(AM_min_hour),
        PM_min_hour=int(PM_min_hour),
        P=row["DemandPk"],
        A=row["DemandAM"],
        M=row["DemandPM"],
    )
    # Reduce memory requirement, perhaps use parquet file?
    scaled_ts_float16 = scaled_ts.values.astype(np.float16)
    results.append(scaled_ts_float16)

FES["Demand"] = results
FES = FES.drop(columns=["DemandPk", "DemandAM", "DemandPM"])

print(type(results[0]))  # should be <class 'list'>
print(len(results[0]))  # should be 8760 (hours in a year)
print(results[0][:5])


FES.to_csv(snakemake.output[0])
