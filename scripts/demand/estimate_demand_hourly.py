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
    
    scaled_ts= scaled_ts.values.astype(np.float16)
    results.append(scaled_ts)

FES['Demand'] = results
FES = FES.drop(columns=['DemandPk', 'DemandAM', 'DemandPM'])

# Store with one data point per hour of year
FES['Demand'] = FES['Demand'].apply(
    lambda x: '[' + ' '.join([f"{val:.6f}" for val in x]) + ']')

FES.to_csv(snakemake.output[0], index=False,)



# Scenario 2 config
S2_HIS = pd.read_csv(snakemake.input["S2_HISdemand"])
S2_FES = pd.read_csv(snakemake.input["FESdemand"])

S2_HIS["SETTLEMENT_DATE"] = pd.to_datetime(S2_HIS["SETTLEMENT_DATE"])
S2_HIS["Hour_yearly"] = (S2_HIS["SETTLEMENT_DATE"].dt.dayofyear - 1) * 24 + (S2_HIS["HOUR"] - 1)
S2_HIS = S2_HIS.sort_values("Hour_yearly")

S2_HIS_peak = S2_HIS["ND"].max()
S2_HIS["ND_norm"] = S2_HIS["ND"] / S2_HIS_peak

S2_peak_hour = S2_HIS.loc[HIS["ND_norm"].idxmax(), "HOUR"]
S2_AM_HIS = S2_HIS[S2_HIS["HOUR"].between(1, 12)]
S2_AM_min_hour = S2_AM_HIS.loc[S2_AM_HIS["ND_norm"].idxmin(), "HOUR"]
S2_PM_HIS = S2_HIS[S2_HIS["HOUR"].between(13, 24)]
S2_PM_min_hour = S2_PM_HIS.loc[S2_PM_HIS["ND_norm"].idxmin(), "HOUR"]

S2_HIS_norm_series = S2_HIS.set_index("Hour_yearly")["ND_norm"]

S2_results = []

for idx, row in S2_FES.iterrows():
    S2_scaled_ts = create_timeseries(
        S2_HIS_norm_series,
        peak_hour=int(S2_peak_hour),
        AM_min_hour=int(S2_AM_min_hour),
        PM_min_hour=int(S2_PM_min_hour),
        P=row["DemandPk"],
        A=row["DemandAM"],
        M=row["DemandPM"],
    )
    
    S2_scaled_ts= S2_scaled_ts.values.astype(np.float16)
    S2_results.append(S2_scaled_ts)

S2_FES['Demand'] = S2_results
S2_FES = S2_FES.drop(columns=['DemandPk', 'DemandAM', 'DemandPM'])

# Store with one data point per hour of year
S2_FES['Demand'] = S2_FES['Demand'].apply(
    lambda x: '[' + ' '.join([f"{val:.6f}" for val in x]) + ']')

S2_FES.to_csv(snakemake.output[1], index=False,)


