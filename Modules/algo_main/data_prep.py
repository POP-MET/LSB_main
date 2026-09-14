import glob

import numpy as np
import pandas as pd
import xarray as xr

MET_VARS = ['wspd_arith_mean', 'wspd_vec_mean', 'wdir_vec_mean',
            'atmos_pressure', 'temp_mean']


def load_met(met_glob, lt_offset_h=10, freq='30min'):
    files = sorted(glob.glob(met_glob))
    ds = xr.open_mfdataset(files, combine='by_coords',
                           data_vars=MET_VARS + [f'qc_{v}' for v in MET_VARS],
                           coords='minimal', compat='override').load()

    df = pd.DataFrame(index=pd.DatetimeIndex(ds['time'].values, name='time'))
    for v in MET_VARS:
        x = ds[v].values.astype(float)
        x[ds[f'qc_{v}'].values != 0] = np.nan
        df[v] = x

    # UTC -> LT (shift the labels, not the values)
    df.index = df.index + pd.Timedelta(hours=lt_offset_h)

    # vector mean direction from u, v
    wd = np.radians(df['wdir_vec_mean'])
    df['u'] = -df['wspd_vec_mean'] * np.sin(wd)
    df['v'] = -df['wspd_vec_mean'] * np.cos(wd)

    out = df.resample(freq).mean()
    out['wd'] = np.degrees(np.arctan2(-out['u'], -out['v'])) % 360
    out['ws'] = out['wspd_arith_mean']
    out['pres'] = out['atmos_pressure'] * 10.0      # kPa -> hPa
    out['temp'] = out['temp_mean']
    return out[['ws', 'wd', 'pres', 'temp']]


def load_sun(sun_csv, lt_offset_h=10):
    sun = pd.read_csv(sun_csv, parse_dates=['date', 'sunrise(UTC)', 'sunset(UTC)'])
    sun = sun.set_index('date')
    lt = pd.Timedelta(hours=lt_offset_h)
    return pd.DataFrame({'lsr': sun['sunrise(UTC)'] + lt,
                         'lss': sun['sunset(UTC)'] + lt})


def load_sst_erddap_csv(sst_csv):
    sst = pd.read_csv(sst_csv, skiprows=[1], parse_dates=['time'])
    daily = sst.groupby(sst['time'].dt.floor('D'))['sst'].mean()
    daily.index = daily.index.tz_localize(None)
    return daily
