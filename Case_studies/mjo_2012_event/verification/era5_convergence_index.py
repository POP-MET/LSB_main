import glob
import warnings

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

matplotlib.use('Agg')
warnings.filterwarnings('ignore')

# ERA5 island convergence index for the whole Feb-Apr 2012 event
# CI = -(divergence over the Manus box minus divergence over the surrounding ocean ring 0.75-2 deg),
#      mean over 1000-925 hPa, 1e-5 s-1, positive = island convergence
# daily index = mean CI 12-16 LT; sea-breeze days from the LSB detector
# test per phase: difference of daily-index means SB - non-SB vs label permutation (10000x)

W = '/home/nma/EXP_1d_diurnal_mse_TWP/claude_work'
LT = pd.Timedelta(hours=10)
PHASES = {'suppressed': ('2012-02-15', '2012-03-07'),
          'transition': ('2012-03-08', '2012-03-15'),
          'active': ('2012-03-16', '2012-04-10')}
PH_COL = {'suppressed': '#eda100', 'transition': '#1baf7a', 'active': '#2a78d6'}
C_SB, C_NON, C_ALL = '#2a78d6', '#eb6834', '#4a5559'
LEVELS = [1000, 975, 950, 925]
BOX = dict(longitude=slice(146.375, 147.625), latitude=slice(-2.625, -1.625))
rng = np.random.default_rng(8)


def open_var(v):
    files = sorted(glob.glob(f'{W}/data/era5/plev/era5pl_{v}_manus_2012*.nc'))
    ds = xr.open_mfdataset(files, combine='by_coords', compat='override', coords='minimal')
    ds = ds.drop_vars([x for x in ['number', 'expver'] if x in ds.variables]).rename({'valid_time': 'time'})
    ds = ds.sel(pressure_level=LEVELS, latitude=slice(-0.25, -4.0), longitude=slice(146.0, 148.5))
    return ds[v].sortby('latitude')


u, v = open_var('u').load(), open_var('v').load()
R = 6.371e6
lat_r = np.radians(u['latitude'])
div = (u.differentiate('longitude') / (R * np.cos(lat_r) * np.radians(1))
       + (v * np.cos(lat_r)).differentiate('latitude') / (R * np.radians(1)) / np.cos(lat_r))
dist = np.hypot((u['longitude'] - 147.0) * np.cos(np.radians(2.05)), u['latitude'] + 2.05)
ring = (dist >= 0.75) & (dist <= 2.0)
ci = -(div.sel(**BOX).mean(['latitude', 'longitude']) - div.where(ring).mean(['latitude', 'longitude']))
ci = (ci.mean('pressure_level') * 1e5).to_series()
ci.index = ci.index + LT
ci.name = 'ci'
ci.to_frame().to_csv(f'{W}/outputs/era5_island_convergence_index_hourly.csv')

import sys
TAG = sys.argv[1] if len(sys.argv) > 1 else ''
DETF = f'{W}/outputs/lsb_manus_v1/lsb_filter10_2012.csv' if TAG == 'v1' else f'{W}/outputs/lsb_filter10_2012.csv'
det = pd.read_csv(DETF, index_col=0)
det.index = pd.to_datetime(det.index).normalize()
daily = ci.between_time('12:00', '16:00').groupby(ci.between_time('12:00', '16:00').index.floor('D')).mean().to_frame('ci_12_16')
daily['phase'] = [next((p for p, (a, b) in PHASES.items() if pd.Timestamp(a) <= d <= pd.Timestamp(b)), None) for d in daily.index]
daily['sea_breeze'] = daily.index.map(det['sea_breeze'])
daily.to_csv(f'{W}/outputs/era5_island_convergence_index_daily{"_" + TAG if TAG else ""}.csv')

rows = []
for ph in PHASES:
    d = daily[daily['phase'] == ph].dropna(subset=['sea_breeze'])
    lab = d['sea_breeze'].astype(bool).values
    x = d['ci_12_16'].values
    obs = x[lab].mean() - x[~lab].mean()
    perm = np.array([(lambda m: x[m].mean() - x[~m].mean())(rng.permutation(lab)) for _ in range(10000)])
    p = float((np.abs(perm) >= abs(obs)).mean())
    med_non = np.median(x[~lab])
    rows.append({'phase': ph, 'sb_days': int(lab.sum()), 'non_sb_days': int((~lab).sum()),
                 'ci_sb_mean': round(float(x[lab].mean()), 2), 'ci_non_mean': round(float(x[~lab].mean()), 2),
                 'diff': round(float(obs), 2), 'p_perm': round(p, 3),
                 'sb_days_above_non_median': f'{int((x[lab] > med_non).sum())}/{int(lab.sum())}',
                 'days_ci>0_sb': f'{int((x[lab] > 0).sum())}/{int(lab.sum())}',
                 'days_ci>0_non': f'{int((x[~lab] > 0).sum())}/{int((~lab).sum())}'})
tab = pd.DataFrame(rows)
tab.to_csv(f'{W}/outputs/era5_island_convergence_index_test{"_" + TAG if TAG else ""}.csv', index=False)
print(tab.to_string(index=False))

# ---- figure ----
fig = plt.figure(figsize=(15, 9.5), constrained_layout=True)
gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.15])
ax0 = fig.add_subplot(gs[0, :])
for ph, (a, b) in PHASES.items():
    ax0.axvspan(pd.Timestamp(a), pd.Timestamp(b) + pd.Timedelta(days=1), color=PH_COL[ph], alpha=0.12, lw=0)
    ax0.text(pd.Timestamp(a) + (pd.Timestamp(b) - pd.Timestamp(a)) / 2, 1.02, ph, transform=ax0.get_xaxis_transform(),
             ha='center', fontsize=10, color='0.25')
ax0.plot(ci.index, ci.values, color='0.75', lw=0.6, label='hourly')
ax0.plot(daily.index + pd.Timedelta(hours=14), daily['ci_12_16'], color=C_ALL, lw=1.4, label='daily mean 12–16 LT')
sb = daily[daily['sea_breeze'] == True]
nsb = daily[daily['sea_breeze'] == False]
ax0.scatter(sb.index + pd.Timedelta(hours=14), sb['ci_12_16'], s=38, color=C_SB, zorder=3, label='detected sea-breeze day')
ax0.scatter(nsb.index + pd.Timedelta(hours=14), nsb['ci_12_16'], s=14, color=C_NON, zorder=3, label='non-sea-breeze day')
ax0.axhline(0, color='0.3', lw=0.8)
ax0.set_ylabel('Island convergence index\n(10⁻⁵ s⁻¹, + = convergence)')
ax0.set_xlim(pd.Timestamp('2012-02-14'), pd.Timestamp('2012-04-12'))
ax0.set_ylim(np.nanpercentile(ci, 0.5), np.nanpercentile(ci, 99.5))
ax0.legend(fontsize=9, frameon=False, ncol=4, loc='upper left')
ax0.grid(color='0.93', lw=0.6)

hour_ci = ci.to_frame()
hour_ci['day'] = hour_ci.index.floor('D')
hour_ci['hour'] = hour_ci.index.hour
tr = tab.set_index('phase')
for j, (ph, (a, b)) in enumerate(PHASES.items()):
    ax = fig.add_subplot(gs[1, j], sharey=None if j == 0 else fig.axes[1])
    days = pd.date_range(a, b)
    for cls, col, sel, lw in [('all days', C_ALL, days, 1.2),
                              ('sea breeze', C_SB, [d for d in days if det.at[d, 'sea_breeze']], 2.2),
                              ('non sea breeze', C_NON, [d for d in days if not det.at[d, 'sea_breeze']], 2.2)]:
        h = hour_ci[hour_ci['day'].isin(sel)]
        piv = h.pivot_table(index='day', columns='hour', values='ci')
        m = piv.mean()
        if cls != 'all days':
            boot = np.array([piv.iloc[rng.integers(0, len(piv), len(piv))].mean().values for _ in range(2000)])
            lo, hi = np.percentile(boot, [2.5, 97.5], axis=0)
            ax.fill_between(m.index, lo, hi, color=col, alpha=0.15, lw=0)
        ax.plot(m.index, m.values, color=col, lw=lw, ls='--' if cls == 'all days' else '-', label=f'{cls} (n={len(piv)})')
    ax.axhline(0, color='0.3', lw=0.8)
    ax.axvspan(12, 16, color='0.9', zorder=0)
    ax.set_xticks(range(0, 24, 3))
    ax.set_xlim(0, 23)
    ax.set_xlabel('Hour (LT)')
    ax.grid(color='0.93', lw=0.6)
    r = tr.loc[ph]
    ax.set_title(f"{ph.capitalize()} — 12–16 LT: SB {r['ci_sb_mean']:+.2f} vs non {r['ci_non_mean']:+.2f}, p = {r['p_perm']:.3f}",
                 fontsize=10.5)
    if j == 0:
        ax.set_ylabel('Convergence index (10⁻⁵ s⁻¹)')
    ax.legend(fontsize=8.5, frameon=False, loc='upper left')
fig.suptitle('ERA5 Manus island convergence index (1000–925 hPa, island box minus surrounding ocean ring)\n'
             f'sea-breeze days: LSB detector {TAG or "v0"}; top: whole event; bottom: diurnal composite by phase (shading: bootstrap 95%; grey band: 12–16 LT daily-index window)', fontsize=12)
fig.savefig(f'{W}/figures/era5_island_convergence_index{"_" + TAG if TAG else ""}.png', dpi=140)
print('figure era5_island_convergence_index.png')
