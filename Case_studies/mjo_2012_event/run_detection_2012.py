import sys

import pandas as pd

sys.path.append('/home/nma/EXP_1d_diurnal_mse_TWP/claude_work/ext/LSB_main/Modules')
from algo_main import LSB_detector, MANUS_V0, load_met, load_sun, load_sst_erddap_csv

MET_GLOB = '/home/nma/EXP_1d_diurnal_mse_TWP/claude_work/data/twpmetC1.b1/twpmetC1.b1.2012*.cdf'
SUN_CSV = '/home/nma/EXP_1d_diurnal_mse_TWP/claude_work/ext/LSB_main/LSB_detectionv2/Datasets/sunrise_sunset_twp.csv'
OUT = '/home/nma/EXP_1d_diurnal_mse_TWP/claude_work/outputs'

# event window (LT), same as the MSE composite
EVENT = ('2012-02-15', '2012-04-10')

met = load_met(MET_GLOB)
sun = load_sun(SUN_CSV).loc[EVENT[0]:EVENT[1]]

det = LSB_detector(MANUS_V0)
f1 = det.run_filter1(met, sun)
f1.to_csv(f'{OUT}/lsb_filter1_2012.csv')

print(f1.head())
print(f1[['n_steps', 'n_valid']].describe().loc[['count', 'min', 'mean', 'max']])
print('window LT:', f1['win_start'].dt.strftime('%H:%M').unique(), '->', f1['win_end'].dt.strftime('%H:%M').unique())

# filter 2
f2 = det.run_filter2(met, sun)
f2.to_csv(f'{OUT}/lsb_filter2_2012.csv')

# phase windows (LT days), same as the MSE composite
PHASES = {'suppressed': ('2012-02-15', '2012-03-07'),
          'transition': ('2012-03-08', '2012-03-15'),
          'active': ('2012-03-16', '2012-04-10')}


def by_phase(res, col='f2'):
    return '  '.join(f"{k[:5]} {int(res.loc[a:b, col].sum())}/{len(res.loc[a:b])}" for k, (a, b) in PHASES.items())


print('filter 2:', by_phase(f2))
for name, (a, b) in PHASES.items():
    p = f2.loc[a:b]
    print(f"  {name:10s} median ws_base {p['ws_base'].median():.2f}  dws_max {p['dws_max'].median():.2f} m/s")
print('first hit LT:', f2['first_hit'].dropna().dt.strftime('%H:%M').value_counts().sort_index().to_dict())

# threshold sensitivity
for thr in [0.5, 0.75, 1.0, 1.25, 1.5]:
    r = LSB_detector({**MANUS_V0, 'ws_change': thr}).run_filter2(met, sun)
    print(f'ws_change {thr:4.2f}:', by_phase(r))

# filter 3
f3 = det.run_filter3(met, sun)
f3.to_csv(f'{OUT}/lsb_filter3_2012.csv')
print('filter 3:', by_phase(f3, 'f3'))
print('onset candidate LT:', f3['onset_cand'].dropna().dt.strftime('%H:%M').value_counts().sort_index().to_dict())

# onshore sector sensitivity (coast normal ~130 deg from the suppressed S1 wind ellipse)
for sec in [(0.0, 180.0), (45.0, 180.0), (40.0, 220.0), (60.0, 200.0)]:
    r = LSB_detector({**MANUS_V0, 'onshore_sector': sec}).run_filter3(met, sun)
    print(f'onshore_sector {sec}:', by_phase(r, 'f3'))

# filter 4
f4 = det.run_filter4(met, sun)
f4.to_csv(f'{OUT}/lsb_filter4_2012.csv')
print('filter 4:', by_phase(f4, 'f4'))
ok = f4[f4['f4']]
print('  passed by calm', int(ok['prev_calm'].sum()), ' offshore', int(ok['prev_offshore'].sum()),
      ' both', int((ok['prev_calm'] & ok['prev_offshore']).sum()))
print('onset LT:', ok['onset'].dt.strftime('%H:%M').value_counts().sort_index().to_dict())

# filter 4 sensitivity
for ref in ['ws', 'pre_onset']:
    for ws in [0.5, 1.0]:
        r = LSB_detector({**MANUS_V0, 'prev_ws_ref': ref, 'prev_ws': ws}).run_filter4(met, sun)
        print(f'{ref:9s} prev_ws {ws}:', by_phase(r, 'f4'))
for h in [-1.0, 0.0]:
    r = LSB_detector({**MANUS_V0, 'pre_onset_start_h': h}).run_filter4(met, sun)
    print(f'pre_onset start LSR{h:+.0f}h:', by_phase(r, 'f4'))
r = LSB_detector({**MANUS_V0, 'prev_ws': 0.0}).run_filter4(met, sun)
print('offshore only:', by_phase(r, 'f4'))

# filters 5-6
f6 = det.run_filter6(met, sun)
f6.to_csv(f'{OUT}/lsb_filter6_2012.csv')
f6['f1_6'] = f6['f4'] & f6['f6']
w5 = det.filter5_LSS(met, sun['lss'].iloc[0])
print('filter 5 window LT:', w5.index[0].strftime('%H:%M'), '->', w5.index[-1].strftime('%H:%M'))
print('filter 6 alone :', by_phase(f6, 'f6'))
print('filters 1-6    :', by_phase(f6, 'f1_6'))
ok = f6[f6['f1_6']]
print('  cessation by calm', int(ok['cess_calm'].sum()), ' offshore', int(ok['cess_offshore'].sum()),
      ' both', int((ok['cess_calm'] & ok['cess_offshore']).sum()))
print('cessation LT:', ok['cessation'].dt.strftime('%H:%M').value_counts().sort_index().to_dict())
for name, (a, b) in PHASES.items():
    p = ok.loc[a:b]
    if len(p):
        print(f"  {name:10s} median cessation {p['cessation'].dt.hour.add(p['cessation'].dt.minute / 60).median():.1f} LT"
              f"  duration {p['duration_h'].median():.1f} h [{p['duration_h'].min():.1f}-{p['duration_h'].max():.1f}]")

# filter 6 sensitivity
for ws in [0.0, 0.5, 1.0, 1.5]:
    r = LSB_detector({**MANUS_V0, 'cess_ws': ws}).run_filter6(met, sun)
    r['f1_6'] = r['f4'] & r['f6']
    first = r.loc[r['f1_6'], 'cessation']
    print(f'cess_ws {ws} ({"offshore only" if ws == 0 else "calm or offshore"}):', by_phase(r, 'f1_6'),
          f' median cessation {first.dt.hour.add(first.dt.minute / 60).median():.1f} LT')
r = LSB_detector({**MANUS_V0, 'cess_offshore_sector': (400.0, 400.0)}).run_filter6(met, sun)
r['f1_6'] = r['f4'] & r['f6']
print('calm only (0.5):', by_phase(r, 'f1_6'))

# filter 5 window sensitivity
for win in [(-1.0, 5.0), (-3.0, 5.0)]:
    r = LSB_detector({**MANUS_V0, 'cess_window_h': win}).run_filter6(met, sun)
    ok = r[r['f4'] & r['f6']]
    c = ok['cessation'].dt.strftime('%H:%M')
    print(f'cess_window {win}:', by_phase(r.assign(f1_6=r['f4'] & r['f6']), 'f1_6'),
          f' at window start {int((ok["cessation"] == ok["cessation"].dt.normalize() + (sun.loc[ok.index, "lss"] - sun.loc[ok.index, "lss"].dt.normalize() + pd.Timedelta(hours=win[0])).dt.ceil("30min")).sum())}/{len(ok)}')

# filter 7 + onshore persistence
f7 = det.run_filter7(met, sun)
f7['f1_7'] = f7['f4'] & f7['f6'] & f7['f7']
f7['f1_7p'] = f7['f1_7'] & f7['fp']
f7.to_csv(f'{OUT}/lsb_filter7_persist_2012.csv')
print('filters 1-6            :', by_phase(f7.assign(x=f7['f4'] & f7['f6']), 'x'))
print('filters 1-7            :', by_phase(f7, 'f1_7'))
print('filters 1-7 + persist  :', by_phase(f7, 'f1_7p'))
ok = f7[f7['f4'] & f7['f6']]
print('cessation LT:', ok['cessation'].dt.strftime('%H:%M').value_counts().sort_index().to_dict())

# ws_max distribution during the sea breeze window, for a Manus threshold
cols = ['onset', 'cessation', 'duration_h', 'ws_max', 'onshore_frac', 'f1_7p']
for name, (a, b) in PHASES.items():
    p = ok.loc[a:b]
    print(f"  {name:10s} ws_max median {p['ws_max'].median():.1f}  range {p['ws_max'].min():.1f}-{p['ws_max'].max():.1f} m/s"
          f" | onshore_frac median {p['onshore_frac'].median():.2f}")
day_max = met['ws'].groupby(met.index.floor('D')).max().loc[EVENT[0]:EVENT[1]]
print('daily max WS all days: p50 %.1f  p90 %.1f  max %.1f m/s' % (day_max.median(), day_max.quantile(0.9), day_max.max()))
print(ok[cols].assign(onset=ok['onset'].dt.strftime('%H:%M'), cessation=ok['cessation'].dt.strftime('%H:%M')).round(2).to_string())

# persistence threshold sensitivity
for fr in [0.5, 0.6, 0.7, 0.8]:
    print(f'onshore_persist {fr}:', by_phase(f7.assign(x=f7['f1_7'] & (f7['onshore_frac'] >= fr)), 'x'))

# filter 8
for how in [None, 's2', 'mean']:
    r = LSB_detector({**MANUS_V0, 'pres_detide': how}).run_filter8(met, sun)
    base = r['f4'] & r['f6'] & r['f7'] & r['fp']
    r['f1_8'] = base & r['f8']
    # also on days failing 1-7 but with onset + cessation, to see if f8 discriminates at all
    both = r['onset'].notna() & r['cessation'].notna()
    print(f'filter 8 pres_detide={how}:', by_phase(r, 'f1_8'),
          f"| rise on 1-7 days median {r.loc[base, 'pres_rise'].median():.2f} [{r.loc[base, 'pres_rise'].min():.2f}-{r.loc[base, 'pres_rise'].max():.2f}] hPa",
          f"| pressure min at LT {r.loc[base, 't_pres_min'].dt.hour.value_counts().sort_index().to_dict()}")
    if how is None:
        f8 = r
f8.to_csv(f'{OUT}/lsb_filter8_2012.csv')

# random-time null: pressure rise from a same-length window ending at the same clock time,
# on all 56 days (does a >= 1 hPa rise happen anyway?)
base = f8['f4'] & f8['f6'] & f8['f7'] & f8['fp']
null = []
for d in base[base].index:
    on, ce = f8.at[d, 'onset'], f8.at[d, 'cessation']
    for dd in sun.index:
        shift = dd - d
        seg = met.loc[on + shift:ce + shift, 'pres'].dropna()
        if len(seg):
            null.append(seg.iloc[-1] - seg.min())
null = pd.Series(null)
print(f'null (same clock windows, all days): rise >= 1 hPa in {100 * (null >= 1).mean():.0f}% of {len(null)} windows, median {null.median():.2f} hPa')

# filter 10, OISST daily mean over the 3x3 cells around Manus
SST_CSV = '/home/nma/EXP_1d_diurnal_mse_TWP/claude_work/data/oisst_manus_2012.csv'
sst = load_sst_erddap_csv(SST_CSV)
f10 = det.run_filter10(met, sun, sst)
f10.to_csv(f'{OUT}/lsb_filter10_2012.csv')
print('SST days', int(f10['sst'].notna().sum()), '/', len(f10), ' range %.2f-%.2f C' % (f10['sst'].min(), f10['sst'].max()))
print('filters 1-7 + persist      :', by_phase(f10.assign(x=f10['first_fail'].isin(['f10', ''])), 'x'))
print('final (+ f10, f8 off)      :', by_phase(f10, 'sea_breeze'))
print('first fail:', f10['first_fail'].replace('', 'PASS').value_counts().to_dict())
cand = f10[f10['first_fail'].isin(['f10', ''])]
print(cand[['onset', 'cessation', 'temp_max_sb', 'sst', 'air_sea_dt', 'f8', 'sea_breeze']]
      .assign(onset=cand['onset'].dt.strftime('%H:%M'), cessation=cand['cessation'].dt.strftime('%H:%M')).round(2).to_string())

# threshold sensitivity
for dt in [0.0, 0.5, 1.0, 1.5, 2.0]:
    r = LSB_detector({**MANUS_V0, 'air_sea_dt': dt}).run_filter10(met, sun, sst)
    print(f'air_sea_dt {dt}:', by_phase(r, 'sea_breeze'))

# null: daily max air T (07:30-19:00 LT) minus SST on all days
tday = met['temp'].between_time('07:30', '19:00')
tmax = tday.groupby(tday.index.floor('D')).max()
dall = (tmax - sst).reindex(sun.index)
for name, (a, b) in PHASES.items():
    q = dall.loc[a:b]
    print(f'  {name:10s} all-days Tmax-SST median {q.median():.2f}  >=0 on {int((q >= 0).sum())}/{len(q)}')
