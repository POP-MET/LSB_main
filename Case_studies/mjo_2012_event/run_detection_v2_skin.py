import os
import sys

import numpy as np
import pandas as pd

sys.path.append('/home/nma/EXP_1d_diurnal_mse_TWP/claude_work/ext/LSB_main/Modules')
from algo_main import LSB_detector, MANUS_V0, load_met, load_sun, load_sst_erddap_csv, load_gndrad

# Manus detector v2 test: filter 10 with ground skin T (ARM GNDRAD) instead of air T,
# threshold sweep scored against the ERA5 island convergence index (daily 12-16 LT CI > 0)
# note: choosing the threshold by agreement with CI makes CI a tuning target, not a fully independent check

W = '/home/nma/EXP_1d_diurnal_mse_TWP/claude_work'
MET_GLOB = f'{W}/data/twpmetC1.b1/twpmetC1.b1.2012*.cdf'
GND_GLOB = f'{W}/data/twpgndrad60sC1.b1/twpgndrad60sC1.b1.2012*.cdf'
SUN_CSV = f'{W}/ext/LSB_main/LSB_detectionv2/Datasets/sunrise_sunset_twp.csv'
SST_CSV = f'{W}/data/oisst_manus_2012.csv'
CI_CSV = f'{W}/outputs/era5_island_convergence_index_daily_v1.csv'
OUT = f'{W}/outputs/lsb_manus_v2'
os.makedirs(OUT, exist_ok=True)
EVENT = ('2012-02-15', '2012-04-10')
PHASES = {'suppressed': ('2012-02-15', '2012-03-07'),
          'transition': ('2012-03-08', '2012-03-15'),
          'active': ('2012-03-16', '2012-04-10')}

met = load_met(MET_GLOB).join(load_gndrad(GND_GLOB), how='left')
sun = load_sun(SUN_CSV).loc[EVENT[0]:EVENT[1]]
sst = load_sst_erddap_csv(SST_CSV)
ci = pd.read_csv(CI_CSV, index_col=0, parse_dates=True)['ci_12_16']
conv = (ci > 0)


def kappa(a, b):
    a, b = np.asarray(a, bool), np.asarray(b, bool)
    po = (a == b).mean()
    pe = a.mean() * b.mean() + (1 - a.mean()) * (1 - b.mean())
    return (po - pe) / (1 - pe) if pe < 1 else np.nan


def score(res):
    rows = {}
    for ph, (a, b) in list(PHASES.items()) + [('all', EVENT)]:
        r = res.loc[a:b]
        c = conv.reindex(r.index).astype(bool)
        sb = r['sea_breeze'].astype(bool)
        rows[ph] = {'n_sb': int(sb.sum()), 'n_days': len(r),
                    'sb_with_ci>0': f'{int((sb & c).sum())}/{int(sb.sum())}',
                    'ci>0_detected': f'{int((sb & c).sum())}/{int(c.sum())}',
                    'agreement': round(float((sb == c).mean()), 2), 'kappa': round(float(kappa(sb, c)), 2),
                    'ci_sb': round(float(ci.reindex(r.index)[sb].mean()), 2) if sb.any() else np.nan,
                    'ci_non': round(float(ci.reindex(r.index)[~sb].mean()), 2) if (~sb).any() else np.nan}
    return rows


rows = []
base = LSB_detector(MANUS_V0).run_filter10(met, sun, sst)
for ph, v in score(base).items():
    rows.append({'filter10': 'air (v1)', 'threshold_K': MANUS_V0['air_sea_dt'], 'phase': ph, **v})
skin_all = None
for thr in [0, 4, 6, 8, 10, 12, 14]:
    res = LSB_detector({**MANUS_V0, 'air_sea_ref': 'skin', 'air_sea_dt': float(thr)}).run_filter10(met, sun, sst)
    if thr == 0:
        skin_all = res
    for ph, v in score(res).items():
        rows.append({'filter10': 'skin', 'threshold_K': thr, 'phase': ph, **v})
tab = pd.DataFrame(rows)
tab.to_csv(f'{OUT}/filter10_skin_threshold_vs_era5_ci.csv', index=False)
pd.set_option('display.width', 250)
print(tab.to_string(index=False))
skin_all[['sea_breeze', 'first_fail', 'onset', 'cessation', 'sst', 'temp_max_sb', 'air_sea_dt']].join(ci).to_csv(f'{OUT}/lsb_filter10_skin0_2012.csv')
print('\nskin T max (onset-cessation) - SST on days reaching filter 10:')
cand = skin_all[skin_all['first_fail'].isin(['f10', ''])]
print(cand[['air_sea_dt']].join(ci).assign(phase=[next(p for p, (a, b) in PHASES.items() if a <= str(d.date()) <= b) for d in cand.index]).round(2).to_string())
