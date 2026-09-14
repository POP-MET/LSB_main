import os
import sys

import numpy as np
import pandas as pd

sys.path.append('/home/nma/EXP_1d_diurnal_mse_TWP/claude_work/ext/LSB_main/Modules')
from algo_main import (LSB_detector, MANUS_V0, MANUS_V2, load_met, load_sun, load_sst_erddap_csv,
                       load_gndrad)

# Manus detector v2 vs v1, consistency check with the ERA5 island convergence index
# v2 = v1 + filter 10 on ground skin T (ARM GNDRAD, >= 0 K); detection uses ARM observations only
# reference (not a target): ERA5 CI 12-16 LT > 0 AND background passes the synoptic screen
# (speed <= 6 m/s and offshore component <= 3 m/s); under strong background CI is not treated as sea-breeze-like

W = '/home/nma/EXP_1d_diurnal_mse_TWP/claude_work'
OUT = f'{W}/outputs/lsb_manus_v2'
os.makedirs(OUT, exist_ok=True)
EVENT = ('2012-02-15', '2012-04-10')
PHASES = {'suppressed': ('2012-02-15', '2012-03-07'),
          'transition': ('2012-03-08', '2012-03-15'),
          'active': ('2012-03-16', '2012-04-10')}
SCREEN = {'speed_max': 6.0, 'offshore_max': 3.0}

met = load_met(f'{W}/data/twpmetC1.b1/twpmetC1.b1.2012*.cdf').join(
    load_gndrad(f'{W}/data/twpgndrad60sC1.b1/twpgndrad60sC1.b1.2012*.cdf'), how='left')
sun = load_sun(f'{W}/ext/LSB_main/LSB_detectionv2/Datasets/sunrise_sunset_twp.csv').loc[EVENT[0]:EVENT[1]]
sst = load_sst_erddap_csv(f'{W}/data/oisst_manus_2012.csv')
_b = pd.read_csv(f'{W}/outputs/era5_synoptic_metrics_2012.csv', index_col=0, parse_dates=True)
bg = pd.DataFrame({'onshore_bg': _b['onshore_day'], 'speed_bg': _b['U_day']})     # ERA5, reference only
ci = pd.read_csv(f'{W}/outputs/era5_island_convergence_index_daily_v1.csv', index_col=0, parse_dates=True)['ci_12_16']

screen_ok = (bg['speed_bg'] <= SCREEN['speed_max']) & (bg['onshore_bg'] >= -SCREEN['offshore_max'])
ref = ((ci > 0) & screen_ok.reindex(ci.index).fillna(False)).rename('ref_sb_like')

v1 = LSB_detector(MANUS_V0).run_filter10(met, sun, sst)
v2 = LSB_detector(MANUS_V2).run_filter10(met, sun, sst)
v2.to_csv(f'{OUT}/lsb_filter10_2012.csv')
# optional v2b: observation-only onshore-background mode for filter 4
v2b = LSB_detector({**MANUS_V2, 'f4_onshore_sunrise': 0.5, 'f4_onshore_rise': 0.5}).run_filter10(met, sun, sst)
v2b.to_csv(f'{OUT}/lsb_filter10_2012_v2b_f4obs.csv')


def kappa(a, b):
    a, b = np.asarray(a, bool), np.asarray(b, bool)
    po = (a == b).mean()
    pe = a.mean() * b.mean() + (1 - a.mean()) * (1 - b.mean())
    return (po - pe) / (1 - pe) if pe < 1 else np.nan


rows = []
for name, res in [('v1', v1), ('v2', v2), ('v2b', v2b)]:
    for ph, (a, b) in list(PHASES.items()) + [('all', EVENT)]:
        r = res.loc[a:b]
        sb = r['sea_breeze'].astype(bool)
        rf = ref.reindex(r.index).fillna(False).astype(bool)
        rows.append({'detector': name, 'phase': ph, 'days': len(r), 'n_sb': int(sb.sum()), 'n_ref': int(rf.sum()),
                     'both': int((sb & rf).sum()), 'sb_only': int((sb & ~rf).sum()), 'ref_only': int((~sb & rf).sum()),
                     'agreement': round(float((sb == rf).mean()), 2), 'kappa': round(float(kappa(sb, rf)), 2)})
tab = pd.DataFrame(rows)
tab.to_csv(f'{OUT}/consistency_v1_v2_vs_era5_ci.csv', index=False)
pd.set_option('display.width', 250)
print(tab.to_string(index=False))

j = v2[['sea_breeze', 'first_fail', 'onset', 'cessation']].join(v1[['sea_breeze', 'onset']], rsuffix='_v1')
j = j.join(ci).join(bg).join(ref)
j['phase'] = [next(p for p, (a, b) in PHASES.items() if a <= str(d.date()) <= b) for d in j.index]
changed = j[(j['sea_breeze'] != j['sea_breeze_v1'])]
print('\nchanged v1 -> v2:')
print(changed[['phase', 'sea_breeze_v1', 'sea_breeze', 'first_fail', 'ci_12_16', 'onshore_bg', 'speed_bg']].round(2).to_string())
dis = j[j['sea_breeze'].astype(bool) != j['ref_sb_like'].astype(bool)]
print('\nremaining disagreements v2 vs reference:')
print(dis[['phase', 'sea_breeze', 'ref_sb_like', 'first_fail', 'ci_12_16', 'onshore_bg', 'speed_bg']].round(2).to_string())
j.to_csv(f'{OUT}/v2_vs_reference_days.csv')

print('\nv2b vs v2 changes:')
k = v2b[['sea_breeze', 'first_fail', 'f4_onshore_bg', 'onshore_sunrise', 'onset']].join(v2[['sea_breeze', 'onset']], rsuffix='_v2').join(ci).join(ref)
print(k[k['sea_breeze'] != k['sea_breeze_v2']].round(2).to_string())
for thr in [0.25, 0.5, 1.0]:
    for rise in [0.25, 0.5, 1.0]:
        r = LSB_detector({**MANUS_V2, 'f4_onshore_sunrise': thr, 'f4_onshore_rise': rise}).run_filter10(met, sun, sst)
        sb = r['sea_breeze'].astype(bool); rf = ref.reindex(r.index).fillna(False).astype(bool)
        print(f'f4_onshore_sunrise {thr} rise {rise}: SB', {p: int(sb.loc[a:b].sum()) for p, (a, b) in PHASES.items()},
              f'agreement {float((sb == rf).mean()):.2f} kappa {kappa(sb, rf):.2f}')
