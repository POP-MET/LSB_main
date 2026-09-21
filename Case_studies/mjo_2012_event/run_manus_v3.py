import argparse
import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(HERE, '..', '..', 'Modules'))
from algo_main import MANUS_V3, LSB_detector, load_gndrad, load_met, load_sst_erddap_csv, load_sun  # noqa: E402

# Manus land-sea breeze detection, final method (MANUS_V3), observations only
# inputs: ARM MET (twpmetC1.b1), ARM GNDRAD (twpgndrad60sC1.b1), sunrise/sunset csv (UTC), OISST daily csv from ERDDAP
# output: one row per local day with every filter flag, onset, cessation, sea_breeze and first_fail
# see README_MANUS_V3.md for the filters and settings
#
# example:
#   python run_manus_v3.py --met "data/twpmetC1.b1/twpmetC1.b1.2012*.cdf" \
#       --gndrad "data/twpgndrad60sC1.b1/twpgndrad60sC1.b1.2012*.cdf" \
#       --sst data/oisst_manus_2012.csv --start 2012-02-15 --end 2012-04-10 --out lsb_manus_v3_2012.csv

p = argparse.ArgumentParser()
p.add_argument('--met', required=True, help='glob of ARM MET b1 files')
p.add_argument('--gndrad', required=True, help='glob of ARM GNDRAD 60 s files')
p.add_argument('--sst', required=True, help='OISST v2.1 ERDDAP csv (cells around the site)')
p.add_argument('--sun', default=os.path.join(HERE, '..', '..', 'LSB_detectionv2', 'Datasets', 'sunrise_sunset_twp.csv'))
p.add_argument('--start', required=True)
p.add_argument('--end', required=True)
p.add_argument('--out', default='lsb_manus_v3.csv')
a = p.parse_args()

met = load_met(a.met).join(load_gndrad(a.gndrad), how='left')      # 30-min LT: ws, wd, pres, temp, skin_T
sun = load_sun(a.sun).loc[a.start:a.end]                           # LT sunrise / sunset
sst = load_sst_erddap_csv(a.sst)                                   # daily mean SST, degC

res = LSB_detector(MANUS_V3).run_filter10(met, sun, sst)
res.to_csv(a.out)
print(f'{int(res["sea_breeze"].sum())} sea-breeze days of {len(res)}; first failing filter:')
print(res['first_fail'].replace('', 'sea breeze').value_counts().to_string())
print('saved', a.out)
