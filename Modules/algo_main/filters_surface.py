#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 23 11:13:56 2025

@author: nalex2023
"""
import numpy as np
import pandas as pd

from .configs import AZORIN_2011


# sector (lo, hi) in degrees, lo > hi means the sector crosses north
def in_sector(wd, sector):
    lo, hi = sector
    if lo <= hi:
        return (wd >= lo) & (wd <= hi)
    return (wd >= lo) | (wd <= hi)


class LSB_detector:

    def __init__(self, input_configs=None):
        self.cfg = {**AZORIN_2011, **(input_configs or {})}

    # filter 1: expected onset period after local sunrise
    def filter1_LSR(self, met, lsr):
        h0, h1 = self.cfg['onset_window_h']
        t0 = (lsr + pd.Timedelta(hours=h0)).ceil('30min')
        t1 = (lsr + pd.Timedelta(hours=h1)).floor('30min')
        return met.loc[t0:t1]

    def run_filter1(self, met, sun):
        rows = []
        for d in sun.index:
            win = self.filter1_LSR(met, sun.at[d, 'lsr'])
            if len(win) == 0:
                continue
            rows.append({'date': d,
                         'lsr': sun.at[d, 'lsr'],
                         'win_start': win.index[0],
                         'win_end': win.index[-1],
                         'n_steps': len(win),
                         'n_valid': int(win['ws'].notna().sum())})
        return pd.DataFrame(rows).set_index('date')

    # filter 2: wind speed increase inside the onset window
    # 'step': 30-min step to step change (paper)
    # 'lsr' : increase over the mean wind around sunrise (LSR -1h to LSR +30min)
    def ws_base_lsr(self, met, lsr):
        t = lsr.floor('30min')
        h0, h1 = self.cfg['ws_base_window_h']
        return met.loc[t + pd.Timedelta(hours=h0):t + pd.Timedelta(hours=h1), 'ws'].mean()

    def filter2_ws_change(self, met, win, lsr):
        if self.cfg['ws_change_ref'] == 'step':
            dws = met['ws'].diff().reindex(win.index)
        else:
            dws = win['ws'] - self.ws_base_lsr(met, lsr)
        return dws, win[dws >= self.cfg['ws_change']]

    def run_filter2(self, met, sun):
        rows = []
        for d in sun.index:
            lsr = sun.at[d, 'lsr']
            win = self.filter1_LSR(met, lsr)
            if len(win) == 0:
                continue
            dws, hits = self.filter2_ws_change(met, win, lsr)
            rows.append({'date': d,
                         'f2': len(hits) > 0,
                         'ws_base': self.ws_base_lsr(met, lsr),
                         'first_hit': hits.index[0] if len(hits) else pd.NaT,
                         'n_hits': len(hits),
                         'dws_max': dws.max(),
                         't_dws_max': dws.idxmax() if dws.notna().any() else pd.NaT})
        return pd.DataFrame(rows).set_index('date')

    # filter 3: wind direction onshore at the filter 2 hits
    # onset candidate = first step that passes both 2 and 3
    def filter3_wd_onshore(self, hits):
        return hits[in_sector(hits['wd'], self.cfg['onshore_sector'])]

    def run_filter3(self, met, sun):
        rows = []
        for d in sun.index:
            lsr = sun.at[d, 'lsr']
            win = self.filter1_LSR(met, lsr)
            if len(win) == 0:
                continue
            dws, hits = self.filter2_ws_change(met, win, lsr)
            onshore = self.filter3_wd_onshore(hits)
            rows.append({'date': d,
                         'f2': len(hits) > 0,
                         'f3': len(onshore) > 0,
                         'onset_cand': onshore.index[0] if len(onshore) else pd.NaT,
                         'wd_onset_cand': onshore['wd'].iloc[0] if len(onshore) else np.nan,
                         'n_onshore_hits': len(onshore),
                         'wd_first_hit': hits['wd'].iloc[0] if len(hits) else np.nan})
        return pd.DataFrame(rows).set_index('date')

    # filter 4: calm or offshore before onset
    # 'step'     : 30-min speed change of the step before onset < prev_ws_change (paper)
    # 'ws'       : wind speed of the step before onset < prev_ws
    # 'pre_onset': any step from LSR + pre_onset_start_h to the step before onset
    #              has wind speed < prev_ws or offshore direction
    def filter4_prev_ws_wd(self, met, t, lsr):
        prev = t - pd.Timedelta(minutes=30)
        if self.cfg['prev_ws_ref'] == 'pre_onset':
            t0 = lsr.floor('30min') + pd.Timedelta(hours=self.cfg['pre_onset_start_h'])
            pre = met.loc[t0:prev].dropna(subset=['ws', 'wd'])
            calm = pre['ws'] < self.cfg['prev_ws']
            offshore = in_sector(pre['wd'], self.cfg['prev_offshore_sector'])
            return bool(calm.any()), bool(offshore.any())

        if prev not in met.index or pd.isna(met.at[prev, 'ws']):
            return False, False
        if self.cfg['prev_ws_ref'] == 'step':
            calm = met['ws'].diff().get(prev, np.nan) < self.cfg['prev_ws_change']
        else:
            calm = met.at[prev, 'ws'] < self.cfg['prev_ws']
        offshore = in_sector(met.at[prev, 'wd'], self.cfg['prev_offshore_sector'])
        return bool(calm), bool(offshore)

    def run_filter4(self, met, sun):
        rows = []
        for d in sun.index:
            lsr = sun.at[d, 'lsr']
            win = self.filter1_LSR(met, lsr)
            if len(win) == 0:
                continue
            dws, hits = self.filter2_ws_change(met, win, lsr)
            onshore = self.filter3_wd_onshore(hits)

            # onset = first onshore hit whose previous step is calm or offshore
            onset, calm, offshore = pd.NaT, False, False
            for t in onshore.index:
                calm, offshore = self.filter4_prev_ws_wd(met, t, lsr)
                if calm or offshore:
                    onset = t
                    break
            prev = onset - pd.Timedelta(minutes=30) if pd.notna(onset) else pd.NaT
            rows.append({'date': d,
                         'f2': len(hits) > 0,
                         'f3': len(onshore) > 0,
                         'f4': pd.notna(onset),
                         'onset': onset,
                         'prev_calm': calm if pd.notna(onset) else False,
                         'prev_offshore': offshore if pd.notna(onset) else False,
                         'ws_prev': met.at[prev, 'ws'] if pd.notna(prev) else np.nan,
                         'wd_prev': met.at[prev, 'wd'] if pd.notna(prev) else np.nan,
                         'ws_onset': met.at[onset, 'ws'] if pd.notna(onset) else np.nan,
                         'wd_onset': met.at[onset, 'wd'] if pd.notna(onset) else np.nan})
        return pd.DataFrame(rows).set_index('date')

    # filter 5: expected cessation period around local sunset
    def filter5_LSS(self, met, lss):
        h0, h1 = self.cfg['cess_window_h']
        t0 = (lss + pd.Timedelta(hours=h0)).ceil('30min')
        t1 = (lss + pd.Timedelta(hours=h1)).floor('30min')
        return met.loc[t0:t1]

    # filter 6: cessation = first step in the filter 5 window that is calm or offshore
    def filter6_cessation(self, cwin):
        cwin = cwin.dropna(subset=['ws', 'wd'])
        calm = cwin['ws'] < self.cfg['cess_ws']
        offshore = in_sector(cwin['wd'], self.cfg['cess_offshore_sector'])
        hit = cwin[calm | offshore]
        if len(hit) == 0:
            return pd.NaT, False, False
        t = hit.index[0]
        return t, bool(calm[t]), bool(offshore[t])

    def run_filter6(self, met, sun):
        f4 = self.run_filter4(met, sun)
        rows = []
        for d in f4.index:
            cwin = self.filter5_LSS(met, sun.at[d, 'lss'])
            cess, calm, offshore = self.filter6_cessation(cwin)
            onset = f4.at[d, 'onset']
            rows.append({'date': d,
                         'f6': pd.notna(cess),
                         'cessation': cess,
                         'cess_calm': calm,
                         'cess_offshore': offshore,
                         'ws_cess': met.at[cess, 'ws'] if pd.notna(cess) else np.nan,
                         'wd_cess': met.at[cess, 'wd'] if pd.notna(cess) else np.nan,
                         'duration_h': (cess - onset) / pd.Timedelta(hours=1) if pd.notna(cess) and pd.notna(onset) else np.nan})
        return f4.join(pd.DataFrame(rows).set_index('date'))

    # filter 7: no strong synoptic wind during the sea breeze
    def filter7_ws_max(self, sb):
        ws_max = sb['ws'].max()
        return ws_max, bool(ws_max <= self.cfg['ws_max'])

    # onshore persistence (not in the paper): fraction of steps from onset
    # to the step before cessation with onshore direction
    def filter_persist(self, sb):
        wd = sb['wd'].dropna()
        if len(wd) == 0:
            return np.nan, False
        frac = float(in_sector(wd, self.cfg['onshore_sector']).mean())
        return frac, bool(frac >= self.cfg['onshore_persist'])

    def run_filter7(self, met, sun):
        res = self.run_filter6(met, sun)
        res['f7'] = False
        res['ws_max'] = np.nan
        res['fp'] = False
        res['onshore_frac'] = np.nan
        for d in res.index:
            onset, cess = res.at[d, 'onset'], res.at[d, 'cessation']
            if pd.isna(onset) or pd.isna(cess):
                continue
            sb = met.loc[onset:cess - pd.Timedelta(minutes=30)]
            res.at[d, 'ws_max'], res.at[d, 'f7'] = self.filter7_ws_max(sb)
            res.at[d, 'onshore_frac'], res.at[d, 'fp'] = self.filter_persist(sb)
        return res

    # pressure used by filter 8
    # None  : raw station pressure (paper)
    # 'mean': minus the mean 30-min time-of-day cycle over the whole record
    # 's2'  : minus the fitted 12-h harmonic only (atmospheric tide)
    def pres_for_filter8(self, met):
        p = met['pres']
        how = self.cfg['pres_detide']
        if how is None:
            return p
        tod = p.index.hour + p.index.minute / 60
        if how == 'mean':
            cyc = p.groupby(tod).transform('mean')
            return p - cyc + p.mean()
        w = 2 * np.pi * tod / 12
        X = np.column_stack([np.ones(len(p)), np.cos(w), np.sin(w)])
        ok = p.notna().values
        coef, *_ = np.linalg.lstsq(X[ok], p.values[ok], rcond=None)
        return p - (X[:, 1:] @ coef[1:])

    # filter 8: pressure rise at cessation over the minimum during the sea breeze
    def filter8_pres_rise(self, pres_sb):
        rise = pres_sb.iloc[-1] - pres_sb.min()
        return rise, bool(rise >= self.cfg['pres_rise'])

    def run_filter8(self, met, sun):
        res = self.run_filter7(met, sun)
        pres = self.pres_for_filter8(met)
        res['f8'] = False
        res['pres_rise'] = np.nan
        res['t_pres_min'] = pd.NaT
        for d in res.index:
            onset, cess = res.at[d, 'onset'], res.at[d, 'cessation']
            if pd.isna(onset) or pd.isna(cess):
                continue
            psb = pres.loc[onset:cess].dropna()
            if len(psb) == 0:
                continue
            res.at[d, 'pres_rise'], res.at[d, 'f8'] = self.filter8_pres_rise(psb)
            res.at[d, 't_pres_min'] = psb.idxmin()
        return res

    # filter 10: land warmer than the sea during the sea breeze
    # dt = max air temperature onset -> cessation minus daily SST
    def filter10_air_sea(self, sb, sst):
        dt = sb['temp'].max() - sst
        return dt, bool(dt >= self.cfg['air_sea_dt'])

    def run_filter10(self, met, sun, sst):
        res = self.run_filter8(met, sun)
        res['f10'] = False
        res['sst'] = sst.reindex(res.index).values
        res['temp_max_sb'] = np.nan
        res['air_sea_dt'] = np.nan
        for d in res.index:
            onset, cess = res.at[d, 'onset'], res.at[d, 'cessation']
            if pd.isna(onset) or pd.isna(cess) or pd.isna(res.at[d, 'sst']):
                continue
            sb = met.loc[onset:cess]
            res.at[d, 'temp_max_sb'] = sb['temp'].max()
            res.at[d, 'air_sea_dt'], res.at[d, 'f10'] = self.filter10_air_sea(sb, res.at[d, 'sst'])

        # final decision, optional filters only if switched on
        gates = ['f2', 'f3', 'f4', 'f6', 'f7', 'fp', 'f10']
        if self.cfg['use_f8']:
            gates.insert(gates.index('f10'), 'f8')
        flags = res[gates].fillna(False).astype(bool)
        res['sea_breeze'] = flags.all(axis=1)
        res['first_fail'] = [next((g for g in gates if not r[g]), '') for _, r in flags.iterrows()]
        return res
