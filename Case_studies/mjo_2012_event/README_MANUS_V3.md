# Manus land–sea breeze detection (MANUS_V3)

Final detection method for the ARM Manus site (TWP C1, 2.06°S, 147.43°E), used for the Feb–Apr 2012 MJO case study. It adapts the filters of Azorin-Molina, Tijm & Chen (2011) to a weak, gradual tropical sea breeze. **Detection uses observations only.** ERA5 is used afterwards to check consistency, not to detect.

## Run

```bash
python Case_studies/mjo_2012_event/run_manus_v3.py \
    --met "data/twpmetC1.b1/twpmetC1.b1.2012*.cdf" \
    --gndrad "data/twpgndrad60sC1.b1/twpgndrad60sC1.b1.2012*.cdf" \
    --sst data/oisst_manus_2012.csv \
    --start 2012-02-15 --end 2012-04-10 \
    --out lsb_manus_v3_2012.csv
```

In Python:

```python
from algo_main import MANUS_V3, LSB_detector, load_met, load_gndrad, load_sun, load_sst_erddap_csv

met = load_met(met_glob).join(load_gndrad(gndrad_glob), how='left')
sun = load_sun('LSB_detectionv2/Datasets/sunrise_sunset_twp.csv').loc[start:end]
sst = load_sst_erddap_csv(sst_csv)
res = LSB_detector(MANUS_V3).run_filter10(met, sun, sst)
```

For 15 Feb–10 Apr 2012 this gives 25 sea-breeze days of 56 (suppressed 15/22, transition 5/8, active 5/26). The earlier driver `run_detection_v2.py` produces the same v3 result together with the older versions and the ERA5 comparison.

## Inputs

| Input | Source | Processing (`Modules/algo_main/data_prep.py`) |
|---|---|---|
| Wind speed, direction, pressure, air T | ARM MET `twpmetC1.b1` (1 min): `wspd_arith_mean`, `wspd_vec_mean`, `wdir_vec_mean`, `atmos_pressure`, `temp_mean` | `load_met`: QC flag ≠ 0 → NaN; UTC → LT (+10 h, labels shifted); 30-min means labelled by interval start; `ws` = arithmetic mean speed, `wd` = direction of the 30-min vector mean |
| Ground skin temperature | ARM GNDRAD `twpgndrad60sC1.b1`: `sfc_ir_temp` (downward-looking IRT) | `load_gndrad`: QC → NaN, K → °C, LT, 30-min mean (`skin_T`) |
| Sunrise / sunset | `LSB_detectionv2/Datasets/sunrise_sunset_twp.csv` (UTC) | `load_sun`: LT (2012 event: sunrise 06:09–06:19, sunset 18:13–18:29 LT) |
| SST | NOAA OISST v2.1 daily via ERDDAP, 3 × 3 cells 147.125–147.625°E, 2.375–1.875°S | `load_sst_erddap_csv`: daily mean of the cells |

## Filters and default settings

Day = local calendar day. Times refer to 30-min steps. Sectors are wind-from directions in degrees. Onshore component is the wind component from `coast_normal` = 130°.

| Step | What it tests | MANUS_V3 setting | Paper (AZORIN_2011) |
|---|---|---|---|
| F1 onset window | search window | `onset_window_h` = (1.0, 7.0) h after sunrise → 07:30–13:00 LT | (1.0, 7.5) |
| F2 wind rise | a step in F1 with `ws` ≥ base + `ws_change`; base = mean `ws` sunrise −1 h … +30 min | `ws_change` = 0.5 m/s, `ws_change_ref` = 'lsr', `ws_base_window_h` = (−1.0, 0.5) | 1.5 m/s vs previous step ('step') |
| F3 onshore | first F2 step with direction in `onshore_sector` = candidate onset | `onshore_sector` = (0, 180) | (45, 180) |
| F4a morning check | if onshore component around sunrise < 0.5 m/s: any step from sunrise −1 h to onset −30 min is calm (`ws` < `prev_ws`) or offshore (`prev_offshore_sector`) | `prev_ws_ref` = 'pre_onset', `pre_onset_start_h` = −1.0, `prev_ws` = 0.5 m/s, `prev_offshore_sector` = (181, 359) | step change < 1.5 m/s |
| F4b onshore-background branch | if onshore component around sunrise ≥ `f4_onshore_sunrise`: onset = start of the first 1-h period in F1 whose mean onshore component ≥ sunrise value + `f4_onshore_rise` | `f4_onshore_sunrise` = 0.5 m/s, `f4_onshore_rise` = 0.5 m/s, `f4_onset_sustain_h` = 1.0 h, `coast_normal` = 130° | off |
| F5 cessation window | search window | `cess_window_h` = (−3.0, 5.0) h around sunset → 15:30–23:00 LT | (−1.0, 5.0) |
| F6 cessation | first step with `ws` < `cess_ws` or direction in `cess_offshore_sector` | `cess_ws` = 0.5 m/s, `cess_offshore_sector` = (181, 359), `cess_sustain_h` = 0 | 1.5 m/s, (226, 44) |
| F7 strong wind | max `ws` onset → cessation ≤ `ws_max` | `ws_max` = 13.9 m/s | same |
| FP persistence (added) | fraction of onshore steps onset → cessation −30 min ≥ `onshore_persist`; interior non-onshore runs ≤ `persist_gap_h` count as onshore | `onshore_persist` = 0.7, `persist_gap_h` = 1.0 h | not in paper |
| F8 pressure | pressure rise min → cessation ≥ 1 hPa; **computed, not a gate** | `use_f8` = False, `pres_detide` = None | gate |
| F9 aerial tide | not implemented (semi-diurnal tide dominates at 2°S) | – | gate |
| F10 land–sea contrast | max `skin_T` onset → cessation − daily SST ≥ `air_sea_dt` | `air_sea_ref` = 'skin', `air_sea_dt` = 0.0 K | air T |

**Sea-breeze day:** F2, F3, F4, F6, F7, FP and F10 all pass (F8 is added only if `use_f8` = True). `first_fail` gives the first gate that failed.

## Output columns (`run_filter10`)

`f2 … f10`, `fp` gate flags; `f4_onshore_bg` (F4b used), `onshore_sunrise`; `onset`, `cessation` (LT timestamps), `duration_h`; `ws_onset`, `wd_onset`, `ws_cess`, `wd_cess`, `ws_max`, `onshore_frac`; `pres_rise`, `t_pres_min`; `sst`, `temp_max_sb` (skin T max), `air_sea_dt`; `sea_breeze`, `first_fail`.

## Why the Manus changes

- **F2 0.5 m/s:** the diurnal onshore wind at the site is only 0.4–0.5 m/s in amplitude and rises gradually. The paper's 1.5 m/s step almost never occurs.
- **F3 0–180°:** the sea lies to the north-east through south-east of the site.
- **F4b:** in the suppressed phase the easterly background is often already onshore at sunrise, so there is no calm or offshore morning. Requiring a sustained 1-h strengthening keeps onsets tied to a real change.
- **F5 −3 h:** some breezes end before sunset −1 h.
- **FP:** removes days whose onshore flow is interrupted, e.g. a background wind shift.
- **F8, F9:** pressure changes follow the semi-diurnal tide. The rise at cessation correlates with cessation hour (r = 0.90) and passes on about half of random same-clock windows.
- **F10 skin T:** 2 m air T stays within ~1 K of SST and failed mainly on rainy afternoons. Skin T measures land heating directly.

## Consistency check (not part of detection)

ERA5 island convergence index, 1000–925 hPa island box minus ocean ring, 12–16 LT > 0, with a weak-background screen: agreement 0.82, Cohen's κ 0.64 over 56 days. See `verification/` and issue #24.

## Known limitations

- Single coastal station at the east tip of Manus; weak breezes under an onshore background can be missed (suppressed phase, κ 0.26 against ERA5).
- A convective outflow can pass the filters: the wind turns onshore at onset, or offshore/calm at cessation (e.g. 17 Mar 2012, heavy rain at onset). Such days were kept because the flow is onshore.
- F7 does not constrain this event (max 6.2 m/s on sea-breeze days); a Manus threshold needs a longer record.
