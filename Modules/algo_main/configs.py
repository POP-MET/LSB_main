
# Azorin-Molina, Tijm & Chen (2011) defaults
AZORIN_2011 = {
    'lt_offset_h': 10,
    'onset_window_h': (1.0, 7.5),       # filter 1, LSR +1h to +7:30h
    'ws_change': 1.5,                   # filter 2
    'ws_change_ref': 'step',            # 'step' (paper) or 'lsr'
    'ws_base_window_h': (-1.0, 0.5),    # sunrise wind, used with 'lsr'
    'onshore_sector': (45.0, 180.0),    # filter 3
    'prev_ws_change': 1.5,              # filter 4
    'prev_ws_ref': 'step',              # 'step' (paper), 'ws' or 'pre_onset'
    'pre_onset_start_h': -1.0,          # used with 'pre_onset', from LSR
    'prev_ws': 1.5,                     # used with 'ws'
    'prev_offshore_sector': (181.0, 44.0),
    'cess_window_h': (-1.0, 5.0),       # filter 5, LSS -1h to +5h
    'cess_ws': 1.5,                     # filter 6
    'cess_offshore_sector': (226.0, 44.0),
    'ws_max': 13.9,                     # filter 7
    'onshore_persist': 0.0,             # not in the paper (0 = off)
    'pres_rise': 1.0,                   # filter 8, hPa
    'pres_detide': None,                # None (paper), 'mean' or 's2'
    'use_f8': True,                     # filter 8 as a gate
    'atide_anom': -0.5,                 # filter 9, hPa
    'air_sea_dt': 0.0,                  # filter 10, C
}

# current Manus choices (issue #6 flowchart defaults)
MANUS_V0 = {
    **AZORIN_2011,
    'onset_window_h': (1.0, 7.0),
    'ws_change': 1.0,                   # weak, gradual onset at Manus
    'ws_change_ref': 'lsr',
    'onshore_sector': (0.0, 180.0),
    'prev_ws_ref': 'pre_onset',         # gradual onset at Manus
    'prev_ws': 0.5,                     # calm before onset
    'prev_offshore_sector': (181.0, 359.0),
    'cess_window_h': (-3.0, 5.0),       # earlier start, breeze can end before LSS -1h
    'cess_ws': 0.5,                     # calm at cessation
    'onshore_persist': 0.7,             # fraction of onshore steps onset -> cessation
    'use_f8': False,                    # flag only, pressure rise is the tide near the equator
    'cess_offshore_sector': (181.0, 359.0),
}
