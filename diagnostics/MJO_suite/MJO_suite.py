# This file is part of the MJO_suite module of the MDTF code package (see LICENSE.txt)
"""MJO_suite driver (Python version; replaces the NCL scripts of the AMWG/MJO CLIVAR suite).

Steps, each in its own module (the NCL script it replaces in parentheses):
  1. daily_netcdf     read inputs, period, 40S-40N, fill missing days  (daily_netcdf.ncl, utils.ncl)
  2. daily_anom       remove the daily climatology                     (daily_anom.ncl, calc_utils.ncl)
  3. mjo_EOF          single-variable seasonal EOFs                    (mjo_EOF.ncl)
  4. mjo_lag_lat_lon  lag-longitude / lag-latitude correlations         (mjo_lag_lat_lon.ncl)
  5. mjo_spectra      wavenumber-frequency spectra                      (mjo_spectra.ncl, my_cam_diags.ncl)
  6. mjo_EOF_cal      combined OLR/U850/U200 EOF -> MJO_PC_INDEX.nc     (mjo_EOF_cal.ncl)
  7. mjo_life_cycle   MJO phase composites                              (mjo_life_cycle.ncl)
  8. mjo_u200_index   U200 band-pass time series                        (mjo.ncl)
Numerical routines are in mjo_utils.py, plotting helpers in mjo_plots.py.
"""
import os
import sys
import time
import traceback

sys.path.insert(0, os.environ["POD_HOME"])

import daily_netcdf  # noqa: E402
import daily_anom  # noqa: E402
import mjo_EOF  # noqa: E402
import mjo_lag_lat_lon  # noqa: E402
import mjo_spectra  # noqa: E402
import mjo_EOF_cal  # noqa: E402
import mjo_life_cycle  # noqa: E402
import mjo_u200_index  # noqa: E402

STEPS = [
    ("OBTAINING DAILY OUTPUT", daily_netcdf, True),
    ("COMPUTING DAILY ANOMALIES", daily_anom, True),
    ("COMPUTING MJO EOF (may take a while for long time samples)", mjo_EOF, False),
    ("MJO lag plots", mjo_lag_lat_lon, False),
    ("MJO spectra", mjo_spectra, False),
    ("MJO combined EOF and PC index", mjo_EOF_cal, False),
    ("MJO life cycle composite", mjo_life_cycle, False),
    ("MJO U200 time series", mjo_u200_index, False),
]

failed = []
for title, module, required in STEPS:
    print(title, flush=True)
    t0 = time.time()
    try:
        module.main()
    except Exception:
        traceback.print_exc()
        failed.append(module.__name__)
        if required:
            print(f"ERROR: MJO_suite step {module.__name__} failed; later steps need its output.")
            break
    print(f"{module.__name__} finished in {time.time() - t0:.0f} s", flush=True)

if failed:
    print("MJO_suite.py finished with errors in: " + ", ".join(failed))
    sys.exit(1)
print("MJO_suite.py finished.")
