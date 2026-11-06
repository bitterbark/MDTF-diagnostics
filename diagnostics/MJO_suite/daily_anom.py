# This file is part of the MJO_suite module of the MDTF code package (see LICENSE.txt)
"""Step 2 (replaces daily_anom.ncl and calc_utils.ncl compute_daily_anom).

Removes the daily climatology from each field and writes
$WORK_DIR/model/netCDF/<CASENAME>.<var>.day.anom.nc. Precipitation is converted to mm/day.
"""
import os

import numpy as np
import xarray as xr

import mjo_utils as U
from daily_netcdf import day_file

VARS = ["pr", "rlut", "u200", "u850", "v200", "v850"]


def anom_file(var):
    return day_file(var).replace(".day.nc", ".day.anom.nc")


def main():
    for var in VARS:
        fin, fout = day_file(var), anom_file(var)
        if not os.path.isfile(fin):
            print(f"ERROR: daily_anom: missing file that should have been made by daily_netcdf: {fin}")
            continue
        print(f"daily_anom: {var}: reading {fin}")
        with xr.open_dataset(fin, use_cftime=True) as ds:
            ds.load()
        anom = U.daily_anomaly(ds[var].astype(np.float64))
        anom = anom.drop_vars("mmdd", errors="ignore")
        anom.attrs = dict(ds[var].attrs)
        if var == "pr":
            anom = anom * 1000.0 * 86400.0
            anom.attrs.update(ds[var].attrs, units="mm/day")
        anom = anom.astype(np.float32)
        anom.name = var
        out = xr.Dataset({var: anom, "gw": ds["gw"], "date": ds["date"]})
        if os.path.exists(fout):
            os.remove(fout)
        out.to_netcdf(fout, encoding={var: {"_FillValue": np.float32(9.96921e36)},
                                      "gw": {"_FillValue": None}, "date": {"_FillValue": None}})
        print(f"daily_anom: wrote {fout}")


if __name__ == "__main__":
    main()
