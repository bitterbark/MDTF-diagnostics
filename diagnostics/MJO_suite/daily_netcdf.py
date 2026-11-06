# This file is part of the MJO_suite module of the MDTF code package (see LICENSE.txt)
"""Step 1 (replaces daily_netcdf.ncl and utils.ncl read_model_file).

Reads the six daily input fields, keeps the analysis period and 40S-40N, fills missing days by
linear interpolation in time, puts every field on the pr grid if needed, and writes
$WORK_DIR/model/netCDF/<CASENAME>.<var>.day.nc with the variable, gw(lat) and date(time).
"""
import os

import numpy as np
import xarray as xr

import mjo_utils as U

LAT_S, LAT_N = -40.0, 40.0
VARS = ["pr", "u850", "v850", "u200", "v200", "rlut"]   # pr first: it sets the common grid


def out_dir():
    d = os.path.join(U.env("WORK_DIR"), "model", "netCDF")
    os.makedirs(d, exist_ok=True)
    return d


def day_file(var):
    return os.path.join(out_dir(), f"{U.env('CASENAME')}.{var}.day.nc")


def time_slice(time):
    """Indices of the first and last time steps in [startdate, enddate] (NCL searched with a
    24-hour tolerance; here the date part is compared)."""
    ymd = U.yyyymmdd(time)
    y1, m1, d1 = U.parse_date(U.env("startdate"))
    y2, m2, d2 = U.parse_date(U.env("enddate"), is_end=True)
    keep = np.nonzero((ymd >= y1 * 10000 + m1 * 100 + d1) & (ymd <= y2 * 10000 + m2 * 100 + d2))[0]
    if keep.size == 0:
        raise ValueError(f"MJO_suite: no data between {U.env('startdate')} and {U.env('enddate')}")
    i1, i2 = int(keep[0]), int(keep[-1])
    if np.any(np.diff(ymd[i1:i2 + 1]) <= 0):
        raise ValueError("MJO_suite daily_netcdf: dates are not monotonic increasing")
    return i1, i2


def read_var(var, i1, i2):
    """Read one field, period i1..i2 and 40S-40N, with standard dim names (time, lat, lon)."""
    tname, latname, lonname = U.env("time_coord", "time"), U.env("lat_coord", "lat"), U.env("lon_coord", "lon")
    fname = U.env(var.upper() + "_FILE")
    vname = U.env(var + "_var", var)
    if not os.path.isfile(fname):
        raise FileNotFoundError(f"MJO_suite: daily {var} input file does not exist: {fname}")
    print(f"daily_netcdf: reading {vname} from {fname}")
    ds = xr.open_dataset(fname, use_cftime=True)
    da = ds[vname].rename({tname: "time", latname: "lat", lonname: "lon"})
    da = da.transpose("time", "lat", "lon").isel(time=slice(i1, i2 + 1))
    gw_full = ds["gw"].values if "gw" in ds else U.lat_weights(ds[latname].values)
    gw = xr.DataArray(gw_full, dims="lat", coords={"lat": ds[latname].values})
    if da.lat[0] < da.lat[-1]:
        da = da.sel(lat=slice(LAT_S, LAT_N))
        gw = gw.sel(lat=slice(LAT_S, LAT_N))
    else:
        da = da.sel(lat=slice(LAT_N, LAT_S))
        gw = gw.sel(lat=slice(LAT_N, LAT_S))
    da = da.load().astype(np.float32)
    da.name = var
    return da, gw


def to_grid(da, lat, lon):
    """Bilinear interpolation to (lat, lon), cyclic in longitude (NCL linint2, fiCyclicX=True)."""
    print(f"daily_netcdf: interpolating {da.name} from {da.sizes['lon']}x{da.sizes['lat']} "
          f"to {lon.size}x{lat.size}")
    lon0 = da.lon.values
    ext = xr.concat([da.isel(lon=[-1]).assign_coords(lon=[lon0[-1] - 360.0]), da,
                     da.isel(lon=[0]).assign_coords(lon=[lon0[0] + 360.0])], dim="lon")
    return ext.interp(lat=lat, lon=lon).astype(np.float32)


def main():
    tname = U.env("time_coord", "time")
    pr_file = U.env("PR_FILE")
    print(f"daily_netcdf: reading {pr_file} for time coordinates "
          "(assuming without checking that all files have the same times)")
    with xr.open_dataset(pr_file, use_cftime=True) as ds:
        time = ds[tname].load()
    i1, i2 = time_slice(time)
    print(f"daily_netcdf: time range indices {i1} {U.yyyymmdd(time[i1:i1+1])[0]} - "
          f"{i2} {U.yyyymmdd(time[i2:i2+1])[0]}")

    grid = None
    for var in VARS:
        fout = day_file(var)
        if os.path.exists(fout):
            os.remove(fout)
        da, gw = read_var(var, i1, i2)
        da = U.fill_missing_time(da)
        if grid is None:
            grid = (da.lat, da.lon, gw)
        elif da.sizes["lon"] != grid[1].size or da.sizes["lat"] != grid[0].size:
            da = to_grid(da, grid[0], grid[1])
        gw = grid[2]
        date = xr.DataArray(U.yyyymmdd(da.time).astype(np.int32), dims="time")
        out = xr.Dataset({var: da, "gw": gw.astype(np.float64), "date": date})
        enc = {var: {"_FillValue": np.float32(9.96921e36)}, "gw": {"_FillValue": None},
               "date": {"_FillValue": None}, "lat": {"_FillValue": None}, "lon": {"_FillValue": None},
               "time": {"_FillValue": None}}
        out.to_netcdf(fout, encoding=enc)
        print(f"daily_netcdf: wrote {fout}")


if __name__ == "__main__":
    main()
