# This file is part of the MJO_suite module of the MDTF code package (see LICENSE.txt)
"""Step 3 (replaces mjo_EOF.ncl): single-variable EOFs of 20-100 day filtered anomalies, 30S-30N,
for Nov-Apr (winter) and May-Oct (summer). Writes <CASENAME>.MJO.EOF.<var>.<season>-0.png (maps)
and -1.png (time series)."""
import numpy as np
import xarray as xr
import cartopy.crs as ccrs

import mjo_utils as U
import mjo_plots as P
from daily_anom import anom_file

VARS = ["pr", "rlut", "u200", "u850", "v200", "v850"]
SEASONS = {"winter": (11, 12, 1, 2, 3, 4), "summer": (5, 6, 7, 8, 9, 10)}
NEOF = 4
LAT_S, LAT_N = -30.0, 30.0


def main():
    casename = U.env("CASENAME")
    wgt = U.mjo_bandpass_weights()
    for var in VARS:
        print(f"mjo_EOF: starting {var}")
        with xr.open_dataset(anom_file(var), use_cftime=True) as ds:
            da = ds[var].sel(lat=slice(min(LAT_S, LAT_N), max(LAT_S, LAT_N))).load()
            if da.lat.size == 0:   # descending latitudes
                da = ds[var].sel(lat=slice(LAT_N, LAT_S)).load()
            gw = ds["gw"].sel(lat=da.lat).values
        units = da.attrs.get("units", "")
        time = da.time
        years = [d.year for d in U._as_dates(time.values)]
        yr1, yr2 = min(years), max(years)
        x = U.wgt_runave(da.values, wgt, axis=0)                      # (time, lat, lon)
        x = x - np.nanmean(x, axis=0, keepdims=True)                  # dim_rmvmean
        mon = U.month_of(time)
        yrfrac = U.year_fraction(time)
        nlat, mlon = da.sizes["lat"], da.sizes["lon"]
        for season, months in SEASONS.items():
            sel = np.isin(mon, months)
            xs = x[sel] * gw[np.newaxis, :, np.newaxis]               # weight by gw (as NCL)
            print(f"mjo_EOF: computing {var}, {season} ({sel.sum()} days)")
            xs2 = xs.reshape(xs.shape[0], -1).T                        # (space, time)
            eof, pcvar, _ = U.eofunc(xs2, NEOF)
            ts = U.eofunc_ts(xs2, eof) / (nlat * mlon)
            eof = eof.reshape(NEOF, nlat, mlon)
            title = f"{var}: {yr1}-{yr2} {season}"
            name = f"{casename}.MJO.EOF.{var}.{season}"
            plot_maps(eof, pcvar, da.lat.values, da.lon.values, title, name + "-0")
            plot_ts(ts, pcvar, yrfrac[sel], units, title, name + "-1")


def plot_maps(eof, pcvar, lat, lon, title, name):
    lev = P.sym_levels(eof)
    fig = P.plt.figure(figsize=(8.5, 9))
    for n in range(eof.shape[0]):
        ax = P.map_axes(fig, eof.shape[0], 1, n + 1, LAT_S, LAT_N, 180.)
        data, clon = P.cyclic(eof[n], lon)
        cf = ax.contourf(clon, lat, data, levels=lev, cmap=P.CMAP_BLWHRE, extend="both",
                         transform=ccrs.PlateCarree())
        ax.set_title(f"EOF {n + 1}", loc="left", fontsize=9)
        ax.set_title(f"{pcvar[n]:5.1f}%", loc="right", fontsize=9)
        if n == 0:
            ax.set_title(title, loc="center", fontsize=10)
    fig.colorbar(cf, ax=fig.axes, orientation="horizontal", fraction=0.04, pad=0.05, aspect=40)
    P.save(fig, name)


def plot_ts(ts, pcvar, yrfrac, units, title, name):
    fig, axes = P.plt.subplots(ts.shape[0], 1, figsize=(8.5, 9), sharex=True)
    for n, ax in enumerate(axes):
        y = ts[n]
        ax.fill_between(yrfrac, 0, y, where=y > 0, color="red", interpolate=True, lw=0)
        ax.fill_between(yrfrac, 0, y, where=y < 0, color="blue", interpolate=True, lw=0)
        ax.plot(yrfrac, y, color="k", lw=0.3)
        ax.axhline(0, color="k", lw=0.5)
        ax.set_ylabel(units, fontsize=8)
        ax.set_title(f"EOF {n + 1}", loc="left", fontsize=9)
        ax.set_title(f"{pcvar[n]:5.1f}%", loc="right", fontsize=9)
        ax.tick_params(labelsize=8)
    fig.suptitle(title)
    P.save(fig, name)


if __name__ == "__main__":
    main()
