# This file is part of the MJO_suite module of the MDTF code package (see LICENSE.txt)
"""Step 4 (replaces mjo_lag_lat_lon.ncl): lag correlations of 20-100 day filtered precipitation
(color) and U850 (lines) with Indian Ocean precipitation (10S-5N, 75-100E): lag-longitude
(10S-10N average) and lag-latitude (80-100E average), for winter, summer and annual.
Writes <CASENAME>.MJO.lag.lat.lon.<season>.png."""
import numpy as np
import xarray as xr
import matplotlib.ticker as mticker

import mjo_utils as U
import mjo_plots as P
from daily_anom import anom_file

SEASONS = ["winter", "summer", "annual"]
LAT_S_IO, LAT_N_IO, LON_L_IO, LON_R_IO = -10., 5., 75., 100.   # Indian Ocean base region
LAT_S_GLOBE, LAT_N_GLOBE = -30., 30.
LAT_BAND = (-10., 10.)        # lat band for the (lag, lon) panel
LON_BAND = (80., 100.)        # lon band for the (lag, lat) panel
MXLAG = 25


def read(var):
    with xr.open_dataset(anom_file(var), use_cftime=True) as ds:
        da = ds[var].sortby("lat").sel(lat=slice(LAT_S_GLOBE, LAT_N_GLOBE)).load()
    return da.astype(np.float64)


def prep(da2d, wgt):
    """dtrend then band-pass filter along time (NCL dtrend_leftdim + wgt_runave_leftdim)."""
    return U.wgt_runave(U.dtrend(da2d, axis=0), wgt, axis=0)


def main():
    casename = U.env("CASENAME")
    wgt = U.mjo_bandpass_weights()
    pr, u850 = read("pr"), read("u850")
    dates = U.yyyymmdd(pr.time)
    if not np.array_equal(dates, U.yyyymmdd(u850.time)):
        raise ValueError("mjo_lag_lat_lon: pr and u850 dates do not match")

    io = pr.sel(lat=slice(LAT_S_IO, LAT_N_IO), lon=slice(LON_L_IO, LON_R_IO))
    wy = np.cos(0.017459 * io.lat.values)                     # same constant as the NCL code
    pio = (io * wy[:, None]).sum(("lat", "lon")) / (wy.sum() * io.sizes["lon"])
    pio = prep(pio.values, wgt)

    def lon_band(da):
        return prep(da.sel(lat=slice(*LAT_BAND)).mean("lat").values, wgt)

    def lat_band(da):
        return prep(da.sel(lon=slice(*LON_BAND)).mean("lon").values, wgt)

    pr_tlon, u_tlon = lon_band(pr), lon_band(u850)
    pr_tlat, u_tlat = lat_band(pr), lat_band(u850)
    lon, lat = pr.lon.values, pr.lat.values
    period = f"{dates.min()}-{dates.max()}"
    for season in SEASONS:
        r = [U.mjo_xcor_lag(pio, y, dates, MXLAG, season) for y in (pr_tlon, u_tlon, pr_tlat, u_tlat)]
        plot(r, lon, lat, f"{season}: {period}", f"{casename}.MJO.lag.lat.lon.{season}")


def plot(r, lon, lat, title, name):
    lags = np.arange(-MXLAG, MXLAG + 1)
    lev = np.round(np.arange(-1.0, 1.0001, 0.1), 2)
    fig, axes = P.plt.subplots(2, 1, figsize=(6.5, 8.5))
    for ax, (cc, cl), x, fmt in [(axes[0], (r[0], r[1]), lon, P.lon_formatter),
                                 (axes[1], (r[2], r[3]), lat, P.lat_formatter)]:
        cc, cl = U.smth9(cc, 0.5, 0.25), U.smth9(cl, 0.5, 0.25)    # optPlot@smth9 = 0.25
        cf = ax.contourf(x, lags, cc, levels=lev, cmap=P.CMAP_BLWHRE, extend="neither")
        lines = [v for v in lev if abs(v) > 1e-6]                  # zero line removed (NCL)
        # single-color contours: matplotlib dashes the negative ones, as NCL did
        cs = ax.contour(x, lags, cl, levels=lines, colors="k", linewidths=0.8)
        ax.clabel(cs, fmt="%.1f", fontsize=6)
        ax.set_ylabel("lag (days)")
        ax.set_title("precip (color)", loc="left", fontsize=9)
        ax.set_title("U (lines)", loc="right", fontsize=9)
        ax.xaxis.set_major_formatter(mticker.FuncFormatter(fmt))
    axes[0].set_xticks(np.arange(0, 360, 60))
    axes[1].set_xticks(np.arange(-20, 21, 20))
    fig.suptitle(title)
    fig.colorbar(cf, ax=axes, orientation="horizontal", fraction=0.04, pad=0.08, ticks=lev[::2])
    P.save(fig, name)


if __name__ == "__main__":
    main()
