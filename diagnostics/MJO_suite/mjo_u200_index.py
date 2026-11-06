# This file is part of the MJO_suite module of the MDTF code package (see LICENSE.txt)
"""Step 8 (replaces mjo.ncl): 10S-10N mean daily U200, its 20-100 day band-pass (weights from the
file bp_filter_20-100) and the 101-day running variance of the filtered series.
Writes <CASENAME>.MJO.png."""
import os

import numpy as np
import xarray as xr

import mjo_utils as U
import mjo_plots as P
from daily_netcdf import day_file


def read_weights(path):
    """Read the weights in the NCL asciiread(..,-1,'float') way: every number in the file."""
    with open(path) as f:
        return np.array([float(t) for t in f.read().split()])


def main():
    casename = U.env("CASENAME")
    with xr.open_dataset(day_file("u200"), use_cftime=True) as ds:
        lat = ds.lat.values
        sel = (lat >= -10) & (lat <= 10)
        u = ds["u200"].isel(lat=sel).load()
        gw = ds["gw"].values[sel]
    y1 = U.parse_date(U.env("startdate"))[0]
    y2 = U.parse_date(U.env("enddate"), is_end=True)[0]
    ndays = min((y2 - y1 + 1) * 365, u.sizes["time"])     # as NCL: noleap years in model
    u = u.isel(time=slice(0, ndays))
    data = ((u.mean("lon") * gw).sum("lat") / gw.sum()).values.astype(np.float64)
    yrfrac = U.year_fraction(u.time)

    wts = read_weights(os.path.join(U.env("POD_HOME"), "bp_filter_20-100"))
    filt = U.wgt_runave(data, wts)
    var = np.full(filt.size, np.nan)
    for n in range(50, filt.size - 50):
        seg = filt[n - 50:n + 51]
        if not np.isnan(seg).any():
            var[n] = np.var(seg, ddof=1)

    fig, axes = P.plt.subplots(3, 1, figsize=(8.5, 10), sharex=True)
    panels = [(data, "k", f"{casename} Mean of Daily 200 mb Zonal Wind (10S-10N)", (-15, 15), "U (m/s)"),
              (filt, "blue", "Band-pass Filtered for 20-100 Day Periods", (-6, 6), "U (m/s)"),
              (var, "red", "101-day Running Variance of Band-pass Filtered Data", (0, 12),
               "variance (m$^2$/s$^2$)")]
    for ax, (y, c, title, ylim, ylab) in zip(axes, panels):
        ax.plot(yrfrac, y, color=c, lw=0.8)
        ax.axhline(0, color="k", lw=0.5)
        ax.set_title(title, fontsize=10)
        ax.set_ylim(*ylim)
        ax.set_ylabel(ylab)
    axes[-1].set_xlim(y1, y2 + 1)
    axes[-1].set_xticks(np.arange(y1, y2 + 2, 2))
    axes[-1].set_xticks(np.arange(y1, y2 + 2, 1), minor=True)
    axes[-1].set_xlabel("Year")
    P.save(fig, f"{casename}.MJO")


if __name__ == "__main__":
    main()
