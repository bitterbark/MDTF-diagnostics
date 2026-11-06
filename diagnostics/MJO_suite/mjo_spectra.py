# This file is part of the MJO_suite module of the MDTF code package (see LICENSE.txt)
"""Step 5 (replaces mjo_spectra.ncl and my_cam_diags.ncl mjo_wavenum_freq_season_plot_panel):
wavenumber-frequency power spectra of 10S-10N averaged anomalies, pooled over 180-day winter
(from Nov 1) and summer (from May 1) segments. Writes <CASENAME>.MJO.wave.freq.<var>.<season>.png."""
import numpy as np
import xarray as xr

import mjo_utils as U
import mjo_plots as P
from daily_anom import anom_file

VARS = ["pr", "rlut", "u200", "u850", "v200", "v850"]
SEASONS = ["winter", "summer"]
LAT_S, LAT_N = -10., 10.
# contour levels and units from mjo_spectra.ncl
LEVELS = {
    "pr": (0.1 * np.array([0.04, 0.08, 0.12, 0.16, 0.2, 0.24, 0.28, 0.32, 0.36, 0.4, 0.44, 0.48]), "(mm/day)"),
    "u850": (0.1 * np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]), "(m/s)"),
    "v850": (0.001 * np.array([0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25,
                               3.5, 3.75, 4.0]), "(m/s)"),
    "u200": (0.01 * np.array([4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68]), "(m/s)"),
    "v200": (0.001 * np.array([4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 48]), "(m/s)"),
    "rlut": (10. * np.array([0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.22, 0.24, 0.26]),
             "(W/m$^2$)"),
}
NW_PLOT, F_MIN, F_MAX = 6, -0.05, 0.05


def main():
    casename = U.env("CASENAME")
    for var in VARS:
        with xr.open_dataset(anom_file(var), use_cftime=True) as ds:
            x = ds[var].sortby("lat").sel(lat=slice(LAT_S, LAT_N)).mean("lat").load()
        dates = U.yyyymmdd(x.time)
        for season in SEASONS:
            print(f"mjo_spectra: {var} {season}")
            power, wave, freq = U.mjo_wavenum_freq_season(x.values, dates, season)
            plot(power, wave, freq, var, season, casename, dates,
                 f"{casename}.MJO.wave.freq.{var}.{season}")


def plot(power, wave, freq, var, season, casename, dates, name):
    iw = (wave >= 0) & (wave <= NW_PLOT)
    jf = (freq >= F_MIN - 1e-9) & (freq <= F_MAX + 1e-9)
    work = power[np.ix_(iw, jf)].copy()
    fsub = freq[jf]
    work[:, np.argmin(np.abs(fsub))] = work.min()          # blank the 0th frequency (as NCL)
    lev, units = LEVELS[var]
    fig, ax = P.plt.subplots(figsize=(6.5, 7))
    cmap, norm = P.precip_cmap(lev)
    cf = ax.contourf(fsub, wave[iw], work, levels=lev, cmap=cmap, norm=norm, extend="both")
    cs = ax.contour(fsub, wave[iw], work, levels=lev, colors="k", linewidths=0.6)
    ax.clabel(cs, fmt=lambda v: f"{v:g}".replace("0.", "."), fontsize=7)
    for f, style in [(1. / 30., "--"), (1. / 90., "--"), (0.0, "-")]:
        ax.plot([f, f], [0, NW_PLOT], "k" + style, lw=1.0)
    ticks = 1. / np.array([-20, -30, -90, 90, 30, 20])
    ax.set_xticks(ticks)
    ax.set_xticklabels(["20d", "30d", "90d", "90d", "30d", "20d"])
    ax.set_xlim(F_MIN, F_MAX)
    ax.set_ylim(0, NW_PLOT)
    ax.set_ylabel("zonal wavenumber")
    ax.set_title(season.upper(), fontsize=16)
    ax.text(0.03, 0.96, casename, transform=ax.transAxes, va="top", fontsize=9,
            bbox=dict(facecolor="white", edgecolor="k", lw=0.5))
    fig.suptitle(f"{var}:{dates.min()} to {dates.max()}")
    cb = fig.colorbar(cf, ax=ax, orientation="horizontal", fraction=0.05, pad=0.08, ticks=lev)
    cb.ax.set_title(units, fontsize=8)
    cb.ax.tick_params(labelsize=7)
    P.save(fig, name)


if __name__ == "__main__":
    main()
