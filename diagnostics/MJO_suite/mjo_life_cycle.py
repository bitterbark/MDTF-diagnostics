# This file is part of the MJO_suite module of the MDTF code package (see LICENSE.txt)
"""Step 7 (replaces mjo_life_cycle.ncl): composites of 20-100 day filtered OLR (color) and 850 hPa
wind (vectors) for the 8 MJO phases on days with PC1^2 + PC2^2 > 1, for May-Oct (summer) and
Nov-Apr (winter). The sign of an EOF is arbitrary, so the figures are also made with PC1 sign
flipped (.v2). Writes <CASENAME>.MJO.life.cycle.<season>[.v2].png."""
import numpy as np
import xarray as xr
import cartopy.crs as ccrs

import mjo_utils as U
import mjo_plots as P
from daily_anom import anom_file
from mjo_EOF_cal import pc_file

LAT_S, LAT_N = -20., 20.
NPHASE = 8
SEASONS = {"summer": ((5, 6, 7, 8, 9, 10), "May to Oct"), "winter": ((11, 12, 1, 2, 3, 4), "Nov to Apr")}


def read_filtered(var, wgt):
    with xr.open_dataset(anom_file(var), use_cftime=True) as ds:
        da = ds[var].sortby("lat").sel(lat=slice(LAT_S, LAT_N)).load()
    return da, U.wgt_runave(da.values.astype(np.float64), wgt, axis=0)


def main():
    casename = U.env("CASENAME")
    wgt = U.mjo_bandpass_weights()
    da, x = read_filtered("rlut", wgt)
    _, u = read_filtered("u850", wgt)
    _, v = read_filtered("v850", wgt)
    lat, lon = da.lat.values, da.lon.values
    years = [d.year for d in U._as_dates(da.time.values)]
    yr1, yr2 = min(years), max(years)

    with xr.open_dataset(pc_file(), use_cftime=True) as f:
        pc1, pc2, index = f["PC1"].values, f["PC2"].values, f["MJO_INDEX"].values
        mon = U.month_of(f.time)
    for tail, sign in (("", 1.0), (".v2", -1.0)):
        ang = np.degrees(np.arctan2(pc2, sign * pc1))
        ang = np.where(ang < 0, ang + 360., ang)
        for season, (months, label) in SEASONS.items():
            fig = P.plt.figure(figsize=(8, 12))
            ndays = []
            for n in range(NPHASE):
                lo, hi = 45. * n, 45. * (n + 1)
                with np.errstate(invalid="ignore"):
                    nt = (index > 1.0) & np.isin(mon, months) & (ang >= lo) & (ang < hi)
                ndays.append(int(nt.sum()))
                ax = P.map_axes(fig, NPHASE, 1, n + 1, LAT_S, LAT_N, 210., lon_step=60, lat_step=15)
                if nt.any():
                    xa, ua, va = (np.nanmean(a[nt], axis=0) for a in (x, u, v))
                    data, clon = P.cyclic(xa, lon)
                    cf = ax.contourf(clon, lat, data, levels=np.arange(-40, 41, 5), extend="both",
                                     cmap=P.CMAP_VIBLGRWHYEORRE, transform=ccrs.PlateCarree())
                    # thin the vectors to about every 4 deg (NCL vcMinDistanceF)
                    dlat = max(1, int(round(4.0 / abs(lat[1] - lat[0]))))
                    dlon = max(1, int(round(4.0 / abs(lon[1] - lon[0]))))
                    sk = (slice(None, None, dlat), slice(None, None, dlon))
                    q = ax.quiver(lon[sk[1]], lat[sk[0]], ua[sk], va[sk], transform=ccrs.PlateCarree(),
                                  scale=60, width=0.0012, headwidth=4)
                ax.text(0.005, 0.96, f"P{n + 1}: {ndays[-1]}", transform=ax.transAxes, va="top",
                        fontsize=7, bbox=dict(facecolor="white", edgecolor="k", lw=0.4, pad=1.5))
                if n < NPHASE - 1:
                    ax.set_xticklabels([])
            if any(ndays):
                ax.quiverkey(q, 0.93, -0.45, 2.0, "2 m/s", labelpos="E", fontproperties={"size": 7})
                fig.colorbar(cf, ax=fig.axes, orientation="horizontal", fraction=0.025, pad=0.04, aspect=50)
            fig.suptitle(f"{yr1}-{yr2}: {label}\nAnomalous: OLR, U850, V850", y=0.93)
            P.save(fig, f"{casename}.MJO.life.cycle.{season}{tail}")


if __name__ == "__main__":
    main()
