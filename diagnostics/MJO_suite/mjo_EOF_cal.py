# This file is part of the MJO_suite module of the MDTF code package (see LICENSE.txt)
"""Step 6 (replaces mjo_EOF_cal.ncl): combined EOF of 20-100 day filtered OLR, U850 and U200
(20S-20N averages, each normalized by its zonal mean variance), as in Wheeler and Hendon (2004).
Writes $WORK_DIR/model/netCDF/MJO_PC_INDEX.nc with PC1, PC2 (each / its standard deviation)
and MJO_INDEX = PC1^2 + PC2^2."""
import os

import numpy as np
import xarray as xr

import mjo_utils as U
from daily_anom import anom_file
from daily_netcdf import out_dir

NEOF = 2
LAT_S, LAT_N = -20., 20.
MXLAG = 25


def pc_file():
    return os.path.join(out_dir(), "MJO_PC_INDEX.nc")


def main():
    wgt = U.mjo_bandpass_weights()
    fields = []
    for var in ("rlut", "u850", "u200"):
        with xr.open_dataset(anom_file(var), use_cftime=True) as ds:
            x = ds[var].sortby("lat").sel(lat=slice(LAT_S, LAT_N)).mean("lat").load()
        time = x.time
        f = U.wgt_runave(x.values.astype(np.float64), wgt, axis=0)    # (time, lon)
        f = f - np.nanmean(f, axis=0)                                  # dim_rmvmean
        zavg_var = np.nanmean(np.nanvar(f, axis=0, ddof=1))           # zonal mean of the variance
        fields.append(f / np.sqrt(zavg_var))
    mlon = fields[0].shape[1]
    cdata = np.concatenate([f.T for f in fields], axis=0)            # (3*mlon, time)
    eof, pcvar, _ = U.eofunc(cdata, NEOF)
    pcs = U.eofunc_ts(cdata, eof)                                     # (neof, time)
    print("mjo_EOF_cal: pcvar (%) " + " ".join(f"{p:.2f}" for p in pcvar))

    ceof_ts = [U.eofunc_ts(fields[i].T, eof[:, i * mlon:(i + 1) * mlon]) for i in range(3)]
    for n in range(NEOF):
        r = [U.escorc(ceof_ts[a][n][None], ceof_ts[b][n][None])[0] for a, b in ((0, 1), (0, 2), (1, 2))]
        print(f"neof={n}  r_rlut_u850={r[0]:4.3f}  r_rlut_u200={r[1]:4.3f}  r_u850_u200={r[2]:4.3f}")

    pc1 = pcs[0] / np.nanstd(pcs[0], ddof=1)
    pc2 = pcs[1] / np.nanstd(pcs[1], ddof=1)
    index = pc1 ** 2 + pc2 ** 2
    ngood = int(np.sum(~np.isnan(index)))
    nstrong = int(np.sum(index >= 1.0))
    print(f"nGood={ngood}   nStrong={nstrong}   nOther={ngood - nstrong}")

    fout = pc_file()
    if os.path.exists(fout):
        os.remove(fout)
    ds = xr.Dataset(
        {"MJO_INDEX": ("time", index.astype(np.float32), {"long_name": "MJO PC INDEX", "info": "(PC1^2 + PC2^2)"}),
         "PC1": ("time", pc1.astype(np.float32), {"long_name": "PC1", "info": "PC1/stddev(PC1)"}),
         "PC2": ("time", pc2.astype(np.float32), {"long_name": "PC2", "info": "PC2/stddev(PC2)"})},
        coords={"time": time}, attrs={"pcvar": pcvar.astype(np.float32)})
    fill = {"_FillValue": np.float32(9.96921e36)}
    ds.to_netcdf(fout, unlimited_dims=["time"], encoding={"MJO_INDEX": fill, "PC1": fill, "PC2": fill})
    print(f"mjo_EOF_cal: wrote {fout}")


if __name__ == "__main__":
    main()
