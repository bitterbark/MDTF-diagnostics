# This file is part of the MJO_suite module of the MDTF code package (see LICENSE.txt)
"""Numerical helpers for the Python version of the MJO_suite POD.

Each function reproduces the NCL built-in or diagnostics_cam.ncl routine named in its docstring,
so the Python POD gives the same numbers as the NCL version it replaces. Python ports of several
of these routines were started in https://github.com/cyschneck/geocat-research/tree/main/mjo
(mjo_wavenum_freq_season, resolveWavesHayashi, the taper and the Lanczos weights); the versions
here were checked against NCL 6.6.2 output (see doc/ncl_to_python.txt).

Conventions: arrays are numpy, time is the FIRST (leftmost) axis unless stated otherwise, and
missing values are NaN.
"""
import os

import cftime
import numpy as np
import scipy.signal
import xarray as xr

# ----------------------------------------------------------------------------------------------
# environment / dates
# ----------------------------------------------------------------------------------------------


def env(name, default=None):
    """os.environ[name], with a clear error if a required variable is missing."""
    val = os.environ.get(name, default)
    if val is None:
        raise KeyError(f"MJO_suite: required environment variable {name} is not set")
    return val


def yyyymmdd(time):
    """Integer YYYYMMDD for each time value (NCL cd_calendar(time,-2)); works for any calendar."""
    t = time.values if hasattr(time, 'values') else np.asarray(time)
    return np.array([d.year * 10000 + d.month * 100 + d.day for d in _as_dates(t)], dtype=np.int64)


def year_fraction(time):
    """Decimal year for each time (NCL cd_calendar(time,4)), using the data's own calendar."""
    out = []
    for d in _as_dates(time.values if hasattr(time, 'values') else time):
        cal = getattr(d, 'calendar', 'standard') or 'standard'
        start = cftime.datetime(d.year, 1, 1, calendar=cal)
        end = cftime.datetime(d.year + 1, 1, 1, calendar=cal)
        dd = cftime.datetime(d.year, d.month, d.day, d.hour, d.minute, d.second, calendar=cal)
        out.append(d.year + (dd - start).total_seconds() / (end - start).total_seconds())
    return np.array(out)


def _as_dates(t):
    """Convert numpy datetime64 or cftime arrays to a list of objects with year/month/day."""
    t = np.asarray(t)
    if np.issubdtype(t.dtype, np.datetime64):
        return [cftime.DatetimeProlepticGregorian(*_split64(x)) for x in t]
    return list(t)


def _split64(x):
    ts = np.datetime64(x, 's').astype(object)
    return ts.year, ts.month, ts.day, ts.hour, ts.minute, ts.second


def parse_date(s, is_end=False):
    """MDTF startdate/enddate string (YYYY, YYYYMMDD, YYYYMMDDHH[MMSS]) -> (year, month, day)."""
    s = str(s).strip()
    if len(s) == 4:
        return (int(s), 12, 31) if is_end else (int(s), 1, 1)
    if len(s) >= 8:
        return int(s[0:4]), int(s[4:6]), int(s[6:8])
    raise ValueError(f"MJO_suite: cannot parse date '{s}' (expect YYYY or YYYYMMDD...)")


def month_of(time):
    return np.array([d.month for d in _as_dates(time.values if hasattr(time, 'values') else time)])


# ----------------------------------------------------------------------------------------------
# filters
# ----------------------------------------------------------------------------------------------


def filwgts_lanczos(nwt, ihp, fca, fcb=None, nsigma=1.0):
    """NCL filwgts_lanczos: Lanczos (Duchon 1979) filter weights.

    ihp = 0 low-pass (cutoff fca), 1 high-pass (cutoff fca), 2 band-pass (fca < fcb).
    Frequencies are in cycles per time step. Returns nwt symmetric weights (nwt odd).
    """
    if nwt % 2 != 1 or nwt < 3:
        raise ValueError("filwgts_lanczos: nwt must be odd and >= 3")
    nw = nwt // 2
    k = np.arange(1, nw + 1, dtype=np.float64)
    # sigma factors; NCL uses nw+1 (not nw) here, so the end weights are not forced to zero
    sigma = (np.sin(np.pi * k / (nw + 1)) / (np.pi * k / (nw + 1))) ** nsigma

    def lowpass(fc):
        w = np.empty(nwt)
        w[nw] = 2.0 * fc
        half = np.sin(2.0 * np.pi * fc * k) / (np.pi * k) * sigma
        w[nw + 1:] = half
        w[:nw] = half[::-1]
        return w / w.sum()          # NCL normalizes each low-pass filter to unit response at f=0

    if ihp == 0:
        return lowpass(fca)
    if ihp == 1:
        w = -lowpass(fca)
        w[nw] += 1.0
        return w
    if ihp == 2:
        return lowpass(fcb) - lowpass(fca)
    raise ValueError("filwgts_lanczos: ihp must be 0, 1 or 2")


def mjo_bandpass_weights():
    """The 20-100 day band-pass weights used throughout MJO_suite (201 weights, sigma = 1)."""
    return filwgts_lanczos(201, 2, 1. / 100., 1. / 20., 1.0)


def wgt_runave(x, wgt, axis=0):
    """NCL wgt_runave / wgt_runave_leftdim with opt = 0 along `axis` (default time axis 0).

    out(t) = sum_k wgt(k) x(t + k - nw). The first and last nw = len(wgt)//2 values are NaN, as is
    any value whose window contains a NaN.
    """
    x = np.asarray(x, dtype=np.float64)
    wgt = np.asarray(wgt, dtype=np.float64)
    nwt = wgt.size
    nw = nwt // 2
    x = np.moveaxis(x, axis, 0)
    nt = x.shape[0]
    out = np.full(x.shape, np.nan)
    if nt < nwt:
        return np.moveaxis(out, 0, axis)
    shp = x.shape
    x2 = x.reshape(nt, -1)
    out2 = out.reshape(nt, -1)
    k = wgt[::-1].reshape(nwt, 1)
    # work on blocks of columns to bound the FFT memory
    step = max(1, int(2e7 // max(nt, 1)))
    for c0 in range(0, x2.shape[1], step):
        xc = x2[:, c0:c0 + step]
        bad = np.isnan(xc)
        # 'valid' convolution: output index j <-> center t = j + nw
        conv = scipy.signal.fftconvolve(np.where(bad, 0.0, xc), k, mode='valid', axes=0)
        if bad.any():
            nbad = scipy.signal.fftconvolve(bad.astype(np.float64), np.ones((nwt, 1)), mode='valid', axes=0)
            conv[nbad > 0.5] = np.nan
        out2[nw:nt - nw, c0:c0 + step] = conv
    return np.moveaxis(out2.reshape(shp), 0, axis)


def runave(x, nave, axis=-1):
    """NCL runave(x, nave, 0): unweighted running mean, ends (nave//2 each side) NaN."""
    return wgt_runave(x, np.full(nave, 1.0 / nave), axis=axis)


def dtrend(x, axis=0):
    """NCL dtrend / dtrend_leftdim (no missing values): remove the least-squares line (and mean)."""
    return scipy.signal.detrend(np.asarray(x, dtype=np.float64), axis=axis, type='linear')


def taper(x, pct=0.1, axis=0):
    """NCL taper / taper_leftdim with iopt = 0: split cosine bell toward the series mean."""
    x = np.moveaxis(np.array(x, dtype=np.float64), axis, 0)
    n = x.shape[0]
    m = max(1, int(pct * n + 0.5) // 2)
    i = np.arange(1, m + 1)
    w = 0.5 - 0.5 * np.cos(np.pi / m * (i - 0.5))
    w = w.reshape((m,) + (1,) * (x.ndim - 1))
    mean = x.mean(axis=0, keepdims=True)
    # NCL tapers toward zero instead when the series is constant
    const = np.all(x == x[0:1], axis=0, keepdims=True)
    mean = np.where(const, 0.0, mean)
    x[:m] = (x[:m] - mean) * w + mean
    x[n - m:][::-1] = (x[n - m:][::-1] - mean) * w + mean
    return np.moveaxis(x, 0, axis)


def smth9(x, p=0.5, q=0.25):
    """NCL smth9(x, p, q, False): 9-point smoother of the two rightmost dims, edges unchanged."""
    x = np.asarray(x, dtype=np.float64)
    out = x.copy()
    c = x[..., 1:-1, 1:-1]
    side = x[..., :-2, 1:-1] + x[..., 2:, 1:-1] + x[..., 1:-1, :-2] + x[..., 1:-1, 2:]
    corner = x[..., :-2, :-2] + x[..., :-2, 2:] + x[..., 2:, :-2] + x[..., 2:, 2:]
    out[..., 1:-1, 1:-1] = c + p / 4.0 * (side - 4.0 * c) + q / 4.0 * (corner - 4.0 * c)
    return out


# ----------------------------------------------------------------------------------------------
# missing-value fill, anomalies, area weights
# ----------------------------------------------------------------------------------------------


def fill_missing_time(da):
    """Fill NaN by linear interpolation in time at each grid point (NCL linmsg_n(x,-1,0)).

    Series ends take the nearest valid value. Used for obs with missing days (e.g. the 30
    all-missing days of the NOAA daily OLR CDR); the band-pass filter and EOFs cannot use them.
    """
    nmiss = int(da.isnull().sum())
    if nmiss == 0:
        return da
    tdim = da.dims[0]
    idx = xr.DataArray(np.arange(da.sizes[tdim]), dims=tdim)
    work = da.assign_coords({tdim: idx})
    work = work.interpolate_na(dim=tdim, method='linear')
    work = work.ffill(tdim).bfill(tdim)
    filled = nmiss - int(work.isnull().sum())
    print(f"MJO_suite: {da.name}: filled {filled} of {nmiss} missing values by linear interpolation in time")
    return work.assign_coords({tdim: da[tdim]})


def lat_weights(lat):
    """NCL MJO_suite utils.ncl get_gw (no gw on file): sin(lat) differences between midpoints."""
    lat = np.asarray(lat, dtype=np.float64)
    n = lat.size
    edges = np.empty(n + 1)
    if lat[0] < lat[1]:
        edges[0], edges[n] = -90.0, 90.0
    else:
        edges[0], edges[n] = 90.0, -90.0
    edges[1:n] = 0.5 * (lat[:-1] + lat[1:])
    return np.abs(np.sin(np.deg2rad(edges[1:])) - np.sin(np.deg2rad(edges[:-1])))


def daily_anomaly(da):
    """NCL calc_utils.ncl compute_daily_anom: subtract the raw (unsmoothed) daily climatology.

    The climatology is the mean over all years of each calendar day (month, day). Feb 29 (real
    calendars only) is grouped with Feb 28, matching NCL, which maps all data to a 365-day year.
    """
    tdim = da.dims[0]
    dates = _as_dates(da[tdim].values)
    mon = np.array([d.month for d in dates])
    day = np.array([d.day for d in dates])
    day = np.where((mon == 2) & (day == 29), 28, day)
    key = xr.DataArray(mon * 100 + day, dims=tdim, coords={tdim: da[tdim]}, name='mmdd')
    clim = da.groupby(key).mean(tdim)
    return da.groupby(key) - clim


# ----------------------------------------------------------------------------------------------
# EOFs
# ----------------------------------------------------------------------------------------------


def eofunc(x, neof):
    """NCL eofunc(x, neof, False) for x(space, time) (time is the RIGHTMOST axis here, as in NCL).

    Covariance-matrix EOFs. Time steps that are NaN at every point (e.g. the ends left by the
    band-pass filter) are dropped. Returns (eof(neof, space), pcvar(%), eval). As in NCL each EOF
    is a unit-length eigenvector; eval are the covariance-matrix eigenvalues (divisor nt-1).
    """
    x = np.asarray(x, dtype=np.float64)
    good_t = ~np.all(np.isnan(x), axis=0)
    xs = x[:, good_t]
    good_s = ~np.any(np.isnan(xs), axis=1)
    a = xs[good_s]
    a = a - a.mean(axis=1, keepdims=True)
    nt = a.shape[1]
    # SVD of the (space, time) anomaly matrix; covariance = a a^T / (nt - 1)
    u, s, _ = np.linalg.svd(a, full_matrices=False)
    lam = s ** 2 / (nt - 1)
    trace = lam.sum()
    eof = np.full((neof, x.shape[0]), np.nan)
    for n in range(neof):
        v = u[:, n]
        # NCL/LAPACK sign is arbitrary; fix it so the sum of the EOF is positive (reproducible)
        if v.sum() < 0:
            v = -v
        eof[n, good_s] = v
    pcvar = 100.0 * lam[:neof] / trace
    return eof, pcvar, lam[:neof]


def eofunc_ts(x, eof):
    """NCL eofunc_ts(x, eof, False) for x(space, time): project the time anomalies onto the EOFs.

    ts(n, t) = sum_s (x(s, t) - xmean(s)) * eof(n, s). Times with any NaN give NaN.
    """
    x = np.asarray(x, dtype=np.float64)
    good_t = ~np.all(np.isnan(x), axis=0)
    xmean = np.nanmean(x[:, good_t], axis=1, keepdims=True)
    e = np.where(np.isnan(eof), 0.0, eof)
    ts = e @ np.nan_to_num(x - xmean)
    ts[:, np.any(np.isnan(x[~np.all(np.isnan(eof), axis=0)]), axis=0)] = np.nan
    return ts


# ----------------------------------------------------------------------------------------------
# correlations
# ----------------------------------------------------------------------------------------------


def esccr(x, y, mxlag):
    """NCL esccr(x, y, mxlag) along the last axis: cross-correlation, y lagging x by 0..mxlag.

    r(k) = sum_{t=0}^{n-1-k} (x_t - xm)(y_{t+k} - ym) / ((n-k) * sx * sy), with population (1/n)
    standard deviations. NaN pairs are ignored (counts are of valid points).
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    x, y = np.broadcast_arrays(x, y)
    n = x.shape[-1]
    xm = np.nanmean(x, axis=-1, keepdims=True)
    ym = np.nanmean(y, axis=-1, keepdims=True)
    xa = x - xm
    ya = y - ym
    nx = np.sum(~np.isnan(xa), axis=-1)
    sx = np.sqrt(np.nansum(xa ** 2, axis=-1) / nx)
    sy = np.sqrt(np.nansum(ya ** 2, axis=-1) / np.sum(~np.isnan(ya), axis=-1))
    out = np.empty(x.shape[:-1] + (mxlag + 1,))
    for k in range(mxlag + 1):
        prod = xa[..., :n - k] * ya[..., k:]
        nk = np.sum(~np.isnan(prod), axis=-1)
        out[..., k] = np.nansum(prod, axis=-1) / nk / (sx * sy)
    return out


def escorc(x, y):
    """NCL escorc along the last axis: Pearson correlation."""
    return esccr(x, y, 0)[..., 0]


def season_starts(dates, season):
    """Start indices and length for NCL's MJO CLIVAR seasons (winter Nov 1 / summer May 1, 180 days;
    annual Jan 1, 365 days)."""
    mmdd = np.asarray(dates) % 10000
    if season == 'winter':
        start, nday = 1101, 180
    elif season == 'summer':
        start, nday = 501, 180
    elif season == 'annual':
        start, nday = 101, 365
    else:
        raise ValueError(f"unknown season {season}")
    return np.nonzero(mmdd == start)[0], nday


def mjo_xcor_lag(xio, y, dates, mxlag, season):
    """diagnostics_cam.ncl mjo_xcor_lag: lagged correlation of the base series xio(time) with
    y(time, N), averaged over the seasonal segments with Fisher's z. Returns ccr(2*mxlag+1, N)
    for lags -mxlag..mxlag (positive lag: y lags xio).

    Differences from NCL (both bugs there): segments are averaged over the number actually used
    (NCL divides by the number of season starts, including a final incomplete segment), and the
    last winter is dropped only when it is incomplete.
    """
    iSea, nday = season_starts(dates, season)
    ntim = len(xio)
    rz = 0.0
    nused = 0
    for i0 in iSea:
        i1 = i0 + nday - 1
        if i1 > ntim - 1:
            break
        xs = np.asarray(xio[i0:i1 + 1], dtype=np.float64)
        ys = np.asarray(y[i0:i1 + 1], dtype=np.float64).T          # (N, nday)
        y_lead_x = esccr(xs[np.newaxis, :], ys, mxlag)             # (N, mxlag+1)
        x_lead_y = esccr(ys, xs[np.newaxis, :], mxlag)
        ccr = np.empty((ys.shape[0], 2 * mxlag + 1))
        ccr[:, :mxlag] = x_lead_y[:, mxlag:0:-1]
        ccr[:, mxlag:] = y_lead_x
        ccr = np.clip(ccr, -1.0, 1.0)
        ccr = np.minimum(ccr, 0.999)
        rz = rz + 0.5 * np.log((1 + ccr) / (1 - ccr))
        nused += 1
    if nused == 0:
        raise ValueError(f"mjo_xcor_lag: no complete {season} season in the data")
    rz = rz / nused
    return np.tanh(rz).T                                            # (lag, N)


# ----------------------------------------------------------------------------------------------
# wavenumber-frequency spectra
# ----------------------------------------------------------------------------------------------


def resolve_waves_hayashi(work):
    """diagnostics_cam.ncl resolveWavesHayashi: rearrange the 2-D FFT coefficients work(mlon, N)
    (complex) into power pee(mlon+1, N+1) with wavenumber -mlon/2..mlon/2 and frequency
    -N/2..N/2 so that eastward waves have +wave/+freq (Hayashi 1971)."""
    mlon, n = work.shape
    m2, n2 = mlon // 2, n // 2
    vst = np.empty((mlon + 1, n + 1), dtype=complex)
    vst[:m2, :n2] = work[m2:0:-1, n2:n]
    vst[:m2, n2:] = work[m2:0:-1, :n2 + 1]
    vst[m2:, :n2 + 1] = work[:m2 + 1, n2::-1]
    vst[m2:, n2 + 1:] = work[:m2 + 1, n - 1:n2 - 1:-1]
    return np.abs(vst) ** 2


def mjo_wavenum_freq_season(x, dates, season, pct=0.10):
    """diagnostics_cam.ncl mjo_wavenum_freq_season for x(time, lon): pooled space-time power spectrum
    over the seasonal segments. Returns (power(wave, freq), wave, freq)."""
    x = np.asarray(x, dtype=np.float64)
    if np.isnan(x).any():
        raise ValueError("mjo_wavenum_freq_season: missing data not allowed")
    ntim, mlon = x.shape
    iSea, nday = season_starts(dates, season)
    x = dtrend(x, axis=0)
    pow_ = np.zeros((mlon + 1, nday + 1))
    k = 0
    for i0 in iSea:
        i1 = i0 + nday - 1
        if i1 > ntim - 1:
            break
        xs = x[i0:i1 + 1, :]
        xs = xs - xs.mean()
        xs = taper(xs, pct, axis=0)
        # FFT in lon (normalized by mlon), then in time (normalized by nday); NCL cfftf = numpy fft
        work = np.fft.fft(xs, axis=1) / mlon
        work = np.fft.fft(work, axis=0) / nday                     # (time-freq, wave)
        pow_ += resolve_waves_hayashi(work.T)
        k += 1
    if k == 0:
        raise ValueError(f"mjo_wavenum_freq_season: no complete {season} season in the data")
    pow_ /= k
    wave = np.arange(-(mlon // 2), mlon // 2 + 1)
    freq = np.linspace(-nday / 2, nday / 2, nday + 1) / nday
    return pow_, wave, freq
