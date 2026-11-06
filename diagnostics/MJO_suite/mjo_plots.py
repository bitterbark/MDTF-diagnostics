# This file is part of the MJO_suite module of the MDTF code package (see LICENSE.txt)
"""Shared matplotlib helpers for the Python MJO_suite (NCL color tables, map axes, file names)."""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.colors as mcolors  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker as mticker  # noqa: E402
import numpy as np  # noqa: E402

import mjo_utils as U  # noqa: E402

# Anchor colors (every 10th entry) of the NCL color tables used by the NCL version of the POD.
_BLWHRE = [(0, 0, 128), (0, 0, 179), (0, 0, 230), (51, 51, 255), (153, 153, 255), (255, 255, 255),
           (255, 153, 153), (255, 51, 51), (230, 0, 0), (179, 0, 0), (128, 0, 0)]
_VIBLGRWHYEORRE = [(87, 0, 136), (29, 0, 215), (3, 56, 180), (8, 168, 30), (156, 225, 159),
                   (255, 255, 255), (255, 255, 153), (255, 249, 0), (255, 193, 0), (255, 110, 0),
                   (255, 0, 0)]
_PRECIP2_17LEV = [(255, 255, 255), (235, 246, 255), (214, 226, 255), (181, 201, 255), (142, 178, 255),
                  (127, 150, 255), (114, 133, 248), (99, 112, 248), (0, 158, 30), (60, 188, 61),
                  (179, 209, 110), (185, 249, 110), (255, 249, 19), (255, 163, 9), (229, 0, 0),
                  (189, 0, 0), (129, 0, 0), (0, 0, 0)]


def _smooth_cmap(name, rgb):
    return mcolors.LinearSegmentedColormap.from_list(name, [tuple(c / 255. for c in x) for x in rgb])


CMAP_BLWHRE = _smooth_cmap("BlWhRe", _BLWHRE)
CMAP_VIBLGRWHYEORRE = _smooth_cmap("ViBlGrWhYeOrRe", _VIBLGRWHYEORRE)


def precip_cmap(levels):
    """precip2_17lev spread over len(levels)+1 colors (NCL gsnSpreadColors); returns (cmap, norm)
    for contourf(..., levels=levels, extend="both")."""
    cols = np.array(_PRECIP2_17LEV) / 255.
    idx = np.round(np.linspace(0, len(cols) - 1, len(levels) + 1)).astype(int)
    cmap = mcolors.ListedColormap(cols[idx[1:-1]])
    cmap.set_under(cols[idx[0]])
    cmap.set_over(cols[idx[-1]])
    return cmap, mcolors.BoundaryNorm(levels, cmap.N)


def fig_path(name):
    """PNG path in $WORK_DIR/model/ (the names MJO_suite.html links to)."""
    d = os.path.join(U.env("WORK_DIR"), "model")
    os.makedirs(d, exist_ok=True)
    return os.path.join(d, name + ".png")


def save(fig, name):
    path = fig_path(name)
    fig.savefig(path, dpi=110, bbox_inches="tight")
    plt.close(fig)
    print(f"MJO_suite wrote {path}")
    return path


def sym_levels(data, nlev=8):
    """Symmetric 'nice' contour levels about zero (NCL symMinMaxPlt)."""
    mx = float(np.nanmax(np.abs(data)))
    if not np.isfinite(mx) or mx == 0:
        mx = 1.0
    return mticker.MaxNLocator(nlev, symmetric=True).tick_values(-mx, mx)


def lon_formatter(x, pos=None):
    x = ((x + 180) % 360) - 180
    if np.isclose(x, 0) or np.isclose(abs(x), 180):
        return f"{abs(int(round(x)))}"
    return f"{abs(int(round(x)))}{'E' if x > 0 else 'W'}"


def lat_formatter(y, pos=None):
    if np.isclose(y, 0):
        return "0"
    return f"{abs(int(round(y)))}{'N' if y > 0 else 'S'}"


def map_axes(fig, nrow, ncol, index, lat_s, lat_n, central_lon, lon_step=60, lat_step=15):
    import cartopy.crs as ccrs
    proj = ccrs.PlateCarree(central_longitude=central_lon)
    ax = fig.add_subplot(nrow, ncol, index, projection=proj)
    ax.set_extent([0, 359.99, lat_s, lat_n], crs=ccrs.PlateCarree())
    ax.coastlines(resolution="110m", linewidth=0.5)
    ax.set_xticks(np.arange(-180, 180, lon_step) + (central_lon % lon_step), crs=ccrs.PlateCarree())
    ax.set_yticks(np.arange(np.ceil(lat_s / lat_step) * lat_step, lat_n + 0.1, lat_step),
                  crs=ccrs.PlateCarree())
    from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter
    ax.xaxis.set_major_formatter(LongitudeFormatter(degree_symbol=""))
    ax.yaxis.set_major_formatter(LatitudeFormatter(degree_symbol=""))
    ax.tick_params(labelsize=7)
    return ax


def cyclic(data, lon):
    """Append the first longitude at +360 so filled maps have no gap at the dateline seam."""
    from cartopy.util import add_cyclic_point
    return add_cyclic_point(data, coord=lon)
