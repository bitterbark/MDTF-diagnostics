MJO Suite Diagnostic Module From NCAR
=====================================
Last update: 03/11/2019

This module computes many of the diagnostics described by the the ​US-CLIVAR Madden-Julian Oscillation (MJO) working group​ and developed by NCAR’s ​Dennis Shea for observational data​. Using daily precipitation, outgoing longwave radiation, zonal wind at 850 and 200 hPa and meridional wind at 200hPa, the module computes anomalies, bandpass-filters for the 20-100 day period, calculates the MJO Index as defined as the running variance over the bandpass filtered data, performs an EOF analysis, and calculates lag cross-correlations, wave-number frequency spectra and composite life cycles of MJO events.

Contact info
------------

- PI: Rich Neale, NCAR
- Current Developer: Dani Coleman (​bundy@ucar.edu​), NCAR
- Contributors: Dennis Shea, Andrew Gettleman, Jack Chen (NCAR)

Open source copyright agreement
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This package is distributed under the LGPLv3 license (see LICENSE.txt).

Functionality
-------------

The POD is written in Python (it was ported from NCL in 2026; the Python version reproduces the
NCL numbers, see the module docstrings). Scripts are in ``diagnostics/MJO_suite``; the driver
``MJO_suite.py`` runs, in order:

1. Read in daily model output: daily_netcdf.py (missing days are filled by linear interpolation in time)
2. Computes daily anomalies: daily_anom.py
3. Calculates EOFs: mjo_EOF.py
4. Creates MJO lag plots: mjo_lag_lat_lon.py
5. Calculates MJO spectra: mjo_spectra.py
6. Calculates principal components from the combined EOF analysis: mjo_EOF_cal.py
7. MJO life cycle composite: mjo_life_cycle.py
8. U200 band-pass time series: mjo_u200_index.py

Numerical routines (Lanczos filter, EOFs, lagged correlations, wavenumber-frequency spectra) are in
``mjo_utils.py``; plotting helpers are in ``mjo_plots.py``. Python ports of several of these
routines were started in https://github.com/cyschneck/geocat-research/tree/main/mjo .

| Preprocessed observational data from NCEP, ERA and TRMM are located in ``mdtf/inputdata/obs_data/MJO_suite``.
| Place your input data at: ``mdtf/inputdata/model/$model_name/day``
| index.html can be found at: ``mdtf/MDTF_$ver/wkdir/MDTF_$model_name``

Required Programing Language and libraries
------------------------------------------

| Python 3 with numpy, scipy, xarray, netCDF4, cftime, matplotlib and cartopy
  (all in the MDTF ``_MDTF_python3_base`` environment). NCL is no longer needed.

Required input data to the module
---------------------------------

The following five 3-D (lat-lon-time) model fields are required with a daily time output

1. Precipitation rate (units: mm/s = kg/m\ |^2|/s) or mm/day with appropriate conversion
2. Outgoing Longwave radiation (units: W/m\ |^2|)
3. Zonal wind at 850hPa (units: m/s)
4. Zonal wind at 200hPa (units: m/s)
5. Meridional wind at 200hPa (units:m/s)

References
----------

   .. _1:

1. MJO Simulation Diagnostics, Waliser et al, 2009, *J.Clim.*, **​22**:​ 3006-3030, https://doi.org/10.1175/2008JCLI2731.1

More About the Diagnostic
-------------------------

An extensive explanation of the figures and techniques used to achieve them can be found on Dennis Shea’s NCL MJO CLIVAR page​. Some examples are provided below:

.. image:: mjo_suite_fig1.png
   :align: center
   :width: 75 %

**Lag Correlation**: Lag correlation diagram using on 20-100 day band-pass filtered daily data. The reference time series is the central Indian Ocean regional precipitation time series, which is correlated with precipitation and zonal wind anomalies in specified regions at different lags. Lag-longitude and lag-latitude plots of correlation values for different regions are shown. Color is for precipitation correlations while the lagged correlations for the zonal winds are the contours. These are analogous the Figures 5 and 6 in the reference article except they are for one season.

.. image:: mjo_suite_fig2.png
   :align: center
   :width: 75 %

**Wavenumber-Frequency Spectra**: The wavenumber - frequency spectra for each season, with vertical reference lines indicating at for 30 and 80 days.


.. image:: mjo_suite_fig3.png
   :align: center
   :width: 100 %


**EOF analysis (univariate)**: Conventional (covariance) univariate EOF analysis for 20-100 day band-pass filtered data of outgoing longwave radiation over 30S to 30N.

.. image:: mjo_suite_fig4.png
   :align: center
   :width: 100 %

**Composite Life-cycles**: The first and second principle components of the EOF analysis are used to derive the appropriate MJO phase category. The size of the reference anomaly wind vector is in the upper right. The phase (eg P3, means "Phase 3") and the number of days used to create the composite are at the lower right.


.. |^2| replace:: \ :sup:`2`\ 
.. |^3| replace:: \ :sup:`3`\ 
.. |^-1| replace:: \ :sup:`-1`\ 
.. |^-2| replace:: \ :sup:`-2`\ 
.. |^-3| replace:: \ :sup:`-3`\ 