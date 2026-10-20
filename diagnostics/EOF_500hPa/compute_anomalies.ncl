; This file is part of the EOF_500hPa module of the MDTF code package (see LICENSE.txt)

; FILE: anomalies.ncl

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
begin

setfileoption("nc", "Format", getenv("MDTF_NC_FORMAT"))

wk_dir = getenv("WORK_DIR")
casename = getenv("CASENAME")
yr1 = getenv("startdate")
yr2 = getenv("enddate")
convention = getenv("convention")

; 1 for native CESM, 0 for GFDL & CMIP
if (str_match_bool_ic(convention, "CMIP") .or. str_match_bool_ic(convention, "CMIP"))
   date_offset=0
else
   date_offset=1
end if


use_hybrid_sigma = (getenv("USE_HYBRID_SIGMA").eq."1")

; bundy 2026-10-20: MDTF v4 passes startdate/enddate as full dates (YYYYMMDD...); use the years
firstyr = stringtointeger(str_get_cols(yr1,0,3))
lastyr = stringtointeger(str_get_cols(yr2,0,3))
;;date_offset = stringtointeger(date_off) ; 1 for native CESM, 0 for GFDL & CMIP
nyrs = lastyr-firstyr+1
nmons = nyrs*12

time_coord = getenv("time_coord")
lat_coord = getenv("lat_coord")
lon_coord = getenv("lon_coord")
lev_coord = getenv("lev_coord")
file_path = getenv("DATADIR")

; Z500 
if( fileexists(wk_dir+"/model/netCDF/"+casename+".Z500.ANOMS.nc") ) then

print(wk_dir+"/model/netCDF/"+casename+".Z500.ANOMS.nc exists!!")

else

; TSJ edit 2021-01-12
if (use_hybrid_sigma) then
   hyam_var = "hyam" ;getenv("hyam_var")
   hybm_var = "hybm" ;getenv("hybm_var")
   zg_var = getenv("zg_hybrid_sigma_var")
   ps_var = getenv("ps_var")
   f_zg = addfile(getenv("ZG_HYBRID_SIGMA_FILE"),"r")
   f_ps = addfile(getenv("PS_FILE"),"r")
else
   zg_var = getenv("zg_var")
   f_zg = addfile(getenv("ZG_FILE"),"r")
end if
; end TSJ edit

time = f_zg->$time_coord$
lat = f_zg->$lat_coord$
lon = f_zg->$lon_coord$
nlat = dimsizes(lat)
nlon = dimsizes(lon)

slat = new((/nlat+1/),"double")
gw = new((/nlat/),"double")
gw!0 = lat_coord
gw&$lat_coord$ = lat
slat(0) = -90.d0
slat(nlat) = 90.d0
do i=1,nlat-1
 slat(i) = (lat(i-1)+lat(i))/2.d0
end do
pi = 4.d0*atan(1.d0)
do i=0,nlat-1
 gw(i) = sin(slat(i+1)/180.d0*pi)-sin(slat(i)/180.d0*pi)
end do

time_all = f_zg->$time_coord$

; bundy 2026-10-20: label each sample with the month it averages (YYYYMM), whatever the
; time-stamp convention. Older CESM stamps a monthly mean at the end of the month (00:00 on the
; 1st of the next month), newer CESM and CMIP in the middle. Use the center of the time bounds
; when the file has them; otherwise take a stamp at 00:00 on the 1st as the end of the previous
; month. Then keep only the full years (January to December) inside the requested range.
if (isfilevar(f_zg,"time_bnds") .or. isfilevar(f_zg,"time_bounds")) then
   if (isfilevar(f_zg,"time_bnds")) then
      tb = f_zg->time_bnds
   else
      tb = f_zg->time_bounds
   end if
   tmid = time_all
   tmid = (/ 0.5*(tb(:,0)+tb(:,1)) /)
   utc = cd_calendar(tmid,0)
   yy = toint(utc(:,0))
   mm = toint(utc(:,1))
   print("Monthly samples labeled by the center of their time bounds")
   delete([/tb, tmid/])
else
   utc = cd_calendar(time_all,0)
   yy = toint(utc(:,0))
   mm = toint(utc(:,1))
   endstamp = toint(utc(:,2)).eq.1 .and. utc(:,3).eq.0 .and. utc(:,4).eq.0
   if (all(endstamp)) then
      print("Monthly samples stamped at the end of the month (no time bounds); using the previous month")
      yy = where(mm.eq.1, yy-1, yy)
      mm = where(mm.eq.1, 12, mm-1)
   else
      print("Monthly samples labeled by their time stamp's month (no time bounds)")
   end if
   delete(endstamp)
end if
date_all = yy*100+mm
delete([/utc, yy, mm/])
do while (firstyr.lt.lastyr .and. .not.any(date_all.eq.firstyr*100+1))
   firstyr = firstyr+1
end do
do while (lastyr.gt.firstyr .and. .not.any(date_all.eq.lastyr*100+12))
   lastyr = lastyr-1
end do
nyrs = lastyr-firstyr+1
nmons = nyrs*12
start_date = firstyr*100+1
end_date = lastyr*100+12
print("Using "+start_date+" to "+end_date)

; in either case, above only used to set following indices
do i=0,dimsizes(date_all)-1
   if( date_all(i).eq.start_date ) then
      i1 = i
   end if
   if( date_all(i).eq.end_date ) then
      i2 = i
   end if
end do

win_anom = new ((/nyrs-1,nlat,nlon/),"float")
win_anom!0 = "year"
win_anom&year = ispan(firstyr+1,lastyr,1) 
win_anom!1 = lat_coord
win_anom&$lat_coord$ = lat
win_anom!2 = lon_coord
win_anom&$lon_coord$ = lon
win_anom@long_name = "Winter (DJFM) 500mb height anomaly"
win_anom@units = "meters"

; TSJ edit 2021-01-12: temporary hack to extract 500 mb level of f_zg (z500)
if (use_hybrid_sigma) then
   ; hybrid sigma level extraction
   hyam = f_zg->$hyam_var$
   hybm = f_zg->$hybm_var$
   z = f_zg->$zg_var$($time_coord$|i1:i2,$lev_coord$|:,$lat_coord$|:,$lon_coord$|:)
   ps = f_ps->$ps_var$($time_coord$|i1:i2,$lat_coord$|:,$lon_coord$|:)
   z500_4D = vinth2p(z,hyam,hybm,500.,ps,2,1000.,1,False) ; (time,plev,lat,lon)
   ; level dimension is now degenerate (length 1), so drop it to get 3D data
   z500_dims = dimsizes(z500_4D)
   z500 = reshape(z500_4D, (/z500_dims(0), z500_dims(2), z500_dims(3)/))
   delete (z)
   delete (z500_4D)
   delete (ps)
else
   rank = dimsizes( filevardimsizes(f_zg,zg_var) )
   if (rank.ne.3 .or. rank.gt.4) then
      print("getVarSlice: only 3D supported: rank="+rank+"D")
      exit
   end if
   ; 3D data; assume we got 500 hPa level
    z500 = f_zg->$zg_var$($time_coord$|i1:i2,$lat_coord$|:,$lon_coord$|:)
end if
; end TSJ edit
anom = (/rmMonAnnCycTLL(z500(:,:,:))/)
delete (z500)
n = 11
do yr = firstyr+1,lastyr
  win_anom(yr-firstyr-1,:,:) = (/(anom(n,:,:)+anom(n+1,:,:)+anom(n+2,:,:)+anom(n+3,:,:))/4./)
  n = n + 12
end do

fo = addfile(wk_dir+"/model/netCDF/"+casename+".Z500.ANOMS.nc","c")
fo->Z500_ANOM = win_anom
fo->gw = gw

if (use_hybrid_sigma) then
   delete(f_ps)
end if
delete(f_zg)
delete(win_anom)

end if


end
