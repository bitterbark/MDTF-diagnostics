; This file is part of the MJO_suite module of the MDTF code package (see LICENSE.txt)

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
print("...L7")
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
print("...L9")
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
print("...L11")
load "$POD_HOME/calc_utils.ncl"
print("...L10")

begin
print("...L13")

nc_format = getenv("MDTF_NC_FORMAT")
print("...L16")
if ((nc_format.eq."NetCDF4_Classic").or.(nc_format.eq."NETCDF4_CLASSIC")) then
    nc_format = "NetCDF4Classic"
    print("...L19")
end if
print("...L21")
setfileoption("nc", "Format", nc_format)
print("...L23")

routine_name = "daily_anom"
casename = getenv("CASENAME")
time_coord = getenv("time_coord")
lat_coord = getenv("lat_coord")
lon_coord = getenv("lon_coord")

print("...L23")
; Use a path relative to WK_DIR (the process cwd, set by MDTF's
; SubprocessRuntimeManager) instead of the full absolute path. NCL 6.6.2's
; isfilepresent() copies its argument into a fixed 256-byte internal buffer
; with no length check and aborts the whole process (glibc catches it as
; "buffer overflow detected") once the path hits 256 characters, which long
; CESM case names reach easily since the case name is repeated multiple
; times in MDTF's generated work-directory path.
file_dir = "model/netCDF/"
datadir = getenv("DATADIR")
yr1 = stringtointeger(getenv("FIRSTYR"))
yr2 = stringtointeger(getenv("LASTYR"))
start_date = yr1*10000+201
end_date = (yr2+1)*10000+101
print("...L30")
delete_existing = False  ;True = overwrite, False = don't overwrite
if ( delete_existing ) then
  print("WARNING: daily_anom.ncl deleting existing anom files")
  system("ls "+file_dir+"/*anom*.nc")
end if

var_names = (/ "pr","rlut","u200","u850","v200","v850" /)
print("...L38")
do i=0,dimsizes(var_names)-1
  var_name = var_names(i)
  print("... "+i+" "+var_name)
  print(routine_name+" starting variable "+var_name)

  file_in   = file_dir +"/"+casename+"."+var_name+".day.nc"
  file_out  = file_dir +"/"+casename+"."+var_name+".day.anom.nc"

  if( isfilepresent(file_out) .and. (.not.delete_existing)) then
    print(file_out+" already exists. To force recalculation, set daily_anom.ncl delete_existing = True")

  else

    if( isfilepresent(file_out) .and. (delete_existing)) then
      print(file_out+" already exists. Overwriting since delete_existing = "+delete_existing)
      system("rm "+file_out)
    end if
 
    if( isfilepresent(file_in) ) then
    	print("reading "+file_in)
	f = addfile(file_in,"r")
	time = f->$time_coord$
	date = cd_calendar(time,-2)

	gw = f->gw
	invar = f->$var_name$($time_coord$|:,$lat_coord$|:,$lon_coord$|:)

	var = compute_daily_anom(invar)
	var!0 = "time"
	var!1 = "lat"
	var!2 = "lon"
	var&time = time
	var&lat = f->$lat_coord$
	var&lon = f->$lon_coord$

        ; this was already done in daily_netcdf
        if (var_name.eq."pr") then
          var = (/var*1000.0*86400.0/) ; mm/day
          var@units = "mm/day"
        end if

	print("writing "+file_out)
	fo = addfile(file_out,"c")

	fo->gw = gw
	fo->date = date
	fo->$var_name$ = var

	delete(f)
	delete(fo)
	delete(invar)
	delete(var)
 
      else
        print("ERROR: missing file that should have been made by daily_netcdf.ncl "+file_in)
      end if ; if( isfilepresent(file_in)
    end if ; isfilepresent(wkdir+ ... anom ...
    
  end do
print("...L102")
end
