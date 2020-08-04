load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "config.txt"

begin
  f = addfile("era_201807_allf.nc","r")
  time = f->time
  utime = cd_calendar(time, -3)
  t = closest_val(yyyymmddhh, utime)

  res = True
  res@tiMainString = dset
  res@gsnCenterString = yyyymmddhh
  res@gsnLeftString = "slp"
  res@gsnRightString = "hPa"
  res@gsnMaskLambertConformal = True
  res@mpProjection = "LambertConformal"
  res@mpMinLonF = lon0
  res@mpMaxLonF = lon1
  res@mpCenterLonF = 0.5*(lon0+lon1)
  res@mpMinLatF = lat0
  res@mpMaxLatF = lat1
  res@mpGridAndLimbOn = True
  res@mpDataBaseVersion = "Ncarg4_1"
  res@cnLineThicknessF = 3
  if (dslp.ne.0) then
    res@cnLevelSelectionMode = "ManualLevels"
    res@cnMinLevelValF = slp0
    res@cnMaxLevelValF = slp1
    res@cnLevelSpacingF = dslp
  end if

  wks = gsn_open_wks(dev, "slp"+yyyymmddhh+dset)

  slp = f->slp(t,:,:);short2flt(f->msl(t,:,:))
  slp = slp * 0.01

  plot = gsn_csm_contour_map(wks, slp, res)
end
