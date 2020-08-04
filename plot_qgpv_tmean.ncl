load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "qgpv.ncl"
load "config.txt"

begin
  f = addfile("bc.nc","r")
  ;f = addfile("era_201807_all.nc","r")
  time = f->time
  utime = cd_calendar(time, -3)
  nt = dimsizes(time)
  vt = closest_val(yyyymmddhh, utime)
  mmdd=str_get_cols(tostring(yyyymmddhh), 4, 7)

  ofile = "qgpv_mean.nc"
  if (isfilepresent(ofile)) then
    system("rm " + ofile)
  end if
  out = addfile(ofile,"c") 

  u = f->u;short2flt(f->u)
  v = f->v;short2flt(f->v)
  T = f->t;short2flt(f->t)
  lat = f->lat
  lon = f->lon 
  plev = f->level
  print(plev)
  latc = 0.5*(lat0+lat1)
  do t = 2, 6
    ;printMinMax(T(t,:,:,:), 0)
    Qlist = qgpv(u(t,:,:,:),v(t,:,:,:),T(t,:,:,:),latc,lat,lon,plev)
    Q = Qlist[0]
    s0 = Qlist[1]
    theta = Qlist[2]
    if (t.eq.2) then
      Q_ = Q
      s_ = s0
      th_ = theta
    else
      Q_ = Q_ + Q
      s_ = s_ + s0
      th_ = th_ + theta
    end if
    delete(Qlist)
    delete(Q)
    delete(s0)
    delete(theta)
  end do
  Q_ = Q_ / 5.
  s_ = s_ / 5.
  th_ = th_ / 5.
  printMinMax(Q_, 0)
;************************************************
; create netcdf4 file
;************************************************
  ndim = dimsizes(Q_)
  nlev = ndim(0)
  nlat = ndim(1)
  nlon = ndim(2)
  print(nlev)
  print(nlat)
  print(nlon)
  dim_names = (/"level","lat","lon"/)
  dim_sizes = (/nlev,nlat,nlon/)
  dim_unlimited = (/False,False,False/)
  filedimdef(out,dim_names,dim_sizes,dim_unlimited)

  atts = True
  atts@description = "quasi-geostrophic potential vorticity"
  atts@timestamp = systemfunc("date")
  print("Defining file attributes")
  fileattdef(out,atts)

  print("Defining file variables")
  filevardef(out,"level","float","level")
  filevardef(out,"lat","float","lat")
  filevardef(out,"lon","float","lon")
  filevardef(out,"qgpv","float",(/"level","lat","lon"/))
  filevardef(out,"theta","float","level")
  filevardef(out,"dthdp","float","level")
  
  print("Defining file variable attributes")
  levatts = True
  levatts@long_name = plev@long_name
  levatts@units = "hPa"
  filevarattdef(out,"level",levatts)

  latatts = True
  latatts@long_name = lat@long_name
  latatts@units = lat@units
  filevarattdef(out,"lat",latatts)

  lonatts = True
  lonatts@long_name = lon@long_name
  lonatts@units = lon@units
  filevarattdef(out,"lon",lonatts)

  vatts = 0.0
  vatts@long_name = "5day mean quasi-geostrophic potential vorticity"
  vatts@units     = "/s"
  vatts@_FillValue = 1e32
  filevarattdef(out,"qgpv",vatts)

  thatts = 0.0
  thatts@long_name = "reference potential temparature"
  thatts@units     = "K"
  thatts@_FillValue = 1e32
  filevarattdef(out,"theta",thatts)

  dtatts = 0.0
  dtatts@long_name = "reference potential temparature p-gradient"
  dtatts@units     = "K/Pa"
  dtatts@_FillValue = 1e32
  filevarattdef(out,"dthdp",dtatts)

  level = plev
  out->level = (/level/)
  out->lat   = (/lat/)
  out->lon   = (/lon/)

  out->qgpv = (/Q_(:,:,:)/)
  out->theta = (/th_(:)/)
  out->dthdp = (/s_(:)/)
  
;************************************************
; create plot
;************************************************
  res = True
  res@gsnCenterString = "mean_0703-0707"
  res@gsnLeftString = ""
  res@gsnRightString = ""
  res@gsnDraw = False
  res@gsnFrame = False
  res@gsnSpreadColors = True
  res@gsnSpreadColorStart = 5
  res@gsnMaskLambertConformal = True
  res@gsnAddCyclic = False
  res@mpProjection = "LambertConformal"
  res@mpMinLonF = lon0
  res@mpMaxLonF = lon1
  res@mpCenterLonF = 0.5*(lon0+lon1)
  res@mpMinLatF = lat0
  res@mpMaxLatF = lat1
  res@mpGridAndLimbOn = True
  res@mpDataBaseVersion = "Ncarg4_1"
  res@cnFillOn = True
  res@cnLinesOn = False
  res@cnLevelSelectionMode = "ManualLevels"
  
  wks = gsn_open_wks(dev, "qgpv"+phgt+"_201807"+dset)
  gsn_define_colormap(wks,"posneg_1")
  plot = new(2, graphic)

  Qlist = qgpv(u(vt,:,:,:),v(vt,:,:,:),T(vt,:,:,:),latc,lat,lon,plev)
  Q = Qlist[0]
  Q = Q - Q_
  s0 = Qlist[1]
  printMinMax(Q, 0)
  print(ind(plev.eq.phgt))
  nlev = dimsizes(plev)
  scale = -s0(ind(plev.eq.phgt))*9.81; -g*dth/dp
  print(scale)
  res@cnMinLevelValF = Q0/scale
  res@cnMaxLevelValF =  Q1/scale
  res@cnLevelSpacingF = dQ/scale
  plot(0) = gsn_csm_contour_map(wks, Q_({phgt},:,:), res)

  res@gsnCenterString = "anomaly_"+mmdd
  plot(1) = gsn_csm_contour_map(wks, Q({phgt},:,:), res)
  
  resP = True
  resP@gsnPanelMainString = dset+"_201807_"+phgt+"hPa"
  gsn_panel(wks, plot, (/1,2/), resP)
end
