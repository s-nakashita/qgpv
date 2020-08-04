load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "qgpv.ncl"
load "config.txt"

begin
  f = addfile("bc.nc","r")
  ;f = addfile("era_201807_all.nc","r")
  f2 = addfile("era_201807.nc","r")
  time = f->time
  utime = cd_calendar(time, -3)
  vt = closest_val(yyyymmddhh, utime)

  ofile = "qgpv.nc"
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
  ;print(plev)
  latc = 0.5*(lat0+lat1)
  ndim = dimsizes(u)
  nlev = ndim(1)
  nlat = ndim(2)
  nlon = ndim(3)
  print(nlev)
  print(nlat)
  print(nlon)
  Q = new((/7,nlev,nlat,nlon/), typeof(u))
  s0 = new((/7,nlev/), typeof(u))
  theta = new((/7,nlev/), typeof(u))
  do t = 0, 6
    Qlist = qgpv(u(t,:,:,:),v(t,:,:,:),T(t,:,:,:),latc,lat,lon,plev)
    Q(t,:,:,:) = Qlist[0]
    s0(t,:) = Qlist[1]
    theta(t,:) = Qlist[2]
    printMinMax(Q, 0)  
  end do
;************************************************
; create netcdf4 file
;************************************************
  dim_names = (/"time","level","lat","lon"/)
  dim_sizes = (/7,nlev,nlat,nlon/)
  dim_unlimited = (/True,False,False,False/)
  filedimdef(out,dim_names,dim_sizes,dim_unlimited)

  atts = True
  atts@description = "quasi-geostrophic potential vorticity"
  atts@timestamp = systemfunc("date")
  print("Defining file attributes")
  fileattdef(out,atts)

  print("Defining file variables")
  filevardef(out,"time","float","time")
  filevardef(out,"level","float","level")
  filevardef(out,"lat","float","lat")
  filevardef(out,"lon","float","lon")
  filevardef(out,"qgpv","float",(/"time","level","lat","lon"/))
  filevardef(out,"theta","float",(/"time","level"/))
  filevardef(out,"dthdp","float",(/"time","level"/))
  
  print("Defining file variable attributes")
  tatts = True
  tatts@standard_name = time@standard_name
  tatts@units = time@units
  filevarattdef(out,"time",tatts)

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
  vatts@long_name = "quasi-geostrophic potential vorticity"
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

  t0 = time(0:6)
  out->time = (/t0/)
  out->level = (/plev/)
  out->lat   = (/lat/)
  out->lon   = (/lon/)

  out->qgpv = (/Q(:,:,:,:)/)
  out->theta = (/theta(:,:)/)
  out->dthdp = (/s0(:,:)/)
  
;************************************************
; create plot
;************************************************

  res = True
  res@gsnCenterString = "qgpv"
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
  
  wks = gsn_open_wks(dev, "qgpv"+phgt+"_"+yyyymmddhh+dset)
  gsn_define_colormap(wks,"posneg_1")
  plot = new(2, graphic)

  print(ind(plev.eq.phgt))
  nlev = dimsizes(plev)
  scale = -s0(vt,ind(plev.eq.phgt))*9.81; -g*dth/dp
  print(scale)
  res@cnMinLevelValF = Q0/0.005
  res@cnMaxLevelValF =  Q1/0.005
  res@cnLevelSpacingF = dQ/0.005
  plot(0) = gsn_csm_contour_map(wks, Q(vt,{phgt},:,:), res)

  pv = short2flt(f2->pv(vt,{phgt},:,:))
  res@gsnCenterString = "pv"
  res@cnMinLevelValF = Q0
  res@cnMaxLevelValF =  Q1
  res@cnLevelSpacingF = dQ
  plot(1) = gsn_csm_contour_map(wks, pv, res)
  
  resP = True
  resP@gsnPanelMainString = dset+"_"+yyyymmddhh+"_"+phgt+"hPa"
  gsn_panel(wks, plot, (/1,2/), resP)
end
