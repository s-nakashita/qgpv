load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "inv_qgpv.ncl"
load "config.txt"

begin
  ;f = addfile("era_201807_all.nc","r")
  f = addfile("qgpv.nc","r")
  fm = addfile("qgpv_mean.nc","r")
  g = addfile("psi.nc","r")
  time = f->time
  utime = cd_calendar(time, -3)
  t = closest_val(yyyymmddhh, utime)
  print(t)

  res = True
;  res@tiMainString = dset
  res@gsnLeftString = ""
  res@gsnRightString = ""
;  res@gsnMaximize = True
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
  res@cnMinLevelValF = -0.5
  res@cnMaxLevelValF = 0.5
  res@cnLevelSpacingF = 0.1

  resz = True
  resz@gsnDraw = False
  resz@gsnFrame = False
  resz@gsnLeftString = ""
  resz@gsnRightString = ""
  resz@cnFillOn = False
  resz@cnLinesOn = True
  resz@gsnContourNegLineDashPattern = 1
  resz@cnLevelSelectionMode = "ManualLevels"
  resz@cnMinLevelValF = z0
  resz@cnMaxLevelValF = z1
  resz@cnLevelSpacingF = dzb

  resv = True
  resv@gsnDraw = False
  resv@gsnFrame = False
  resv@gsnLeftString = ""
  resv@gsnRightString = ""
  resv@vcRefAnnoSide = "Right"
  resv@vcRefAnnoString2On = False
  resv@vcRefMagnitudeF = vref
  resv@vcRefLengthF = 0.045
  resv@vcRefAnnoOrthogonalPosF = -1.0
  resv@vcRefAnnoArrowLineColor = "black"
  resv@vcRefAnnoArrowUseVecColor = False
  resv@vcGlyphStyle = "CurlyVector"
  resv@vcLineArrowColor = "black"
  resv@vcLineArrowThicknessF = 3
  resv@vcMinDistanceF = 0.02

  wks = gsn_open_wks(dev, "fort_inv_qgpv"+phgt+"_"+yyyymmddhh+dset)
  gsn_define_colormap(wks,"posneg_1")
  plot = new(2, graphic)

  ;u = short2flt(f->u(:,:,:,:))
  ;v = short2flt(f->v(:,:,:,:))
  ;T = short2flt(f->t(:,:,:,:))
  ;phi = short2flt(f->z(:,:,:,:))
  ;z = phi
  ;z = phi / g_earth
  Q     = f->qgpv(t,:,:,:)
  dthdp = f->dthdp(t,:)
  theta = f->theta(t,:)
  ;printMinMax(Q, 0)
  Q_ = fm->qgpv
  Q = Q - Q_
  ;lat = f->lat
  ;lon = f->lon 
  plev = f->level
  nlev = dimsizes(plev)
  Np = p_stab(theta,plev)
  print(Np)
  ;resp = True
  ;resp@tiMainString = "static stability"
  ;resp@trYReverse   = True
  ;plot(1) = gsn_csm_xy(wks, Np, plev, resp)
;; ncl  
  ;T_ = T(2,:,:,:)
  ;phi_ = phi(2,:,:,:)
  ;do it=3,6
  ;  T_ = T_ + T(it,:,:,:)
  ;  phi_ = phi_ + phi(it,:,:,:)
  ;end do
  ;T_ = T_ / 5.
  ;phi_ = phi_ / 5.
  ;T_ = T(t,:,:,:)-T_
  ;phi_ = phi(t,:,:,:)-phi_
  
  ;latc = 0.5*(lat0+lat1)
  ;vlist = inv_qgpv(Q, phi_,T_,theta,dthdp,latc,lat,lon,plev)
  ;ub = vlist[0]
  ;vb = vlist[1]
  ;zb = vlist[2]
;; fortran
  psi = g->psi(0,:,:,:)
  printVarSummary(psi)
  phi = psi
  phi@units = "m**2 s**-2"
  zb = psi
  zb@units = "m"
  vb = psi
  vb@units = "m/s"
  ub = psi
  ub@units = "m/s"
  Tb = psi
  Tb@units = "K"
  ;print(z&latitude)

  lat = dble2flt(g->lat)
  lon = dble2flt(g->lon)
  fcori = lat
  fcori = 2 * Omega * sin(lat*deg2rad)
  phi = conform(psi,fcori,1) * psi
  zb = phi / g_earth
  Tb = center_finite_diff_n(phi, conform(phi,plev,0), False, 0, 0) 
  Tb = -Rd*Tb/conform(Tb,plev,0)
  ;gradsf(psi, vb, ub) ;global
  ;ub = -ub
  grad = grad_latlon_cfd(psi, lat, lon, False, False)
  ub    = -grad[0]
  vb    = grad[1]
  ilev = ind(plev.eq.phgt)
  res@gsnCenterString = "balance"
  plott = gsn_csm_contour_map(wks, Tb(ilev,:,:), res)
  plotz = gsn_csm_contour(wks, zb(ilev,:,:), resz)
  overlay(plott,plotz)
  plotv = gsn_csm_vector(wks, ub(ilev,:,:), vb(ilev,:,:), resv)
  overlay(plott, plotv)
  plot(0) = plott

  res@gsnCenterString = "qgpv"
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnMinLevelValF = Q0*10^2
  res@cnMaxLevelValF =  Q1*10^2
  res@cnLevelSpacingF = dQ*10^2
  plot(1) = gsn_csm_contour_map(wks, Q({phgt},:,:), res)
  
  resP = True
  resP@gsnPanelMainString = dset+"_"+yyyymmddhh+"_"+phgt+"hPa"
  gsn_panel(wks, plot, (/1,2/), resP)
end
