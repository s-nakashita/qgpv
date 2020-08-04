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
  res@gsnCenterString = yyyymmddhh+"_"+phgt+"hPa"
  res@gsnLeftString = "wind"
  res@gsnRightString = "z T"
  res@gsnDraw = False
  res@gsnFrame = False
  res@gsnSpreadColors = True
  res@gsnSpreadColorStart = 5
  res@gsnMaskLambertConformal = True
  res@mpProjection = "LambertConformal"
  res@mpMinLonF = lon0
  res@mpMaxLonF = lon1
  res@mpCenterLonF = 0.5*(lon0+lon1)
  res@mpMinLatF = lat0
  res@mpMaxLatF = lat1
  res@mpGridAndLimbOn = True
  res@mpDataBaseVersion = "Ncarg4_1"
  res@mpGeophysicalLineThicknessF = 2
  res@cnFillOn = True
  res@cnLinesOn = False
  if (dT.ne.0) then
    res@cnLevelSelectionMode = "ManualLevels"
    res@cnMinLevelValF = T0
    res@cnMaxLevelValF = T1
    res@cnLevelSpacingF =  dT
  end if

  resc = True
  resc@gsnLeftString = ""
  resc@gsnRightString = ""
  resc@gsnDraw = False
  resc@gsnFrame = False
  resc@cnLevelSpacingF = dz
  resc@cnLineThicknessF = 3
  resc@cnInfoLabelOn = False

  resv = True
  resv@gsnDraw = False
  resv@gsnFrame = False
  resv@gsnLeftString = ""
  resv@gsnRightString = ""
  resv@vcRefAnnoSide = "Right"
  resv@vcRefAnnoString2On = False
  resv@vcRefMagnitudeF = 20.0
  resv@vcRefLengthF = 0.045
  resv@vcRefAnnoOrthogonalPosF = -1.0
  resv@vcRefAnnoArrowLineColor = "black"
  resv@vcRefAnnoArrowUseVecColor = False
  resv@vcGlyphStyle = "CurlyVector"
  resv@vcLineArrowColor = "black"
  resv@vcLineArrowThicknessF = 1
  resv@vcMinDistanceF = 0.02

  wks = gsn_open_wks(dev, "z-T-UV_"+phgt+"_"+yyyymmddhh+"_"+dset)
  gsn_define_colormap(wks,"posneg_1")

  z = f->z(t,{phgt},:,:);short2flt(f->z(t,{phgt},:,:))
  z = z / 9.81
  printMinMax(z, 0)
  T = f->t(t,{phgt},:,:);short2flt(f->t(t,{phgt},:,:))
  u = f->u(t,{phgt},:,:);short2flt(f->u(t,{phgt},:,:))
  v = f->v(t,{phgt},:,:);short2flt(f->v(t,{phgt},:,:))
  plot = gsn_csm_contour_map(wks, T, res)
  plotc = gsn_csm_contour(wks, z, resc)
  plotv = gsn_csm_vector(wks, u, v, resv)
  overlay(plot, plotc)
  overlay(plot, plotv)
  draw(plot)
  frame(wks)
end
