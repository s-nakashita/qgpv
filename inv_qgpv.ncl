external SOR "./sor90.so"

pi = acos(-1.)
g_earth = 9.81; m/s^2
r_earth = 6.37*10^6; m
dayinsec = 86400.; s
Omega = 2 * pi / dayinsec
deg2rad = pi / 180.
Rd = 287.; J/kg/K
cp = 1004.; J/kg/K
kp = Rd / cp
p0 = 1000.*10^2; Pa

function inv_qgpv(q:numeric, phi:numeric, T:numeric, theta:numeric, dthdp:numeric, lat0:numeric, lat:numeric, lon:numeric, lev:numeric)
; input
;   q(:,:,:) (/s), phi(:,:,:) (m^2/s^2), T(:,:,:) (K), theta(:) (K), dthdp(:) (K), 
;   lat0 (deg), lat(:) (deg,S->N), lon(:) (deg), lev(:) (hPa,down)
; output
;   Qlist[0] u
;   Qlist[1] v
;   Qlist[2] z
; NB. gradient and divergence is computed using SPHEREPACK functions
;     no missing values are allowed
local plev,Np,exner,psi,fcori,grad
begin
  
  n = dimsizes(q)
  nlev = n(0)
  nlat = n(1)
  nlon = n(2)
  Qlist = NewList("lifo")
  ;print(n)
  ;lat = new(nlat,float)
  ;lat = (/lat1(::-1)/)
  ;dx = 2 * pi * r_earth * 0.75 / 360.
  ;dy = new(nlat, float)
  ;do j = 0, nlat-1
  ;  dy(j) = 2 * pi * r_earth * cos(lat(j)*deg2rad) * 0.75 / 360.
  ;end do
  plev = lev*100
  plev@units = "Pa"
  ;print(plev)
  ;dp = new(nlev,float)
  ;dp = 0.0
  ;do k = 1, nlev-2
  ;  dp(k) = plev(k+1)-plev(k)
  ;end do
;
; calcurate p-coordinate static_stability
;
  ;print(theta)
  ;print(dthdp)
  Np = new(nlev,float)
  exner = new(nlev,float)
  exner = cp*(plev/p0)^kp
  Np = -plev * dthdp / (kp * exner * theta^2)
;  print(dx)
;  print(dy)
;  print(dp) 
;  print(kp)
;  print(exner)
;  print(Np)
; print(lat)
;; invert qgpv
;  psi = new((/nlev,nlat,nlon/),float) 
;  SOR::SOR(nlev,nlat,nlon,lat,lat0,plev,Np,T,phi,q,psi) 
  f = addfile("psi.nc","r")
  psi = f->psi(0,:,:,:)
;  print(dimsizes(psi))
; calcurate variables
  z = q
  z@units = "m"
  v = q
  v@units = "m/s"
  u = q
  u@units = "m/s"
  ;print(z&latitude)

  fcori = lat
  fcori = 2 * Omega * sin(lat*deg2rad)
  z = conform_dims((/nlev,nlat,nlon/),fcori,1) * psi / g_earth 
  ;z = z(:,::-1,:)
  ListPush(Qlist,z)

  gradsf(psi, v, u) ;global
  u = -u
  ;grad = grad_latlon_cfd(z, lat, lon, False, False) ;regional
  ;v    = grad[0]
  ;u    = -grad[1] 
  ListPush(Qlist,v)
  ListPush(Qlist,u)
  return Qlist
 end

function p_stab(theta:numeric,plev:numeric)
; input
;   theta(:) (K,reference), plev(:) (Pa or hPa)
; output
;   Np(:) (m**4 s**2 kg**-2,p-coordinate static stability)
local nlev,p,exner,dthdp
begin
  nlev = dimsizes(theta)
  p = plev
  if(isatt(plev,"units").and.(plev@units.eq."hPa".or.plev@units.eq."millibars")) then
    p = plev*100
  end if
  Np = new(nlev,float)
  exner = new(nlev,float)
  exner = cp*(p/p0)^kp
  dthdp = new(nlev,float)
  dthdp = center_finite_diff_n(theta, p, False, 0, 0)
  Np = -p * dthdp / (kp * exner * theta^2)
  Np@standard_name = "p-coordinate static stability"
  Np@units         = "m**4 s**2 kg**-2"
  return Np
end