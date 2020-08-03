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

function qgpv(u:numeric, v:numeric, T:numeric, latc:numeric, 
              lat:numeric, lon:numeric, lev:numeric)
; input
;   u(:,:,:) (m/s), v(:,:,:) (m/s), T(:,:,:) (K), latc (deg), 
;   lat(:) (deg), lon(:) (deg), lev(:) (hPa)
; output
;   Qlist[0] qgpv
;   Qlist[1] d(theta_ref)/dp
;   Qlist[2] theta_ref
; NB. gradient and divergence is computed using SPHEREPACK functions
;     no missing values are allowed
local theta,dtheta,Tr,Ta,s1,s2,zeta,plev,dp
begin
  n = dimsizes(T)
  nlev = n(0)
  nlat = n(1)
  nlon = n(2)
  Qlist = NewList("lifo")
  ;print(n)
  
  f0   = 2.0 * Omega * sin(latc*deg2rad)

  plev = lev*100
  plev@units = "Pa"
  printMinMax(plev,0)
  psfc = max(plev)
  ptop = min(plev)
  dp = dpres_plevel_Wrap(plev, psfc, ptop, 0)
  ;print(dp)
; calcurate potential temperature
  printMinMax(T,0)
  theta = pot_temp(plev, T, 0, 1)
  printMinMax(theta, 0)
  ;do k = 0, nlev-1
  ;  do j = 0, nlat-1
  ;    do i = 0, nlon-1
  ;      theta(k,j,i) = T(k,j,i)*(plev(k)/p0)^(-kp)
  ;    end do
  ;  end do
  ;end do  
; calcurate static_stability
  Tr = dim_avg_n(T,(/1,2/))
  Ta = T
  Ta = Ta - conform_dims((/nlev,nlat,nlon/), Tr, 0);anomaly
  printMinMax(Ta,0)
  Tr@units = "K"
  s1 = static_stability(plev, Tr, 0, 1); s(0) = -(T/theta)*d(theta)/dp, 
                                         ; s(1) = theta_r, 
                                         ; s(2) = d(theta)/dp
  ;print(s1[1])
  ListPush(Qlist,s1[1])
  ListPush(Qlist,s1[2])
  theta = theta - conform_dims((/nlev,nlat,nlon/), s1[1], 0);anomaly
  theta = theta / conform_dims((/nlev,nlat,nlon/), s1[2], 0);theta'/(d(theta_r)/dp)
  dtheta = center_finite_diff_n(theta, conform(theta,plev,0), False, 0, 0)
  ;Ta = Ta / conform_dims((/nlev,nlat,nlon/), s1[2], 0);T'/(d(theta_r)/dp)
  ;printMinMax(Ta,0)
  ;s2 = static_stability(plev, Ta, 0, 1)
  ;dtheta = s2[2]
  printMinMax(dtheta,0)
; calcurate relative vorticity
  zeta = uv2dvF(v, -u)  ;global
  ;zeta = uv2dv_cfd(v, -u, lat, lon, 2) ;regional
  
; calcurate pv
  q = T
  q@long_name = "quasi geostrophic potential vorticity"
  q@units = "/s"
  do k = 0, nlev-1
    do j = 0, nlat-1  
      do i = 0, nlon-1
        q(k,j,i) = 2.0 * Omega * sin(lat(j)*deg2rad) + zeta(k,j,i) + f0 * dtheta(k,j,i)
      end do
    end do
  end do
  ListPush(Qlist, q)
  ;printVarSummary(q)
 ; Q(2,:,:) = smth9(Q(2,:,:),0.5,0.25,True)
  return Qlist
 end

