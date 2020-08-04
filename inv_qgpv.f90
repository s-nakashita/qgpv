program inv_qgpv
    use read_qgpv
    use sor_m
    implicit none
    integer,parameter :: slon=161,elon=221,slat=141,elat=201
    integer,parameter :: nlon=elon-slon+1,nlat=elat-slat+1
    real :: dx(nlat),dy,dp(kmax),lat(nlat),lat0
    real :: Np(kmax),exn(kmax),N
    real :: phi(imax,jmax,kmax),phi_(imax,jmax,kmax)
    real :: q(imax,jmax,kmax),q_(imax,jmax,kmax)
    real :: theta(kmax),dthdp(kmax)
    real :: qpart(nlon,nlat,kmax),phip(nlon,nlat,kmax)
    real :: psi(nlon,nlat,kmax)
    real :: buf4(nlon,nlat)
    character*100 :: qfile='./qgpv.nc',qmfile='./qgpv_mean.nc',vfile='./bc.nc',ofile='./psi.grd'
    integer :: j,k,irec,it

    call fread_t(qfile,"theta",5,theta)
    call fread_t(qfile,"dthdp",5,dthdp)
    call fread_q(qfile,"qgpv",5,q)
    call fread_q(qmfile,"qgpv",0,q_)
    !print *, "lon",rlon
    !print *, "lat",rlat
    !print *, "level",rlev
    !print *, "time",rtime
    !print *, "theta",theta
    !print *, "dthdp",dthdp
    !print *, "qgpv",q(1,1,:)
    !print *, "qgpv(mean)",q_(1,1,:)

    qpart = q(slon:elon,slat:elat,:) - q_(slon:elon,slat:elat,:)
    
    phi_ = 0.0
    do it = 3, 7 
        call fread_v(vfile,"z",it,phi)
        phi_ = phi_ + phi
    end do
    phi_ = phi_/5.0
    call fread_v(vfile,"z",5,phi)
    phip = phi(slon:elon,slat:elat,:) - phi_(slon:elon,slat:elat,:)
    print *, "phi",minval(phip),maxval(phip)
    
    lat = rlat(slat:elat)
    print *, "lat", lat
    dy = 2*pi*r_earth*0.75/360.0
    do j = 1,nlat
        dx(j) = 2*pi*r_earth*cos(lat(j)*deg2rad)*0.75/360.0
    end do
    do k = 1,kmax-1
        dp(k) = (rlev(k+1)-rlev(k))*100.0 !Pa
    end do
    print *, dp
    lat0 = 45.0
!
! p-coordinate static_stability
!    
    do k = 1,kmax
        exn(k) = cp*(rlev(k)/p0)**kp
    end do
    print *, "kp",kp
    print *, "exner",exn
    do k = 1,kmax
        Np(k) = -rlev(k)*100*dthdp(k)/(kp*exn(k)*theta(k)**2)
    end do
    print *, "Np",Np
    !N = 2.0e-6
    !print *, "N",N
    call SOR(kmax,nlat,nlon,dx,dy,dp(1),lat0,lat,Np,qpart,phip,psi)
    !call SOR_s(kmax,nlat,nlon,dy,dy,dp(1),N,lat0,qpart,psi)
    print *, "psi",minval(psi),maxval(psi)
    open(11,file=ofile,status='new',&
    &    access='direct',convert='big_endian',&
    &    form='unformatted',recl=4*nlon*nlat)
    irec=1
    do k = 1, kmax
        buf4 = psi(:,:,k)
        WRITE(11,rec=irec) buf4
        irec = irec + 1
    end do
    close(11)

end program inv_qgpv