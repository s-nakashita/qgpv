module read_qgpv
  implicit none
  include '/opt/local/include/netcdf.inc'

  integer,parameter :: ntime=7,imax=480,jmax=241,kmax=19!37
  real,save :: rlon(imax),rlat(jmax),rlev(kmax),rtime(ntime) 
  integer,save :: lev(kmax),time(ntime)  
contains

! for theta, dthdp
  subroutine fread_t(fname,vname,ip,z)
    implicit none
    integer :: istart(2),icount(2)
    integer :: ic,ncid,idvar,idlev,idtime
    integer,intent(in) :: ip
    real,intent(out) :: z(kmax)
    character,intent(in) :: fname*100,vname*5

    ic=NF_OPEN(fname,0,ncid)

    ic=NF_INQ_VARID(ncid,vname,idvar)
    ic=NF_INQ_VARID(ncid,'level',idlev)
    ic=NF_INQ_VARID(ncid,'time',idtime)
    print *, idvar
    
    ic=NF_GET_VAR_REAL(ncid,idlev,rlev)
    ic=NF_GET_VAR_REAL(ncid,idtime,rtime)

    istart(1) = 1
    istart(2) = ip

    icount(1) = kmax
    icount(2) = 1

    ic=NF_GET_VARA_REAL(ncid,idvar,istart,icount,z)

    ic=NF_CLOSE(ncid)

    return
  end subroutine fread_t

! for qgpv
  subroutine fread_q(fname,vname,ip,z)
    implicit none
    integer,allocatable :: istart(:),icount(:)
    integer :: ic,ncid,idvar,idlon,idlat,idlev,idtime
    integer,intent(in) :: ip
    real,intent(out) :: z(imax,jmax,kmax)
    character,intent(in) :: fname*100,vname*4

    if(ip.ne.0) then
      allocate(istart(4))
      allocate(icount(4))
    else
      allocate(istart(3))
      allocate(icount(3))
    endif
    ic=NF_OPEN(fname,0,ncid)

    ic=NF_INQ_VARID(ncid,vname,idvar)
    ic=NF_INQ_VARID(ncid,'lon',idlon)
    ic=NF_INQ_VARID(ncid,'lat',idlat)
    ic=NF_INQ_VARID(ncid,'level',idlev)
    if(ip.ne.0) ic=NF_INQ_VARID(ncid,'time',idtime)

    ic=NF_GET_VAR_REAL(ncid,idlon,rlon)
    ic=NF_GET_VAR_REAL(ncid,idlat,rlat)
    ic=NF_GET_VAR_REAL(ncid,idlev,rlev)
    if(ip.ne.0) ic=NF_GET_VAR_REAL(ncid,idtime,rtime)

    istart(1) = 1
    istart(2) = 1
    istart(3) = 1
    if(ip.ne.0) istart(4) = ip

    icount(1) = imax
    icount(2) = jmax
    icount(3) = kmax
    if(ip.ne.0) icount(4) = 1

    ic=NF_GET_VARA_REAL(ncid,idvar,istart,icount,z)

    ic=NF_CLOSE(ncid)

    return
  end subroutine fread_q
! for T,phi
  subroutine fread_v(fname,vname,ip,z)
    implicit none
    integer :: istart(4),icount(4)
    integer :: ic,ncid,idvar,idlon,idlat,idlev,idtime
    integer,intent(in) :: ip
    real,intent(out) :: z(0:imax-1,jmax,kmax)
    character(len=*),intent(in) :: fname,vname

    ic=NF_OPEN(fname,0,ncid)
    
    ic=NF_INQ_VARID(ncid,vname,idvar)
    ic=NF_INQ_VARID(ncid,'lon',idlon)
    ic=NF_INQ_VARID(ncid,'lat',idlat)
    ic=NF_INQ_VARID(ncid,'level',idlev)
    ic=NF_INQ_VARID(ncid,'time',idtime)

    ic=NF_GET_VAR_REAL(ncid,idlon,rlon)
    ic=NF_GET_VAR_REAL(ncid,idlat,rlat)
    ic=NF_GET_VAR_REAL(ncid,idlev,lev)
    ic=NF_GET_VAR_INT(ncid,idtime,time)

    istart(1) = 1
    istart(2) = 1
    istart(3) = 1
    istart(4) = ip

    icount(1) = imax
    icount(2) = jmax
    icount(3) = kmax
    icount(4) = 1

    ic=NF_GET_VARA_REAL(ncid,idvar,istart,icount,z)

    ic=NF_CLOSE(ncid)

    return
  end subroutine fread_v
!
! calcurate time steps between 2 dates
!
  subroutine calc_steps(idate,edate,dt,nt)
    implicit none
    integer,intent(in) :: idate,edate !format:yyyymmddhh
    integer,intent(in) :: dt ! timestep
    integer,intent(out) :: nt
    integer :: diff,ny,nm,nd,nh
    integer :: iy,im,id,ih
    integer :: ey,em,ed,eh

    print *, dt
    diff = idate
    iy = diff/1000000
    diff = diff - iy*1000000
    im = diff/10000
    diff = diff - im*10000
    id = diff/100
    diff = diff - id*100
    ih = diff
    print*, ih,id,im,iy

    diff = edate
    ey = diff/1000000
    diff = diff - ey*1000000
    em = diff/10000
    diff = diff - em*10000
    ed = diff/100
    diff = diff - ed*100
    eh = diff
    print*, eh,ed,em,ey

    ny = ey - iy
    if(em.ge.im)then
      nm = em - im
    else
      nm = em+12 - im
      ny = ny - 1
    endif
    if(ed.ge.id)then
      nd = ed - id
    else
      if(em.eq.3)then
        if(mod(ey,4).ne.0)then
          nd = ed+28 - id
          nm = nm - 1
        else
          nd = ed+29 - id
          nm = nm -1
        endif
      elseif(em.eq.1.or.em.eq.2.or.em.eq.4.or.em.eq.6.or.em.eq.9.or.em.eq.11)then
        nd = ed+31 - id
        nm = nm - 1
      else
        nd = ed+30 - id
        nm = nm - 1
      endif
    endif
    if(eh.ge.ih)then
      nh = ceiling(real(eh - ih))/dt
    else
      nh = ceiling(real(eh+24 - ih))/dt
      nd = nd - 1
    endif
    print*, nh,nd,nm,ny
    nt = nh + nd*(24/dt) + nm*(24/dt)*30 + ny*(24/dt)*365
    return
  end subroutine calc_steps
!
! calcurate date after any timesteps
!
  subroutine calc_date(idate,nt,edate)
    implicit none
    integer,intent(in) :: idate, nt !format:yyyymmddhh
    integer,intent(out) :: edate
    integer :: dt=6 !default timestep(6h)
    integer :: ey,em,ed,eh

    edate = idate + dt*nt
    eh = mod(edate,100)
    if (eh==24) then
      eh = 0
      edate = edate + 100 - 24
    endif
    ed = mod(edate,10000)/100
    em = mod(edate,1000000)/10000
    ey = (edate - em*10000 - ed*100 - eh)/1000000
    if (ed.gt.31.and.&
    &(em.eq.1.or.em.eq.3.or.em.eq.5.or.em.eq.7&
    &.or.em.eq.8.or.em.eq.10.or.em.eq.12))then
      edate = edate + 10000 - 3100
    elseif (ed.gt.30.and.&
      &(em.eq.4.or.em.eq.6.or.em.eq.9.or.em.eq.11))then
      edate = edate + 10000 - 3000
    elseif (em==2.and.ed.gt.28) then
      if (mod(ey,4)==0.and.ed.gt.29) then
        edate = edate + 10000 - 2900
      elseif(mod(ey,4).ne.0) then
        edate = edate + 10000 - 2800
      endif
    endif
    em = mod(edate,1000000)/10000
    if (em.gt.12) then
      edate = edate + 1000000 - 120000
    endif

    return
  end subroutine calc_date
!
! calcurate specific humidity from relative humidity
!
  SUBROUTINE calc_q(t,rh,p,q)
    IMPLICIT NONE
    double precision,PARAMETER :: t0=273.15d0
    double precision,PARAMETER :: e0c=6.11d0
    double precision,PARAMETER :: al=17.3d0
    double precision,PARAMETER :: bl=237.3d0
    double precision,PARAMETER :: e0i=6.1121d0
    double precision,PARAMETER :: ai=22.587d0
    double precision,PARAMETER :: bi=273.86d0
    real,INTENT(IN) :: t,rh,p !rh=%,p=hPa
    real,INTENT(OUT) :: q
    real :: e,es,tc

    tc = t-t0
    IF(tc >= 0.0d0) THEN
      es = e0c * exp(al*tc/(bl+tc))
    ELSE IF(tc <= -15.d0) THEN
      es = e0i * exp(ai*tc/(bi+tc))
    ELSE
      es = e0c * exp(al*tc/(bl+tc)) * (15.0d0+tc)/15.0d0 &
        + e0i * exp(ai*tc/(bi+tc)) * (-tc) / 15.0d0
    END IF

    e = rh * 0.01 * es
  
    q = e * 0.622d0 / (p - e * 0.378d0)
   
    RETURN
  END SUBROUTINE calc_q
end module read_qgpv
