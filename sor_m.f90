module sor_m
    !$ USE OMP_LIB
    IMPLICIT NONE
    REAL,PARAMETER :: pi=atan(1.0)*4.0, omega=2*pi/86400.0, deg2rad=pi/180.0
    REAL,PARAMETER :: r_earth=6.37*10**6, Rd=287., cp=1004., kp=Rd/cp, p0=1000.0!hPa
    
    contains
SUBROUTINE SOR(nlev,nlat,nlon,dx,dy,dp,lat0,lat,Np,q,phi,psi)
    IMPLICIT NONE      
    INTEGER,INTENT(IN) :: nlev,nlat,nlon
    REAL,INTENT(IN) :: dx(nlat),dy,dp
    REAL,INTENT(IN) :: lat0,lat(nlat),Np(nlev) 
    REAL,INTENT(IN) :: q(nlon,nlat,nlev),phi(nlon,nlat,nlev)
      
    REAL,INTENT(OUT) :: psi(nlon,nlat,nlev)

    REAL,PARAMETER :: w=1.0, epsilon=10.0
    REAL :: f0
    REAL :: lam2(nlat),Ndp(0:nlev+1)
    REAL :: c(nlat,nlev),diff
    REAL :: pold(nlon,nlat,nlev),pbc(0:nlon+1,0:nlat+1,0:nlev+1)
    REAL :: rmse
    INTEGER :: I,J,K,iter
    DOUBLE PRECISION :: t1,t2

    t1 = omp_get_wtime()
    f0 = 2*omega*sin(lat0*deg2rad)
!$OMP PARALLEL DO
    DO J=1,nlat
       lam2(J) = dx(J)*dx(J)/dy/dy
    ENDDO
!$OMP END PARALLEL DO
    print *, "lam2",lam2
    Ndp(0) = 1.0/Np(1)
!$OMP PARALLEL DO
    DO K=1,nlev
        Ndp(K) = 1.0/Np(K) !N^(-2)(k)
    ENDDO
!$OMP END PARALLEL DO
    Ndp(nlev+1) = 1.0/Np(nlev)
    print *, "Ndp",Ndp
!                    
! calcurate coeff.
!
!$OMP PARALLEL DO PRIVATE(J)
    do K = 1, nlev
        do J = 1, nlat
            c(J,K) = 2.0*(1.0+lam2(J)) + &
            &   f0*f0*dx(J)*dx(J)*(Ndp(K+1)+2*Ndp(K)+Ndp(K-1))/2/dp/dp
        end do
    end do
!$OMP END PARALLEL DO
!
! Initialization
!
    psi = 0.0
!
! Start iteration
!
    DO iter = 1, 20000
!            
! (0) copy old version
! 
    pold = psi
!
! (1) boundary
!   
    pbc(1:nlon,1:nlat,1:nlev) = psi
    !vertical
    pbc(1:nlon,1:nlat,0) = psi(:,:,2)
    pbc(1:nlon,1:nlat,nlev+1) = psi(:,:,nlev-1)
    ! horizontal
    psi(1,:,:) = phi(1,:,:)/f0
    psi(nlon,:,:) = phi(nlon,:,:)/f0
    psi(:,1,:) = phi(:,1,:)/f0
    psi(:,nlat,:) = phi(:,nlat,:)/f0
!$OMP PARALLEL DO PRIVATE(I,J,diff)
    DO K = 1, nlev
        DO J = 2, nlat-1
            DO I = 2, nlon-1
!
! (2) calcurate diff
!
    diff = pbc(I+1,J,K)+pbc(I-1,J,K)+lam2(J)*(pbc(I,J+1,K)+pbc(I,J-1,K)) &  
    &    + f0*f0*dx(J)*dx(J)*((Ndp(K+1)+Ndp(K))*pbc(I,J,K+1) + (Ndp(K)+Ndp(K-1))*pbc(I,J,K-1))/2/dp/dp &
    &    - dx(J)*dx(J)*q(I,J,K)
!     
! (3) update
!
    psi(I,J,K) = (1-w)*pold(I,J,K) + w*diff/c(J,K)
            ENDDO
        ENDDO
    ENDDO
!$OMP END PARALLEL DO
!
! (4) converge check
!
    rmse = 0.0
!$OMP PARALLEL DO PRIVATE(I,J) REDUCTION(+:rmse)
    DO K = 1, nlev
        DO J = 1, nlat
            DO I = 1, nlon
                rmse = rmse + (psi(I,J,K)-pold(I,J,K))**2
            ENDDO
        ENDDO
    ENDDO
!$OMP END PARALLEL DO
    !print *, sqrt(rmse)
    IF(sqrt(rmse) .lt. epsilon) THEN 
        print*,"END",sqrt(rmse)
        EXIT
    ENDIF
    IF(mod(iter,100)==0)print*, sqrt(rmse)
    ENDDO
    t2 = omp_get_wtime()
    print *, "time",t2-t1

    RETURN

END SUBROUTINE SOR
! simplify(dx=dy,dp=const.,Np=const.,b.c.=Neumann)
SUBROUTINE SOR_s(nlev,nlat,nlon,dx,dy,dp,Np,lat0,q,phi,psi)
    IMPLICIT NONE      
    INTEGER,INTENT(IN) :: nlev,nlat,nlon
    REAL,INTENT(IN) :: dx,dy,dp,Np,lat0
    REAL,INTENT(IN) :: q(nlon,nlat,nlev),phi(nlon,nlat,nlev)
      
    REAL,INTENT(OUT) :: psi(nlon,nlat,nlev)

    REAL,PARAMETER :: w=1.0, epsilon=10.
    REAL :: f0
    REAL :: fcori(nlat),lam(nlat),Ndp(nlev)
    REAL :: c,diff
    REAL :: pold(nlon,nlat,nlev)
    REAL :: pbc(0:nlon+1,0:nlat+1,0:nlev+1)
    REAL :: rmse
    INTEGER :: I,J,K,iter
    DOUBLE PRECISION :: t1,t2

    t1 = omp_get_wtime()
    f0 = 2*omega*sin(lat0*deg2rad)
    c = 4 + 2*f0*f0*dy*dy/Np/dp/dp
    print *, c, f0,Np,dy,dp,2*f0*f0*dy*dy/Np/dp/dp
!
! Initialization
!
    psi = 0.0
!
! Start iteration
!
    DO iter = 1, 1000
!            
! (0) copy old version
! 
    pold = psi 
!
! (1) boundary
!   
    pbc(1:nlon,1:nlat,1:nlev) = psi
    !vertical
    pbc(1:nlon,1:nlat,0) = psi(:,:,2)
    pbc(1:nlon,1:nlat,nlev+1) = psi(:,:,nlev-1)
    ! horizontal
    psi(1,:,:) = phi(1,:,:)/f0
    psi(nlon,:,:) = phi(nlon,:,:)/f0
    psi(:,1,:) = phi(:,1,:)/f0
    psi(:,nlat,:) = phi(:,nlat,:)/f0
!$OMP PARALLEL DO PRIVATE(I,J,diff)
    DO K = 1, nlev
        DO J = 1, nlat
            DO I = 1, nlon
!
! (2) calcurate diff
!
    diff = (pbc(I+1,J,K)+pbc(I-1,J,K)+pbc(I,J+1,K)+pbc(I,J-1,K)) &
    &    + f0*f0*dy*dy*(pbc(I,J,K+1) + pbc(I,J,K-1))/Np/dp/dp &
    &    - dy*dy*q(I,J,K)
!     
! (3) update
!
    psi(I,J,K) = (1-w)*pold(I,J,K) + w*diff/c
            ENDDO
        ENDDO
    ENDDO
!$OMP END PARALLEL DO
!
! (4) converge check
!
    rmse = 0.0
!$OMP PARALLEL DO PRIVATE(I,J) REDUCTION(+:rmse)
    DO K = 1, nlev
        DO J = 1, nlat
            DO I = 1, nlon
                rmse = rmse + (psi(I,J,K)-pold(I,J,K))**2
            ENDDO
        ENDDO
    ENDDO
!$OMP END PARALLEL DO
    !print *, sqrt(rmse)
    IF(sqrt(rmse) .lt. epsilon) THEN 
        print*,"END",sqrt(rmse)
        EXIT
    ENDIF
    IF(mod(iter,100)==0)print*, sqrt(rmse)
    ENDDO
    t2 = omp_get_wtime()
    print *, "time",t2-t1

    RETURN

END SUBROUTINE SOR_s
end module sor_m