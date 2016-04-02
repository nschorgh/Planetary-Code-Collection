program exospherebody
!***********************************************************************
! surface temperatures of airless body with H2O exosphere
!***********************************************************************
  use grid
  use body, only: solarDay, dtsec, Rbody
  implicit none

  integer, parameter :: Np=40000  ! maximum number of computational particles

  integer n, i, nequil, k
  integer cc(6), cc_prod, cc_prod_total, ccc(4)
  integer :: idum=-92309   ! random number seed
  real(8) tmax, time, tequil
  real(8), dimension(veclen) :: Tsurf, Qn
  real(8) residencetime, HAi

  real(8), dimension(np,2) :: p_r ! longitude(1) and latitude(2)
  integer, dimension(np) :: p_s ! status 0=on surface, 1=inflight, <0= destroyed or trapped
  real(8), dimension(np) :: p_t ! time
  integer, dimension(np) :: p_n ! # of hops (diagnostic only)

  integer, external :: inbox, totalnr
  real(8), external :: residence_time
  real(8), external :: flux_noatm
  
  ! equilibration time for thermal model in solar days
  tequil=15.

  ! run time for hopping+thermal model in solar days
  tmax=4442.*1.5
  !tmax=2./29.53
  !tmax = 2*86400./solarDay

  ! set some constants
  nequil=int(tequil*solarDay/dtsec)

  print *,'Model parameters:'
  print *,'Time step=',dtsec,'sec'
  print *,'Equilibration time',tequil,' solar days'
  print *,'Maximum time',tmax,' solar days',tmax*solarDay/(86400.*365.242),' years'
  print *,'grid size',nlon,'x',nlat,'veclen=',veclen
  print *,'number of molecules',np

  HAi = 0.  ! noon

  ccc(:)=0
  cc_prod_total=0
  p_n(:) = 0

  ! initial configuration
  !p_t(:)=0.;  p_s(:)=0    ! all particles on surface
  p_t(:)=1d99; p_s(:)=-9   ! no particles
  !do i=1,np
     !p_r(i,1)=360.*ran2(idum); 
     !p_r(i,1)=0.
     !p_r(i,2)=0.; 
  !enddo

  open(unit=20,file='Tsurface',status='unknown',action='write')
  open(unit=21,file='Tsurface_end',status='unknown',action='write')
  open(unit=30,file='series',status='unknown',action='write')
  open(unit=50,file='particles',status='unknown',action='write')
  open(unit=51,file='particles_end',status='unknown',action='write')

  !-------loop over time steps 
  do n=-nequil,1000000
     time =(n+1)*dtsec   ! time at n+1 
     if (time>tmax*solarDay) exit

     !-- Temperature
     call SurfaceTemperature(dtsec,HAi,time,Tsurf,Qn) ! model T
     if (n<0) cycle   ! skip remainder

     ! some output
     call totalnrs(np,p_s,cc)
     print *,time/3600.,'Ten minute call',sum(cc(1:2))
     write(30,*) time/3600.,cc(1:2),ccc(1:4),cc_prod_total

     ! create new particles
     call production(Np,p_r,p_s,p_n,idum,Tsurf,cc_prod)
     !cc_prod = 0
     cc_prod_total = cc_prod_total + cc_prod

     ! update residence times with new temperature
     do i=1,np
        if (p_s(i)==0) then ! on surface
           k = inbox(p_r(i,:))
           residencetime = residence_time(Tsurf(k))
           p_t(i) = residencetime   ! relative to time
        endif
     enddo
     
     if (n==0) then ! write out initial distribution
        !call writeparticles(50,np,p_r,p_s,p_t,p_n)
        call writeglobe(20,Tsurf)
     endif

     ! 1 hour of hopping
     call montecarlo(np,idum,p_r,p_s,p_t,p_n,Tsurf,dtsec,ccc,Qn)
     !call destruction(np,p_r,p_s,p_t,idum,dtsec,veclen,sigma)

     call totalnrs(Np,p_s,cc)

  enddo
!50 continue
  call writeparticles(51,Np,p_r,p_s,p_t,p_n)
  call writeglobe(21,Tsurf)

  call totalnrs(Np,p_s,cc)
  print *,'# particles on surface',cc(1)
  print *,'# particles in flight',cc(2)
  print *,'# particles destroyed, photo',ccc(1)
  print *,'# particles destroyed, escape',ccc(2)
  print *,'# particles coldtrapped',ccc(3),ccc(4),ccc(3)+ccc(4)
  print *,'# particles produced',cc_prod_total
  print *,'# active particles',sum(cc(1:6)),Np
  !print *,'# produced / surface density ',cc_prod_total,cc_prod_total/(4*pi*Rbody**2)

  close(20); close(21)
  close(30)
  close(40); close(41)
end program exospherebody


subroutine SurfaceTemperature(dtsec,HAi,time,Tsurf,Qn)
  ! surface temperature model
  use body, only: solarDay, zmax, Fgeotherm, semia, albedo, emiss
  use grid, only: VECLEN
  implicit none
  real(8), intent(IN) :: dtsec, HAi, time
  real(8), intent(INOUT) :: Tsurf(veclen)
  real(8), intent(OUT) :: Qn(veclen)

  integer, parameter :: NMAX=1000
  real(8), parameter :: pi=3.1415926535897932, d2r=pi/180.
  real(8), parameter :: sigSB=5.6704d-8
  real(8), parameter :: zero=0.

  real(8), save :: ti(NMAX), rhocv(NMAX), z(NMAX)
  real(8), save :: T(NMAX,veclen)
  logical, save :: FirstCall = .TRUE.

  integer, parameter :: nz=30
  real(8), parameter :: zfac=1.07d0

  integer k, i
  real(8) thIn, rhoc, delta, decl, Qnp1, sunR
  real(8) lon,lat,HA,geof,Tmean,Fsurf
  real(8) eps, ecc, omega, Ls
  parameter(ecc = 0.075822766, eps = 4.*d2r, omega=301.*d2r) ! Ceres
  real(8), external :: flux_noatm

  call generalorbit(time/86400.,semia,ecc,omega,eps,Ls,decl,sunR)

  ! initialization
  if (FirstCall) then
     !thIn= 50.;  rhoc=1000000.  
     thIn= 15.;  rhoc=1200.*500.  ! Ceres
     delta = thIn/rhoc*sqrt(solarDay/pi)  ! skin depth

     print *,'Thermal model parameters:'
     print *,'nz=',nz,' zmax=',zmax,' zfac=',zfac
     print *,'Thermal inertia=',thIn,' rho*c=',rhoc
     print *,'Geothermal flux=',Fgeotherm
     print *,'Diurnal skin depth=',delta
     print *,'Albedo=',albedo

     ! set up depth grid
     call setgrid(nz,z,zmax,zfac)
     if (z(6)>delta) then
        print *,'WARNING: less than 6 points within diurnal skin depth'
     endif
     do i=1,nz
        if (z(i)<delta) cycle
        print *,i-1,' grid points within diurnal skin depth'
        exit
     enddo
     if (z(1)<1.e-5) print *,'WARNING: first grid point is too shallow'
     
     ti(1:nz) = thIn
     rhocv(1:nz) = rhoc

     do k=1,veclen
        call k2lonlat(k,lon,lat)
        geof = max(cos(lat),cos(lat+eps/2),cos(lat-eps/2))/pi
        Tmean=(1370*(1.-albedo)*geof/sigSB)**0.25 - 10.
        if (.not. Tmean>0.) Tmean=(Fgeotherm/sigSB)**0.25  ! fixes nan and zero
        T(1:nz,k) = Tmean
        Tsurf(k) = Tmean
        HA=2.*pi*mod((time-dtsec)/solarDay+(lon-HAi)/360.,1.d0)    ! hour angle
        Qn(k)=(1-albedo)*flux_noatm(sunR,decl,lat,HA,zero,zero)
     enddo
  end if

  do k=1,veclen
     call k2lonlat(k,lon,lat)
     ! longitude of morning terminator = -time/solarDay + lon + HAi ??
     HA=2.*pi*mod(time/solarDay+(lon-HAi)/360.,1.d0)    ! hour angle
     Qnp1=(1-albedo)*flux_noatm(sunR,decl,lat,HA,zero,zero)
     call conductionQ(nz,z,dtsec,Qn(k),Qnp1,T(:,k),ti,rhocv,emiss, &
          & Tsurf(k),Fgeotherm,Fsurf)
     Qn(k)=Qnp1
  enddo

  FirstCall = .FALSE.
end subroutine SurfaceTemperature
