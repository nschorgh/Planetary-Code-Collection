PROGRAM allofmoon
!***********************************************************************
! surface temperatures of airless body with H2O exosphere
!***********************************************************************
  use grid
  use body, only: solarDay, dtsec, Rmoon => Rbody
  implicit none

  integer, parameter :: Np=2000000  ! maximum number of computational particles

  integer n, i, nequil, k
  integer cc(6), cc_prod, cc_prod_total, ccc(4)
  integer :: idum=-92309   ! random number seed
  real(8) tmax, time, tequil
  real(8), dimension(veclen) :: Tsurf, Qn !, sigma 
  real(8) residencetime, HAi
  character(10) ext

  real(8), dimension(np,2) :: p_r ! longitude(1) and latitude(2)
  integer, dimension(np) :: p_s ! status 0=on surface, 1=inflight, <0= destroyed or trapped
  real(8), dimension(np) :: p_t ! time
  integer, dimension(np) :: p_n ! # of hops (diagnostic only)

  integer, external :: inbox, totalnr
  real(8), external :: residence_time
  real(8), external :: flux_noatm, ran2
  
  ! if used with Diviner surface temperature maps as input
  !integer, parameter :: NT=708  ! tied to timestep of 708/(29.53*86400.) = 3596 sec
  !real(8) lon(nlon),lat(nlat),Tbig(nlon,nlat,24)

  ! equilibration time for thermal model in solar days
  tequil=10.  ! model T
  !tequil=0.  ! Diviner T

  ! run time for hopping+thermal model in solar days
  tmax=4.
  !tmax=2./29.53
  !tmax = 2*86400./solarDay

  ! set some constants
  nequil=int(tequil*solarDay/dtsec)

  print *,'Model parameters:'
  print *,'Time step=',dtsec,'sec'
  print *,'Equilibration time',tequil,' lunar days'
  print *,'Maximum time',tmax,' lunar days',tmax*solarDay/(86400.*365.242),' years'
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
  open(unit=30,file='series',status='unknown',action='write')  ! time series
  open(unit=50,file='particles',status='unknown',action='write')
  open(unit=51,file='particles_end',status='unknown',action='write')

  ! if you work with real coldtrap contours
  !call readcontours
  
  ! if Diviner surface temperatures are read in
  !call readTmaps(nlon,nlat,lon,lat,Tbig)

  !-------loop over time steps 
  do n=-nequil,1000000
     time =(n+1)*dtsec   ! time at n+1 
     if (time>tmax*solarDay) exit

     !-- Temperature
     call SurfaceTemperature(dtsec,HAi,time,Tsurf,Qn) ! model T
     !call interpTmaps(nlon,nlat,NT,n,lon,Tbig,Tsurf,Qn) ! Diviner T
     if (n<0) cycle   ! skip remainder

     ! some output
     call totalnrs(np,p_s,cc)
     print *,time/3600.,'One hour call',sum(cc(1:2))
     write(30,*) time/3600.,cc(1:2),ccc(1:4),cc_prod_total

     ! create new particles
     call production(Np,p_r,p_s,p_n,idum,cc_prod,Tsurf)
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

     write(ext,'(i0.4)') n
     call deblank(ext)
     !open(unit=27,file='particles.'//ext,status='unknown',action='write')
     !call writeparticles(27,np,p_r,p_s,p_t,p_n)
     !close(27)
     !open(unit=28,file='tsurf.'//ext,status='unknown',action='write')
     !call writeglobe(28,Tsurf)
     !close(28)
     !if (time>(tmax-1.)*solarDay) then
     !   call particles2sigma(Np,p_r,p_s,sigma)
     !endif

     ! 1 hour of hopping
     call montecarlo(np,idum,p_r,p_s,p_t,p_n,Tsurf,dtsec,ccc,Qn)

     call totalnrs(Np,p_s,cc)

  enddo
!50 continue
  call writeparticles(51,Np,p_r,p_s,p_t,p_n)
  !call particles2sigma(Np,p_r,p_s,sigma)
  call writeglobe(21,Tsurf)

  call totalnrs(Np,p_s,cc)
  print *,'# particles on surface',cc(1)
  print *,'# particles in flight',cc(2)
  print *,'# particles destroyed, photo',ccc(1)
  print *,'# particles destroyed, escape',ccc(2)
  print *,'# particles coldtrapped',ccc(3),ccc(4),ccc(3)+ccc(4)
  print *,'# particles produced',cc_prod_total
  print *,'# active particles',sum(cc(1:6)),Np
  print *,'# produced / surface density ',cc_prod_total,cc_prod_total/(4*pi*Rmoon**2)

  close(20); close(21)
  close(30)
  close(40); close(41)
END PROGRAM allofmoon



subroutine deblank(chr)
! deblank string
! posted by James Giles at comp.lang.fortran on 3 August 2003.
  implicit none
  character(*), intent(inout) :: chr
  integer i, j, istart
  
  istart = index(chr, " ")
  if (istart == 0) return   ! no blanks in the argument
  
  j = istart-1  ! just before the first blank.
  
  do i = istart, len(trim(chr))
     if (chr(i:i) /= " ") then
        j = j+1
        chr(j:j) = chr(i:i)
     endif
  end do
  
  if(j < len(chr)) chr(j+1:) = " "  ! clear the rest of the string
  return
end subroutine deblank



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
  !real(8) eps, ecc, omega, Ls

  real(8), external :: flux_noatm

  ! toy orbit
  decl = 0.
  sunR = semia
  ! better orbit 
  !eps = 1.54*d2r    ! lunar obliquity to ecliptic 
  !ecc=0.; omega=0.
  !call generalorbit(time/86400.,semia,ecc,omega,eps,Ls,decl,sunR)

  ! initialization
  if (FirstCall) then
     thIn= 50.;  rhoc=1000000.
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
        geof = cos(lat)/pi
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



subroutine particles2sigma(Np, p_r, p_s, sigma)
  ! calculates areal density sigma from particle coordinates
  use grid, only: veclen
  use body, only: Rbody
  implicit none
  integer, intent(IN) :: Np, p_s(Np)
  real(8), intent(IN) :: p_r(Np,2)
  real(8), intent(OUT) :: sigma(veclen)
  integer i, k, nr0(veclen), nr1(veclen), totalnr0, totalnr1
  real(8) dA(veclen)
  logical, save :: FirstCall = .TRUE.
  integer, external :: inbox

  nr0(:)=0;  nr1(:)=0
  do i=1,Np
     if (p_s(i)==0) then
        k = inbox(p_r(i,:))
        nr0(k) = nr0(k)+1
     endif
     if (p_s(i)==1) then
        k = inbox(p_r(i,:))
        nr1(k) = nr1(k)+1
     endif
  enddo

  call areas(dA)
  dA = dA*Rbody**2

  sigma(:) = nr0(:)/dA(:)
  !sigma(:) = nr1(:)/dA(:)  
  !print *,'total area',sum(dA)/(4*pi*Rbody**2)  ! test
  totalnr0 = sum(nr0(:))
  totalnr1 = sum(nr1(:))
  !print *,'total # particles in particles2sigma:',Np,totalnr0  ! for checking purposes

  !where(sigma>maxsigma) maxsigma=sigma
  
  if (FirstCall) then 
     open(unit=40,file='sigma.dat',action='write')
     FirstCall = .FALSE.
     !maxsigma=0.
  else
     open(unit=40,file='sigma.dat',action='write',position='append')
  endif
  write(40,'(999999(1x,g11.5))') sigma
  close(40)
end subroutine particles2sigma



subroutine fallmap(unit, fall)
  ! optional analysis
  use body, only: Rmoon => Rbody
  use grid
  implicit none
  integer, intent(IN) :: unit, fall(veclen)
  integer i,j,k
  real(8) sigma(veclen), dA(veclen)
  real(8) longitude(nlon), latitude(nlat)

  call areas(dA)
  dA = dA*Rmoon**2

  sigma(:) = fall(:)/dA(:)

  call lonlatgrid(longitude,latitude)

  write(unit,110) 0.,90.,sigma(1)
  do j=1,nlat
     do i=1,nlon
        k = 1 + i + (j-1)*nlon 
        write(unit,110) longitude(i),latitude(j),sigma(k)
     enddo
  enddo
  write(unit,110) 0.,-90.,sigma(veclen)
110 format (f5.1,1x,f6.2,1x,g10.4)

  print *,'Total area',sum(dA)
end subroutine fallmap

