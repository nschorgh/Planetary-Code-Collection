module miscparameters
  real(8), parameter :: pi=3.1415926535897932, d2r=pi/180.
  real(8), parameter :: earthDay=86400., marsDay=88775.244, solsperyear=668.60
  real(8), parameter :: sigSB=5.6704e-8, Lco2frost=6.0e5, co2albedo=0.60
end module miscparameters


program tempr_driver
!***********************************************************************
! Temperatures on Mars over the past 21 Myr
!
! written by Norbert Schorghofer 2007
! modernized 2020 -norbert
!***********************************************************************
  use miscparameters
  implicit none
  integer, parameter :: NP=1      ! number of sites
  !integer, parameter :: earliest=21000+1  ! start time [kyr]
  integer, parameter :: earliest=5000+1  ! start time [kyr]
  integer, parameter :: nz=80     ! number of subsurface grid points
  real(8), parameter :: zfac=1.05 ! progressive spacing of grid points
  integer i, k, iargc, ierr
  real(8) zmax, delta, z(nz), icetime
  real(8) p0(earliest)
  real(8), dimension(NP) :: latitude, albedo, thIn, rhoc
  real(8) ecc, omega, eps
  real(8), dimension(NP) :: Tb, Tmean1, Tmean3, Qmean
  real(8), dimension(earliest) :: lasktime, laskecc, laskomega, laskeps
  character(10) ext

  if (iargc() /= 1) stop 'USAGE: a.out ext'
  call getarg( 1, ext )

  thIn(:) = 250. ! will be overwitten by input
  albedo(:) = 0.2  ! will be overwritten by input

  open(unit=21,file='lats.'//ext,action='read',status='old',iostat=ierr)
  if (ierr /= 0) then
     print *,'File lats.'//ext,'not found'
     stop
  endif
  do k=1,NP
     read(21,*) latitude(k),albedo(k),thIn(k)
  enddo
  close(21)

  ! this optional rescaling makes skin depth independent of thIn, and
  ! hence the same zmax can be used for different thIn
  rhoc(:) = 1300.*550.*thIn(:)/200.

  ! set eternal grid
  zmax = 6.
  call setgrid(nz,z,zmax,zfac)
  open(unit=30,file='z.'//ext,action='write',status='unknown')
  write(30,'(999(f8.5,1x))') (z(i),i=1,nz)
  close(30)

  ! ecc = 0.0934;  eps = 25.19*d2r;  omega = 250.87*d2r  ! today
  ! Laskar orbital solution http://vo.imcce.fr/insola/earth/online/mars/mars.html
  open(20,file='INSOLN.LA2004.MARS.ASC',action='read',status='old')
  p0(:) = 600.
  do i = 1,earliest
     read(20,*) lasktime(i),laskecc(i),laskeps(i),laskomega(i)
  enddo
  close(20)

  print *,'RUNNING MARS ICE AGES - TEMPERATURES ONLY'
  print *,'Global model parameters'
  print *,'nz=',nz,' zfac=',zfac,'zmax=',zmax
  print *,'Starting at time',lasktime(earliest)*1000.,'years'
  print *,'Mean total pressure=',sum(p0(1:earliest))/earliest
  print *,'Number of sites=',NP
  do k=1,NP
     print *,'  Latitude (deg)',latitude(k)
     print *,'  rho*c (J/m^3/K)',rhoc(k),'   thermal inertia',thIn(k)
     delta = thIn(k)/rhoc(k)*sqrt(marsDay/pi)
     do i=1,nz
        if (z(i)<delta) cycle
        print *,'  ',i-1,' grid points within diurnal skin depth'
        exit
     enddo
     print *,'  ',zmax/(sqrt(solsperyear)*delta),'times seasonal skin depth'
  enddo

  ! initializations 
  Tb(:) = -9999.

  icetime = lasktime(earliest)*1000
  ecc = laskecc(earliest)
  omega = laskomega(earliest) + pi
  eps = laskeps(earliest)

  open(unit=35,file='out_geo.'//ext,action='write',status='unknown')
  open(unit=37,file='out_subsurf.'//ext,action='write',status='unknown')

  ! equilibrate Tb 
  call icesheet(nz,NP,latitude,albedo,thIn,rhoc,z, &
       &   ecc,omega,eps,p0(earliest),Tb,Tmean1,Tmean3,Qmean)

  ! History begins here
  do i = earliest,1,-1
     ecc = laskecc(i)
     omega = mod(laskomega(i) + pi,2*pi)
     eps = laskeps(i)
     call icesheet(nz,NP,latitude,albedo,thIn,rhoc,z, &
          &        ecc,omega,eps,p0(i),Tb,Tmean1,Tmean3,Qmean)
     icetime = lasktime(i)*1000.
     write(35,'(f11.0,2x,f6.2,2x,f7.5,2x,f5.1)') icetime,eps/d2r,ecc,omega/d2r
     do k=1,NP
        write(37,'(f11.0,2x,f6.2,2(2x,f6.2))') &
             & icetime,latitude(k),Tmean1(k),Tmean3(k)
     enddo
     print *,icetime
  enddo

  close(35); close(37)
end program tempr_driver




subroutine icesheet(nz,NP,latitude,albedo,thIn,rhoc,z,ecc,omega,eps, &
     & p0,Tb,Tmean1,Tmean3,Qmean)
  ! pass-through subroutine for similarity to big iceage code
  use miscparameters
  implicit none
  integer, intent(IN) :: nz, NP
  real(8), intent(IN) :: latitude(NP), albedo(NP), thIn(NP), rhoc(NP), z(nz)
  real(8), intent(IN) :: ecc, omega, eps, p0
  real(8), intent(INOUT) :: Tb(NP)
  real(8), intent(OUT), dimension(NP) :: Tmean1, Tmean3, Qmean
  integer k
  real(8) fracIR, fracDust, ti(nz), rhocv(nz)

  fracIR=0.04*p0/600.; fracDust=0.02*p0/600.
  do k=1,NP 
     ti(1:nz) = thIn(k)
     rhocv(1:nz) = rhoc(k)
     call ajsub(latitude(k)*d2r, albedo(k), nz, z, ti, rhocv, &
          &     fracIR, fracDust, p0, ecc, omega, eps, Tb(k), &
          &     Tmean1(k), Tmean3(k), Qmean(k))
  enddo 
end subroutine icesheet




subroutine ajsub(latitude, albedo0, nz, z, ti, rhocv, &
     &     fracIR, fracDust, p0, ecc, omega, eps, Tb, &
     &     Tmean1, Tmean3, Qmean)
!***********************************************************************
!  ajsub: A 1D thermal model
!***********************************************************************
  use miscparameters
  implicit none
  integer, intent(IN) :: nz
  real(8), intent(IN) :: latitude  ! in radians
  real(8), intent(IN) :: albedo0, z(nz)
  real(8), intent(IN) :: ti(nz), rhocv(nz), fracIR, fracDust, p0
  real(8), intent(IN) :: ecc, omega, eps
  real(8), intent(INOUT) :: Tb
  real(8), intent(OUT) :: Tmean1, Tmean3, Qmean

  real(8), parameter :: Fgeotherm = 0., dt = 0.02
  real(8), parameter :: a = 1.52366 ! Mars semimajor axis [a.u.]
  real(8), parameter :: emiss = 1.  ! emissivity
  integer, parameter :: EQUILTIME = 15. ! [Mars years]

  integer nsteps, n, nm
  real(8) T(nz), tmax, time, Qn, Qnp1, tdays
  real(8) marsR, marsLs, marsDec, HA
  real(8) Tsurf, Tco2frost, albedo, Fsurf, m, dE
  real(8) Told(nz), Fsurfold, Tsurfold, Tmean0
  real(8), external :: flux_mars77, tfrostco2
  
  tmax = EQUILTIME*solsperyear
  nsteps = int(tmax/dt)     ! calculate total number of timesteps

  ! Tco2frost = 150.
  Tco2frost = tfrostco2(p0)

  if (Tb<=0.) then
     !Tmean0 = 210.         ! black-body temperature of planet
     Tmean0 = (589.*(1.-albedo0)*cos(latitude)/pi/sigSB)**0.25 ! better estimate
  else
     Tmean0 = Tb
  endif
  
  albedo = albedo0
  if (Tmean0 < Tco2frost) Tmean0=Tco2frost
  T(1:nz) = Tmean0
  Tsurf = T(1)
  m=0.; Fsurf=0.

  Qmean = 0.
  nm=0; Tmean1=0.; Tmean3=0.

  time=0.
  call generalorbit(0.d0,a,ecc,omega,eps,marsLs,marsDec,marsR)
  HA = 2.*pi*time            ! hour angle
  Qn = flux_mars77(marsR,marsDec,latitude,HA,albedo,fracir,fracdust)
  !-----loop over time steps 
  do n=0,nsteps-1
     time = (n+1)*dt         ! time at n+1 
     tdays = time*(marsDay/earthDay) ! parenthesis may improve roundoff
     call generalorbit(tdays,a,ecc,omega,eps,marsLs,marsDec,marsR)
     HA = 2.*pi*mod(time,1.d0)  ! hour angle
     Qnp1 = flux_mars77(marsR,marsDec,latitude,HA,albedo,fracir,fracdust)
     Qmean = Qmean + Qnp1
     Tsurfold = Tsurf
     Fsurfold = Fsurf
     Told(1:nz) = T(1:nz)

     if (m<=0.) then
        call conductionQ(nz,z,dt*marsDay,Qn,Qnp1,T,ti,rhocv,emiss, &
             &           Tsurf,Fgeotherm,Fsurf)
     endif
     if (Tsurf<Tco2frost .or. m>0.) then ! CO2 frost
        T(1:nz) = Told(1:nz)
        call conductionT(nz,z,dt*marsDay,T,Tsurfold,Tco2frost,ti, &
             &              rhocv,Fgeotherm,Fsurf) 
        Tsurf = Tco2frost
        dE = (- Qn - Qnp1 + Fsurfold + Fsurf + &
             &           emiss*sigSB*(Tsurfold**4+Tsurf**4))/2.
        m = m + dt*marsDay*dE/Lco2frost
     endif
     if (Tsurf>Tco2frost .or. m<=0.) then
        albedo = albedo0
     else
        albedo = co2albedo
     endif
     Qn = Qnp1
     
     if (time >= tmax-nint(solsperyear)) then
        Tmean1 = Tmean1+Tsurf
        Tmean3 = Tmean3+T(nz)
        nm = nm+1
     endif

  enddo  ! end of time loop
  
  Tmean1 = Tmean1/nm; Tmean3 = Tmean3/nm
  Qmean = Qmean/nsteps

  Tb = T(nz)

end subroutine ajsub

