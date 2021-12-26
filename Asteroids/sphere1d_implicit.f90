PROGRAM sphere1d_implicit
! 1D spherically symmetric heat equation and ice retreat
! solved with semi-implicit method
  implicit none
  integer, parameter :: nz = 100    ! number of grid points
  integer i, ierr

  real(8) semia, ecc, Q, time0
  real(8), parameter :: So = 1365    ! solar constant [W/m^2]
  real(8), parameter :: A = 0.05     ! albedo
  real(8), parameter :: emiss = 0.96 ! infrared emissivity
  
  logical :: init = .true.
  real(8) z(nz), dr, kappa, Radius, dtsec, time, oldtime, Deltat
  real(8) T(nz), Tsurf, Tsurfm1, Teff

  real(8) zT  ! depth of ice table
  real(8) Tatz, Elatent
  real(8), external :: interp1, flux2T

  zT = 0.
  time = 0.  ! earliest time in input file
  
  Radius = 500.  ! body radius [m]
  dr = Radius/nz
  do i=1,nz
     z(i) = i*dr
  end do

  kappa = 0.25/(0.5*2600*500)  ! thermal diffusivity kappa=k/rhoc [m^2/s]

  write(*,*) 'albedo=',A,'emissivity=',emiss
  write(*,*) 'Radius=',Radius,'spatial resolution=',dr
  write(*,*) 'thermal diffusivity=',kappa
  write(*,*) 'initial ice depth=',zT

  ! example input file provided by Henry Hsieh
  open(unit=20,file='orbit_A0000240.dat',action='read',iostat=ierr,status='old')
  if (ierr>0) stop 'input file not found'

  if (nz>1000) stop 'tridag is only set up for N<=1000'
  open(unit=31,file='depths_sphere.dat',action='write')

  read(20,*) ! skip headerline
  
  do  ! time loop
     read(20,*,iostat=ierr) time0,semia,ecc
     if (ierr<0) stop 'reached end of file'
     if (ierr>0) stop 'read error'

     oldtime = time
     time = time0
     Deltat = time0-oldtime  ! usually 1000 yrs
     !print *,Deltat
     if (Deltat<0. .or. Deltat>1001.) stop 'Deltat'
     dtsec = Deltat*86400*365.24  ! yr -> sec
     
     !--orbital elements -> surface temperature
     Q = So/semia**2/sqrt(1-ecc**2) ! annual mean insolation
     !Q = Q*faintsun(1e8-time)
     if (.not.init) Tsurfm1 = Tsurf
     Teff = flux2T(Q/4.,A,emiss)

     ! options for surface temperature
     Tsurf = Teff ! upper bound
     !call insolonly1(latitude,semia,omega,ecc,obliq,Q0mean,Qmean,Q4mean)
     !print *,Tsurf,flux2T(Qmean,A,emiss),flux2T(Q4mean,A,emiss)
     !Tsurf = flux2T(Q4mean,A,emiss)  ! cold end-member
     
     if (init) then ! initialization
        print *,'initializing particle',Tsurf
        T(:) = Tsurf
        Tsurfm1 = Tsurf
        init = .false.
     end if
     
     !--thermal evolution
     call conductionT_sphere(nz,dr,dtsec,T(:),Tsurfm1,Tsurf,kappa)
     
     !--ice evolution
     if (zT>=0. .and. zT<Radius) then
        Tatz = interp1(real(0.,8),z,Tsurf,T(:),zT,nz)
        call retreat_s(Tatz,zT,Radius,dtsec,Elatent)
        if (zT>Radius) zT=-9999.
     else
        Tatz = -9999.
     end if

     ! lots of output
     write(31,'(f11.0,1x,f7.4,1x,f6.4,3(1x,f5.1),1x,f8.3)') &
          & time,semia,ecc,Tsurf,T(nz),Tatz,zT
     
  end do
  close(20)
  
  close(31)

  ! radial temperature profile at end of run
  open(unit=40,file='z.dat',action='write')
  do i=1,nz
     write(40,*) z(i),T(i)
  end do
  close(40)
END PROGRAM sphere1d_implicit



subroutine conductionT_sphere(nz,dr,dt,T,Tsurf,Tsurfp1,kappa)
!***********************************************************************
!   conductionT_sphere: solve heat equation in a
!                       spherically symmetric geometry
!   Crank-Nicolson scheme, flux conservative
!
!   Eqn: T_t = (kappa/r^2)*d/dr(r^2*dT_r)
!   BC (z=0, r=R): T=T(t)
!   BC (z=R, r=0): center of sphere, dT_r=0
!
!   kappa = thermal diffusivity [m^2/s]
!   T = radial temperature profile [K]
!   Tsurf, Tsurfp1 = surface temperatures at times n and n+1
!
!   Grid: surface is at z=0
!         T(nz) is at center of sphere
!***********************************************************************
  implicit none
  integer, intent(IN) :: nz
  real(8), intent(IN) :: dr, dt, Tsurf, Tsurfp1, kappa
  real(8), intent(INOUT) :: T(nz)
  integer i
  real(8) cp(nz), cm(nz), cn(nz)
  real(8) alpha, a(nz), b(nz), c(nz), r(nz)
  
  alpha = kappa*dt/dr**2

  do i=1,nz-1
     r(i) = (nz-i)*dr
     cp(i) = (r(i)-dr/2)**2/r(i)**2  ! inward of r(i)
     cm(i) = (r(i)+dr/2)**2/r(i)**2  ! outward of r(i)
     cn(i) = (cp(i)+cm(i))/2.
  enddo
  
  ! elements of tridiagonal matrix
  do i=1,nz-1
     a(i) = -alpha/2*cm(i)  !  a(1) is not used
     b(i) = 1. + alpha*cn(i)
     c(i) = -alpha/2*cp(i)  !  c(nz) is not used
  enddo 
  a(nz) = -3*alpha 
  b(nz) = 1. + 3*alpha
  
  ! Set RHS         
  r(1) = alpha/2*cp(1)*T(2) + (1.-alpha*cn(1))*T(1) + alpha/2*cm(1)*(Tsurf+Tsurfp1)
  do i=2,nz-1
     r(i) = alpha/2*cp(i)*T(i+1) + (1.-alpha*cn(i))*T(i) + alpha/2*cm(i)*T(i-1)
  enddo
  r(nz) = 3*alpha*T(nz-1) + (1.-3*alpha)*T(nz)

  ! Solve for T at n+1
  call tridag(a,b,c,r,T,nz) ! update by tridiagonal inversion
  
end subroutine conductionT_sphere



subroutine retreat_s(T,zT,Radius,dt,Elatent)
  ! retreat of ice table, 1D spherically symmetric
  implicit none
  real(8), intent(IN) :: T  ! scalar temperature [K]
  real(8), intent(IN) :: Radius  ! radius of body [m]
  real(8), intent(IN) :: dt  ! time step [sec]
  real(8), intent(INOUT) :: zT  ! depth of ice table below surface
  real(8), intent(OUT) :: Elatent  ! latent heat, for diagnostics
  real(8), parameter :: pi = 3.1415926535897932
  real(8), parameter :: Lh2o = 2.834e6 ! latent heat of sublimation [J/kg]
  real(8), parameter :: R = 8314.5 ! universal gas constant
  real(8), parameter :: rhoice = 930. ! [kg/m^3]
  real(8) D, rhos, buf, diam, porosity
  real(8), external :: psv, vapordiffusivity
  real(8) zTold, dV
  
  diam = 100d-3; porosity = 0.5
  D=vapordiffusivity(diam,porosity,T)
  
  rhos = psv(T)*18/(R*T) ! p = n k T
  buf = D*rhos/(rhoice*porosity)/(1-zT/Radius)
  zTold = zT
  zT = sqrt(zT**2 + 2*buf*dt) ! fixes divergence at surface
  if (zT>=0. .and. zTold>=0.) then
     dV = (zT-zTold)*4*pi*(Radius-(zT+zTold)/2)**2
     !dV = 4*pi/3*((Radius-zTold)**3 - (Radius-zT)**3) ! more accurate but roundoff sensitive
     Elatent = dV*porosity*rhoice*Lh2o  ! [J]
  else
     Elatent = -999.
  end if
end subroutine retreat_s


