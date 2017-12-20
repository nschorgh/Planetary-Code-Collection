! Mars thermal model with horizons from 3D topography

module miscparams
  real(8), parameter :: pi=3.1415926535897932, d2r=pi/180.
  real(8), parameter :: sigSB = 5.6704e-8
  real(8), parameter :: Tco2frost=145., Lco2frost=6.0e5  ! Mars
  real(8), parameter :: zero = 0.

  real(8), parameter :: fracIR=0.04, fracDust=0.02
  real(8), parameter :: solsy = 668.60 ! solar days per Mars year
  real(8), parameter :: solarDay = 88775.244, Fgeotherm = 0.028  ! Mars
  integer, parameter :: nz=70
end module miscparams


program cratersQ_Mars
  !use omp_lib
  use filemanager
  use allinterfaces
  use miscparams
  use newhorizons
  implicit none

  real(8) ecc, eps, omega
  real(8) edays, marsR, marsLs, marsDec
  real(8) m(NSx,NSy), dE
  real(8), parameter :: albedo0=0.15, co2albedo=0.65
  real(8), parameter :: a = 1.52366 ! Mars semimajor axis in A.U.
  real(8), parameter :: emiss = 1.

  integer nsteps, n, i, j, nm, i0, j0, k
  real(8) tmax, dt, latitude, dtsec, buf
  real(8) HA, sdays, azSun, emax, sinbeta
  real(8), allocatable, dimension(:,:) :: h, surfaceSlope, azFac
  real(8), allocatable, dimension(:,:) :: Qn   ! incoming
  real(8), allocatable, dimension(:,:) :: Tsurf, albedo, Fsurf
  real(8), allocatable, dimension(:,:) :: Qmean, Qmax, Tmean, Tmaxi, Tbottom
  real(8), allocatable, dimension(:,:) :: mmax, frosttime, maxfrosttime, Qnm1
  real(8), allocatable :: T(:,:,:)  ! subsurface
  real(8) Tsurfold
  logical, parameter :: subsurface=.false.  ! control panel
  integer, parameter, dimension(4) :: i00=(/ 41, 42, 44, 74/), j00=(/108, 108, 109, 155/)
  integer, parameter :: MARGIN=20  ! must be at least 1, saves time

  allocate(h(NSx,NSy), surfaceSlope(NSx,NSy), azFac(NSx,NSy))
  allocate(Qn(NSx,NSy), Tsurf(NSx,NSy), albedo(NSx,NSy), Fsurf(NSx,NSy))
  allocate(Qmean(NSx,NSy), Qmax(NSx,NSy), Tmean(NSx,NSy), Tmaxi(NSx,NSy), Tbottom(NSx,NSy))
  allocate(mmax(NSx,NSy), frosttime(NSx,NSy), maxfrosttime(NSx,NSy), Qnm1(NSx,NSy))

  ecc = 0.0934;  eps = 25.19*d2r;  omega = 250.87*d2r   ! today
  
  dt=0.02; 
  tmax = solsy+1.
  !tmax = solsy*10.5
  tmax = 2.
  latitude = -41.6
  albedo(:,:) = albedo0

  ! set some constants
  nsteps=int(tmax/dt)       ! calculate total number of timesteps
  dtsec = dt*solarDay
  
  write(*,*) 'Time step=',dt,' Max number of steps=',nsteps
  write(*,*) 'Calculations performed for latitude=',latitude
  write(*,*) 'Nx=',NSx,'Ny=',NSy,'File=',fileext
  write(*,*) 'Mean albedo=',sum(albedo)/size(albedo),'Emissivity=',emiss
  write(*,*) 'Reflections:',.FALSE.,'Subsurface:',subsurface

  ! setenv OMP_NUM_THREADS 4
  !write (*,'(a,i8)') 'The number of processors available = ', omp_get_num_procs()
  !write (*,'(a,i8)') 'The number of threads available    = ', omp_get_max_threads()

  call readdem(h)
  call difftopo(NSx,NSy,h,dx,dy,surfaceSlope,azFac)

  latitude=latitude*d2r
  Tsurf=200.
  Qmean=0.; Tmean=0.; Tbottom=0.; nm=0
  Qmax=0.; Tmaxi=0.; mmax=0.
  frosttime=0.; maxfrosttime=0.
  m=0.; Fsurf=0.
  
  print *,'...reading horizons file...'
  call readhorizons('Data2/horizons.'//fileext)

  if (subsurface) then
     allocate(T(NSx,NSy,1000))
     call subsurfaceconduction_mars(T(1,1,:),buf,zero,zero,zero,zero,buf,buf,.true.)
  end if

  open(unit=22,file='timeseries.dat',status='unknown',action='write')

  print *,'...calculating...'
  ! loop over time steps 
  do n=0,nsteps-1
     sdays = (n+1)*dtsec/solarDay
     edays = (n+1)*dtsec/86400.
     
     call generalorbit(edays,a,ecc,omega,eps,marsLs,marsDec,marsR)
     HA=2.*pi*mod(sdays,1.)   ! hour angle
     call equatorial2horizontal(marsDec,latitude,HA,sinbeta,azSun)

     if (mod(n,10)==0) print *,n,sdays,HA
     
     !$OMP PARALLEL DO private(emax)
     do i=1+MARGIN,NSx-MARGIN
        do j=1+MARGIN,NSy-MARGIN
           emax = getonehorizon(i,j,azSun)
           Qn(i,j)=(1.-albedo(i,j))*flux_mars( &
                & marsR,marsDec,latitude,HA,0.d0,fracir,fracdust,surfaceSlope(i,j),azFac(i,j),emax)
        enddo
     enddo
     !$OMP END PARALLEL DO

     if (n==0) Qnm1(:,:) = Qn(:,:)
     if (subsurface) then
        if (n==0) Tsurf(:,:) = -9. ! no-value, max 3 digits
        !$OMP PARALLEL DO
        do i=1+MARGIN,NSx-MARGIN
           do j=1+MARGIN,NSy-MARGIN
              call subsurfaceconduction_mars(T(i,j,:),Tsurf(i,j), &
                   & dtsec,Qnm1(i,j),Qn(i,j),emiss,m(i,j),Fsurf(i,j),.false.)
           enddo
        enddo
        !OMP END PARALLEL DO

     else  ! no subsurface conduction
        !$OMP PARALLEL DO private(dE, Tsurfold)
        do i=1+MARGIN,NSx-MARGIN
           do j=1+MARGIN,NSy-MARGIN
              Tsurf(i,j) = (Qn(i,j)/emiss/sigSB)**0.25
              Tsurfold = (Qnm1(i,j)/emiss/sigSB)**0.25
              if (Tsurf(i,j)<Tco2frost.or.m(i,j)>0.) then   ! CO2 condensation
                 Tsurf(i,j)=Tco2frost
                 dE = - Qn(i,j) + emiss*sigSB*(Tsurf(i,j)**4 + Tsurfold**4)/2.
                 m(i,j) = m(i,j) + dtsec*dE/Lco2frost
              endif
           enddo
        enddo
        !$OMP END PARALLEL DO
     endif
     Qnm1(:,:) = Qn(:,:)

     where (Tsurf>Tco2frost.or.m<=0.)
        albedo = albedo0
     elsewhere
        albedo = co2albedo
     end where

     !if (sdays > tmax-1) then
     if (sdays > tmax-solsy) then
        Qmean(:,:) = Qmean(:,:) + Qn
        where (Qn>Qmax) Qmax=Qn
        Tmean = Tmean + Tsurf
        where (Tsurf>Tmaxi) Tmaxi=Tsurf
        where (m>mmax) mmax=m
        if (subsurface) Tbottom=Tbottom+T(:,:,nz)
        nm=nm+1

        do k=1,4
           i0=i00(k); j0=j00(k)
           write(22,'(f9.3,2(1x,i4),2x,f6.1,1x,f5.1,1x,f6.1)') &
                & sdays,i0,j0,Qn(i0,j0),Tsurf(i0,j0),m(i0,j0)
        enddo
     endif
     if (sdays > tmax-2*solsy) then  ! longest continuous period below 200K
        where (Tsurf<200.) 
           frosttime=frosttime+dt
        elsewhere
           frosttime=0.
        end where
        where (frosttime>maxfrosttime) maxfrosttime=frosttime
     endif

  enddo  ! end of time loop

  close(22)
  if (subsurface) deallocate(T)

  Qmean=Qmean/nm
  Tmean=Tmean/nm; Tbottom=Tbottom/nm
  
  open(unit=21,file='qmean.dat',status='unknown',action='write')
  do i=2,NSx-1
     do j=2,NSy-1
        write(21,'(2(i4,1x),f9.2,2x,f6.4,2(1x,f6.1),2(1x,f5.1),1x,f7.1,1x,f6.1)') &
             & i,j,h(i,j),surfaceSlope(i,j),Qmean(i,j),Qmax(i,j), &
             & Tmean(i,j),Tmaxi(i,j),mmax(i,j),maxfrosttime(i,j)
     enddo
  enddo
  close(21)
end program cratersQ_Mars



subroutine subsurfaceconduction_mars(T,Tsurf,dtsec,Qn,Qnp1,emiss,m,Fsurf,init)
  use allinterfaces, only : conductionQ, conductionT
  use miscparams
  implicit none
  integer, parameter :: NMAX=1000
  real(8), intent(INOUT) :: T(NMAX), Tsurf, m, Fsurf
  real(8), intent(IN) :: dtsec,Qn,Qnp1,emiss
  logical, intent(IN) :: init
  integer i
  real(8) zmax, zfac, Tinit, delta
  real(8) Fsurfold, dE, Tsurfold, Told(1:nz)
  real(8), save :: ti(NMAX), rhocv(NMAX), z(NMAX)

  if (init) then ! initialize grid
     ti(:) = 500.;  rhocv(:) = 1200.*800.  ! adjust
     zmax=8.; zfac = 1.05  ! adjust

     delta = ti(1)/rhocv(1)*sqrt(solarDay/pi)  ! skin depth

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
     open(unit=30,file='z',status='unknown');
     write(30,*) (z(i),i=1,nz)
     close(30)

     write(*,*) 'Subsurface model parameters'
     write(*,*) '   nz=',nz,' zmax=',zmax,' zfac=',zfac
     write(*,*) '   Thermal inertia=',ti(1),' rho*c=',rhocv(1)
     print *,'   Diurnal skin depth=',delta,' Geothermal flux=',Fgeotherm

     return
  endif
  
  if (Tsurf<=0.) then  ! initialize temperature profile
     if (Tsurf<=0.) Tinit=200.
     T(1:nz) = Tinit
     Tsurf = Tinit
  endif

  Tsurfold=Tsurf
  Fsurfold=Fsurf
  Told(1:nz)=T(1:nz)
  if (Tsurf>Tco2frost.or.m<=0.) then
     call conductionQ(nz,z,dtsec,Qn,Qnp1,T,ti,rhocv,emiss,Tsurf,Fgeotherm,Fsurf)
  endif
  if (Tsurf<Tco2frost.or.m>0.) then   ! CO2 condensation                                              
     T(1:nz)=Told
     call conductionT(nz,z,dtsec,T,Tsurfold,Tco2frost,ti, &
          &              rhocv,Fgeotherm,Fsurf)
     Tsurf=Tco2frost
     dE = (- Qn - Qnp1 + Fsurfold + Fsurf + &
          &           emiss*sigSB*(Tsurfold**4+Tsurf**4))/2.
     m = m + dtsec*dE/Lco2frost;
  endif

end subroutine subsurfaceconduction_mars



elemental function evap_ingersoll(T)
  ! Returns evaporation rate (kg/m^2/s)
  ! unused
  use allinterfaces, only : psv
  implicit none
  real(8) evap_ingersoll
  real(8), intent(IN) :: T
  real(8) p0,psat,D,Gbuf,rho,R,rhow,nu,drhooverrho,g
  !real(8), external :: psv

  p0=520  ! atmospheric pressure
  psat=psv(T)
  R=8314
  D=20e-4 ! in m^2/s 
  g=3.7
  rhow = psat*18/(R*T)
  rho = p0*44/(R*T)
  drhooverrho=(44-18)*psat/(44*p0-(44-18)*psat) ! Ingersoll (1970)
  !drhooverrho=(44-18)*psat/(44*(p0-e)) ! diverges at p0=e
  !nu=7e-4  ! kinematic viscosity of CO2
  nu=16e-4*(T/200)**1.5
  Gbuf=(drhooverrho*g/nu**2)**(1./3.);
  evap_ingersoll=0.17*D*rhow*Gbuf

end function evap_ingersoll
 

