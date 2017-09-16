program cratersQ_earth
!***********************************************************************
!   cratersQ: program to calculate insolation with shadowing and
!             reflections  
!             Earth orbit, with optional Mauna Kea atmosphere and
!             optional 1D subsurface conduction  
!
!   written by Norbert Schorghofer 2010-2016
!***********************************************************************
  use filemanager
  use allinterfaces
  use dateformat
  use newhorizons
  implicit none
  real(8), parameter :: pi=3.1415926535897932, d2r=pi/180.
  real(8), parameter :: earthDay=86400.
  real(8), parameter :: sigSB = 5.6704e-8  

  integer nsteps, n, i, j, nm, k, CCMAX, iii, jjj
  real(8) tmax, edays, dtmin, latitude
  real(8) R, dZenithAngle, dAzimuth, longitude
  real(8) azSun, sinbeta, emax, emiss, v
  real(8), dimension(NSx,NSy) :: h, surfaceSlope, azFac
  real(8), dimension(NSx,NSy) :: Qn, QIR, Qrefl   ! incoming
  integer, dimension(NSx,NSy) :: cc
  integer(2), dimension(:,:,:), allocatable :: ii,jj
  real(4), dimension(:,:,:), allocatable :: dO12
  real(8), dimension(NSx,NSy) :: Tsurf, Qvis, skysize, Qabs, albedo, QIRin, QIRre
  real(8) Qmeans(NSx,NSy,4), Qmax(NSx,NSy), Tmean(NSx,NSy), Qref, Qrefm
  real(8) I0,D0,S0,unsd  ! atmosphere
  real(8), allocatable :: T(:,:,:), Qnm1(:,:)  ! subsurface conduction
  type(cTime) udtTime
  real(8), parameter :: zero=0.
  logical, parameter :: reflection=.true., subsurface=.false., atmosphere=.true.
  real(8), parameter :: Qother = 100.  ! (W/m^2) background (LW+downwelling)
  
  dtmin=15.
  tmax = 365.+1
  !tmax = 1.+1
  latitude = 19.82; longitude = -155.4683  ! Mauna Kea summit
  ! azimuth in degrees east of north, 0=north facing
  albedo(:,:) = 0.05
  emiss = 0.95

  ! set some constants
  nsteps=int(tmax*1440./dtmin)       ! calculate total number of timesteps

  ! start time in UTC = HST-10
  udtTime = cTime(2012,3,1,0.,0.,0.) ! 14 HST 
  
  write(*,*) 'Starting time',udtTime
  write(*,*) 'Time step=',dtmin,'(min)  Max number of steps=',nsteps
  write(*,*) 'Calculations performed for latitude=',latitude
  write(*,*) 'Nx=',NSx,'Ny=',NSy,'File=',fileext
  write(*,*) 'Mean albedo=',sum(albedo)/size(albedo),'Emissivity=',emiss
  write(*,*) 'Reflections',reflection,'Subsurface',subsurface,'Atmosphere',atmosphere
  
  call readdem(h)
  call difftopo(NSx,NSy,h,dx,dy,surfaceSlope,azFac)

  Tsurf=0.; Qrefl=0.; QIRre=0.  
  Qmax=0.; Qmeans=0.; Tmean=0.; Qrefm=0.
  nm=0   
  
  print *,'...reading horizons file...'
  call readhorizons('horizons.'//fileext)

  if (reflection) then
     print *,'...reading huge fieldofviews file...'
     call getmaxfieldsize(NSx,NSy,fileext,CCMAX)
     print *,'... max field of view size=',CCMAX
     allocate(ii(NSx,NSy,CCMAX), jj(NSx,NSy,CCMAX), dO12(NSx,NSy,CCMAX))
     call getfieldofview(NSx,NSy,fileext,cc,ii,jj,dO12,skysize,CCMAX)
  endif
  if (subsurface) allocate(T(NSx,NSy,nz), Qnm1(NSx,NSy))
  
  print *,'...calculating...'
  ! loop over time steps 
  do n=0,nsteps-1
     edays = (n+1)*dtmin/1440.
     udtTime%dMinutes = udtTime%dMinutes + dtmin
     call addtime(udtTime)
     if (mod(n,96*3)>=96) cycle  ! skip days for speedup
     print *,edays

     call sunpos(udtTime, longitude, latitude, dZenithAngle, dAzimuth, R)
     sinbeta = cos(dZenithAngle*d2r)
     azSun = dAzimuth*d2r

     do i=2,NSx-1
        do j=2,NSy-1
           emax = getonehorizon(i,j,azSun)
           Qn(i,j)=flux_wshad(R,sinbeta,azSun,surfaceSlope(i,j),azFac(i,j),emax)
        enddo
     enddo
     Qref=flux_wshad(R,sinbeta,azSun,zero,zero,zero)

     if (atmosphere) then
        unsd = mk_atmosphere(dZenithAngle*d2r,I0,D0)
        S0=1365./R**2  ! must be the same as in flux_wshad
        if (reflection) then
           forall(i=2:NSx-1, j=2:NSy-1)
              Qn(i,j)=Qn(i,j)*I0 + S0*D0*(1.-skysize(i,j)/(2*pi))
           end forall
        else
           Qn(:,:) = Qn(:,:)*I0 + S0*D0
        endif
        Qref = Qref*I0 + S0*D0
     endif
     
     if (reflection) then
        Qvis(:,:) = Qn + Qrefl
        QIRin(:,:) = QIR + QIRre
        do i=2,NSx-1
           do j=2,NSy-1
              QIR(i,j)=0.; Qrefl(i,j)=0.; QIRre(i,j)=0.
              do k=1,cc(i,j)
                 iii = ii(i,j,k); jjj = jj(i,j,k)
                 v = viewing_angle(i,j,iii,jjj,h)
                 Qrefl(i,j) = Qrefl(i,j) + dO12(i,j,k)/pi*cos(v)*albedo(iii,jjj)*Qvis(iii,jjj)
                 QIR(i,j) = QIR(i,j) + dO12(i,j,k)/pi*cos(v)*emiss*sigSB*Tsurf(iii,jjj)**4
                 QIRre(i,j) = QIRre(i,j) + dO12(i,j,k)/pi*cos(v)*(1-emiss)*QIRin(iii,jjj)
              enddo
           enddo
        enddo
        Qabs(:,:)=(1.-albedo)*(Qn+Qrefl)+emiss*(QIR+QIRre)  ! Q absorbed
        forall(i=2:NSx-1, j=2:NSy-1) !!!!
           Qabs(i,j)=Qabs(i,j)+Qother*(1.-skysize(i,j)/(2*pi))
        end forall
     else
        Qabs(:,:)=(1.-albedo)*Qn  + Qother
     endif
     
     if (subsurface) then
        if (n==0) then
           Qnm1 = Qabs
           Tsurf = -275. ! negative of initialization temperature
        endif
        do i=2,NSx-1
           do j=2,NSy-1
              call subsurfaceconduction(T(i,j,:),Tsurf(i,j),dtmin*60.,Qnm1(i,j),Qabs(i,j),emiss)
           enddo
        enddo
        Qnm1 = Qabs
     else
        Tsurf = (Qabs/sigSB/emiss)**0.25
     endif
     
     if (edays>tmax-365) then
     !if (edays>tmax-1.) then
        Qmeans(:,:,1) = Qmeans(:,:,1) + Qn
        Qmeans(:,:,2) = Qmeans(:,:,2) + Qabs
        where (Qn>Qmax) Qmax=Qn
        Qmeans(:,:,3) = Qmeans(:,:,3) + QIR
        Qmeans(:,:,4) = Qmeans(:,:,4) + Qrefl
        Tmean = Tmean + Tsurf
        Qrefm = Qrefm + Qref
        nm=nm+1
     endif

  enddo  ! end of time loop

  if (reflection) deallocate(ii, jj, dO12)
  if (subsurface) deallocate(T, Qnm1)

  Qmeans=Qmeans/nm
  Tmean=Tmean/nm
  
  print *,'Qref=',Qrefm/nm  ! contains neither albedo nor subsurface conduction
  open(unit=21,file='qmean.dat',status='unknown',action='write')
  do i=2,NSx-1
     do j=2,NSy-1
        write(21,'(2(i4,1x),f9.2,2x,f6.4,5(1x,f6.1),1x,f5.1)') &
             & i,j,h(i,j),surfaceSlope(i,j),Qmeans(i,j,1),Qmax(i,j),Qmeans(i,j,2:4),Tmean(i,j)
        !write(21,'(2(i4,1x),f9.2,2x,f6.4,5(1x,f6.1),1x,f5.1)') &  ! instanteneous values
        !     & i,j,h(i,j),surfaceSlope(i,j),Qn(i,j),Qmax(i,j), &
        !     & Qabs(i,j),Qir(i,j),Qrefl(i,j),Tsurf(i,j)
     enddo
  enddo
  close(21)

end program cratersQ_earth



subroutine subsurfaceconduction(T,Tsurf,dtsec,Qn,Qnp1,emiss)
  ! 1d subsurface conduction
  use allinterfaces, only : conductionQ
  implicit none
  integer, parameter :: NMAX=1000, Ni=5
  real(8), parameter :: pi=3.1415926535897932
  real(8), parameter :: solarDay = 86400., Fgeotherm = 0.065 
  integer, parameter :: nz=40

  real(8), intent(INOUT) :: T(NMAX), Tsurf
  real(8), intent(IN) :: dtsec,Qn,Qnp1,emiss
  integer i
  real(8) zmax, zfac, Fsurf, Tinit, delta
  real(8) Tsurfold, Told(1:nz)
  real(8), save :: ti(NMAX), rhocv(NMAX), z(NMAX)
  logical, save :: first = .true.

  if (first) then ! initialize grid
     ti(:) = 1000.;  rhocv(:) = 1200.*800.  ! adjust
     zmax=2.; zfac = 1.05  ! adjust

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

     first = .false.
  endif
  
  if (Tsurf<=0.) then  ! initialize temperature profile
     if (Tsurf==0.) Tinit=273.
     if (Tsurf<0.) Tinit=-Tsurf
     T(1:nz) = Tinit
     Tsurf = Tinit
  endif

  Tsurfold=Tsurf
  Told(1:nz)=T(1:nz)
  call conductionQ(nz,z,dtsec,Qn,Qnp1,T,ti,rhocv,emiss,Tsurf,Fgeotherm,Fsurf)
  
end subroutine subsurfaceconduction



