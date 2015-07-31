program cratersQ_Moon
!****************************************************************************
!   cratersQ: program to calculate insolation with shadowing and reflections
!
!             zero thermal inertia, Earth orbit
!  
!   written by Norbert Schorghofer 2010-2015  
!****************************************************************************
  !use omp_lib
  use filemanager
  use allinterfaces
  implicit none
  real(8), parameter :: pi=3.1415926535897932, d2r=pi/180.
  real(8), parameter :: sigSB = 5.6704e-8  

  integer nsteps, n, i, j, nm, k, CCMAX, iii, jjj, nm2
  real(8) tmax, dt, latitude, dtsec
  real(8) R, Decl, HA, sdays
  real(8) azSun, sinbeta, smax, emiss, v
  real(8), dimension(NSx,NSy) :: h, surfaceSlope, azFac
  real(8), dimension(NSx,NSy) :: Qn, QIR, Qrefl   ! incoming
  integer, dimension(NSx,NSy) :: cc
  integer(2), dimension(:,:,:), allocatable :: ii,jj
  real(4), dimension(:,:,:), allocatable :: dO12  
  real(8), dimension(NSx,NSy) :: Tsurf, Qvis, skysize, Qabs, albedo, QIRin, QIRre
  real(8) Qmeans(NSx,NSy,4)
  real(8), dimension(NSx,NSy) :: Qmax, Tmean, Tmaxi, Tbottom, Tb, Tmaxi2
  real(8), allocatable :: T(:,:,:), Qnm1(:,:)  ! subsurface
  logical, parameter :: reflection=.true., subsurface=.false.
  
  dt=0.01; 
  tmax = 2.
  latitude = 87.
  ! azimuth in degrees east of north, 0=north facing
  albedo(:,:) = 0.12
  emiss = 1.

  ! set some constants
  nsteps=int(tmax/dt)       ! calculate total number of timesteps
  dtsec = dt*solarDay
  
  write(*,*) 'Time step=',dt,' Max number of steps=',nsteps
  write(*,*) 'Calculations performed for latitude=',latitude
  write(*,*) 'Nx=',NSx,'Ny=',NSy,'File=',fileext
  write(*,*) 'Mean albedo=',sum(albedo)/size(albedo),'Emissivity=',emiss
  write(*,*) 'Reflections:',reflection,'Subsurface:',subsurface

  ! setenv OMP_NUM_THREADS 4
  !write (*,'(a,i8)') 'The number of processors available = ', omp_get_num_procs()
  !write (*,'(a,i8)') 'The number of threads available    = ', omp_get_max_threads()

  call readdem(NSx,NSy,h,fileext)
  call difftopo(NSx,NSy,h,dx,dy,surfaceSlope,azFac)

  latitude=latitude*d2r
  Tsurf=0.; Qrefl=0.; QIRre=0.  
  Qmeans=0.; Tmean=0.; Tbottom=0.; nm=0
  Qmax=0.; Tmaxi=0.; Tmaxi2=0.
  Tb=0.; nm2=0
  
  print *,'...reading horizons file...'
  call gethorizon(0,0,azSun,smax,.TRUE.)

  if (reflection) then
     print *,'...reading huge fieldofviews file...'
     call getmaxfieldsize(NSx,NSy,fileext,CCMAX,1)
     print *,'... max field of view size=',CCMAX
     allocate(ii(NSx,NSy,CCMAX), jj(NSx,NSy,CCMAX), dO12(NSx,NSy,CCMAX))
     call getfieldofview(NSx,NSy,fileext,cc,ii,jj,dO12,skysize,CCMAX)
  endif
  if (subsurface) allocate(T(NSx,NSy,1000), Qnm1(NSx,NSy))

  print *,'...calculating...'
  ! loop over time steps 
  do n=0,nsteps-1
     sdays = (n+1)*dtsec/solarDay
     !print *,sdays
     !print *,sdays,Tsurf(20,20)
     Decl=0.; R=1.

     ! for small landscapes these can be here
     HA=2.*pi*mod(sdays,1.)   ! hour angle
     call equatorial2horizontal(Decl,latitude,HA,sinbeta,azSun)

     print *,sdays,HA,n
     
     do i=2,NSx-1
        do j=2,NSy-1
           call gethorizon(i,j,azSun,smax,.FALSE.)
           Qn(i,j)=flux_wgeom(R,sinbeta,azSun,surfaceSlope(i,j),azFac(i,j),smax)
        enddo
     enddo

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
        Qabs(:,:)=(1.-albedo(:,:))*(Qn+Qrefl)+emiss*(QIR+QIRre)  ! Q absorbed
     else
        Qabs(:,:)=(1.-albedo(:,:))*Qn
     endif

     if (subsurface) then
        if (n==0) then ! initialization
           Qnm1 = Qabs
           !Tsurf = 0.
           Tsurf = -200. ! negative of initialization temperature
        endif

        ! for faster convergence
        if (n>=10 .and. n<10+nint(1/dt)) then
           Tb = Tb+Tsurf
           nm2 = nm2+1
        endif
        if (n==10+nint(1/dt)) then
           Tb = Tb/nm2
           Tsurf = -Tb  ! re-initialize
        endif

        ! main part
        do i=2,NSx-1
           do j=2,NSy-1
              call subsurfaceconduction(T(i,j,:),Tsurf(i,j),dtsec,Qnm1(i,j),Qabs(i,j),emiss)
           enddo
        enddo
        Qnm1 = Qabs
     else
        Tsurf = (Qabs/sigSB)**0.25
     endif
     
     if (sdays > tmax-1.) then
        Qmeans(:,:,1) = Qmeans(:,:,1) + Qn
        Qmeans(:,:,2) = Qmeans(:,:,2) + Qabs
        where (Qabs>Qmax) Qmax=Qabs
        Qmeans(:,:,3) = Qmeans(:,:,3) + QIR
        Qmeans(:,:,4) = Qmeans(:,:,4) + Qrefl
        Tmean = Tmean + Tsurf
        where (Tsurf>Tmaxi)
           Tmaxi2=Tmaxi
           Tmaxi=Tsurf
        elsewhere (Tsurf>Tmaxi2)
           Tmaxi2=Tsurf
        end where
        if (subsurface) Tbottom=Tbottom+T(:,:,nz)
        nm=nm+1
     endif

  enddo  ! end of time loop

  if (reflection) deallocate(ii, jj, dO12)
  if (subsurface) deallocate(T, Qnm1)
  
  Qmeans=Qmeans/nm
  Tmean=Tmean/nm; Tbottom=Tbottom/nm

  open(unit=21,file='qmean.dat',status='unknown',action='write')
  do i=2,NSx-1
     do j=2,NSy-1
        write(21,'(2(i4,1x),f9.2,2x,f6.4,5(1x,f6.1),4(1x,f5.1))') &
             & i,j,h(i,j),surfaceSlope(i,j),Qmeans(i,j,1),Qmax(i,j),Qmeans(i,j,2:4), &
             & Tmean(i,j),Tmaxi(i,j),Tmaxi2(i,j),Tbottom(i,j)
        write(22,'(2(i4,1x),f9.2,2x,f6.4,5(1x,f6.1),1x,f5.1)') &  ! instanteneous values
             & i,j,h(i,j),surfaceSlope(i,j),Qn(i,j),Qmax(i,j), &
             & Qabs(i,j),Qir(i,j),Qrefl(i,j),Tsurf(i,j)
     enddo
  enddo
  close(21)

end program cratersQ_Moon
 
