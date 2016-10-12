program cratersQ_earth_static
!***********************************************************************
!   program to calculate insolation with shadowing and reflections  
!             with Mauna Kea atmosphere
!             static solution with mean insolation
!
!   written by Norbert Schorghofer 2016
!***********************************************************************
  use filemanager
  use allinterfaces
  use dateformat
  implicit none
  real(8), parameter :: pi=3.1415926535897932, d2r=pi/180.
  real(8), parameter :: earthDay=86400.
  real(8), parameter :: sigSB = 5.6704e-8  

  integer nsteps, n, i, j, nm, k, CCMAX, iii, jjj
  real(8) tmax, edays, dtmin, latitude
  real(8) R, dZenithAngle, dAzimuth, longitude
  real(8) azSun, sinbeta, smax, emiss, v
  real(8), dimension(NSx,NSy) :: h, surfaceSlope, azFac
  real(8), dimension(NSx,NSy) :: Qn, QIR, Qrefl   ! incoming
  integer, dimension(NSx,NSy) :: cc
  integer(2), dimension(:,:,:), allocatable :: ii,jj
  real(4), dimension(:,:,:), allocatable :: dO12
  real(8), dimension(NSx,NSy) :: Tsurf, Qvis, skysize, Qabs, albedo, QIRin, QIRre, Qabsold
  real(8) Qmeans(NSx,NSy,4), Qmax(NSx,NSy), Tmean(NSx,NSy), Qref, Qrefm, Tsurfm
  real(8) I0,D0,S0,unsd  ! atmosphere
  type(cTime) udtTime
  real(8), parameter :: zero=0.
  logical, parameter :: reflection=.true., atmosphere=.true.
  real(8), parameter :: Qother = 0.  ! (W/m^2) background (LW+downwelling)
  
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
  write(*,*) 'Reflections',reflection,'Atmosphere',atmosphere
  
  call readdem(h)
  call difftopo(NSx,NSy,h,dx,dy,surfaceSlope,azFac)

  Qmax=0.; Qmeans=0.; Tmean=0.; Qrefm=0.
  nm=0   
  
  print *,'...reading horizons file...'
  call gethorizon(0,0,azSun,smax,.TRUE.)

  print *,'...calculating...'
  ! loop over time steps 
  do n=0,nsteps-1
     edays = (n+1)*dtmin/1440.
     print *,edays
     udtTime%dMinutes = udtTime%dMinutes + dtmin
     call addtime(udtTime)

     call sunpos(udtTime, longitude, latitude, dZenithAngle, dAzimuth, R)
     sinbeta = cos(dZenithAngle*d2r)
     azSun = dAzimuth*d2r

     do i=2,NSx-1
        do j=2,NSy-1
           call gethorizon(i,j,azSun,smax,.FALSE.)
           Qn(i,j)=flux_wshad(R,sinbeta,azSun,surfaceSlope(i,j),azFac(i,j),smax)
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
     Qabs(:,:)=(1.-albedo)*Qn
     
     Tsurf = (Qabs/sigSB/emiss)**0.25
     
     if (edays>tmax-365) then
     !if (edays>tmax-1.) then
        Qmeans(:,:,1) = Qmeans(:,:,1) + Qn
        Qmeans(:,:,2) = Qmeans(:,:,2) + Qabs
        where (Qn>Qmax) Qmax=Qn
        Tmean = Tmean + Tsurf
        Qrefm = Qrefm + Qref
        nm=nm+1
     endif

  enddo  ! end of time loop

  Qmeans=Qmeans/nm
  Tmean=Tmean/nm
  Qrefm = Qrefm/nm

  Tsurfm = (((1.-albedo(1,1))*Qrefm+Qother)/sigSB/emiss)**0.25
  print *,'Qref=',Qrefm,Tsurfm  ! contains no albedo

  ! static part

  Tsurf=Tmean; Qn=Qmeans(:,:,1)
  Qrefl=0.; QIRre=0.  
  Tmean=0.
  nm=0
  
  if (reflection) then
     print *,'...reading huge fieldofviews file...'
     call getmaxfieldsize(NSx,NSy,fileext,CCMAX,1)
     print *,'... max field of view size=',CCMAX
     allocate(ii(NSx,NSy,CCMAX), jj(NSx,NSy,CCMAX), dO12(NSx,NSy,CCMAX))
     call getfieldofview(NSx,NSy,fileext,cc,ii,jj,dO12,skysize,CCMAX)

     Qabs = Qmeans(:,:,2);

     do n=1,5 ! loop over relaxation iteration
        
        print *,'...static iteration',n
        Qabsold = Qabs;
        
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
        forall(i=2:NSx-1, j=2:NSy-1)
           Qabs(i,j)=Qabs(i,j)+Qother*(1.-skysize(i,j)/(2*pi))
        end forall
        
        Tsurf = (Qabs/sigSB/emiss)**0.25
        print *,'Difference',sum(abs(Qabs-Qabsold))/size(Qabs)
     enddo  ! end of iteration loop
     
     deallocate(ii, jj, dO12)
  else ! no reflection
     Qabs(:,:)=(1.-albedo)*Qn  + Qother
     Tsurf = (Qabs/sigSB/emiss)**0.25
  endif
  
  Qmeans(:,:,2) = Qabs
  Qmeans(:,:,3) = QIR
  Qmeans(:,:,4) = Qrefl
  Tmean = Tsurf
  
  open(unit=21,file='qmean.dat',status='unknown',action='write')
  do i=2,NSx-1
     do j=2,NSy-1
        write(21,'(2(i4,1x),f9.2,2x,f6.4,5(1x,f6.1),1x,f5.1)') &
             & i,j,h(i,j),surfaceSlope(i,j),Qmeans(i,j,1),Qmax(i,j),Qmeans(i,j,2:4),Tmean(i,j)
     enddo
  enddo
  close(21)

end program cratersQ_earth_static

