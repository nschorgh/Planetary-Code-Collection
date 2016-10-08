program insol_earth
!***********************************************************************
!   insol_earth: program to calculate direct insolation with shadowing, 
!                Earth orbit, and Mauna Kea atmosphere
!
!   written by Norbert Schorghofer 2010-2016
!***********************************************************************
  use filemanager
  use allinterfaces
  use dateformat
  implicit none
  real(8), parameter :: pi=3.1415926535897932, d2r=pi/180.

  integer nsteps, n, i, j, nm
  real(8) tmax, edays, dtmin, latitude
  real(8) R, dZenithAngle, dAzimuth, longitude
  real(8) azSun, sinbeta, smax
  real(8), dimension(NSx,NSy) :: h, surfaceSlope, azFac
  real(8), dimension(NSx,NSy) :: Qn   ! incoming
  real(8), dimension(NSx,NSy) :: skysize
  real(8) Qmeans(NSx,NSy), Qmax(NSx,NSy), Qref, Qrefm
  real(8) I0, D0, S0, unsd  ! atmosphere
  type(cTime) udtTime
  real(8), parameter :: zero=0.
  logical, parameter :: atmosphere=.true.

  !dtmin=15. 
  dtmin=10.
  tmax = 365.+1
  !tmax = 60.
  latitude = 19.821; longitude = -155.468  ! Mauna Kea summit
  ! azimuth in degrees east of north, 0=north facing

  ! set some constants
  nsteps=int(tmax*1440./dtmin)       ! calculate total number of timesteps

  ! start time in UTC = HST-10
  udtTime = cTime(2013,1,1,0.,0.,0.) ! 14:00 HST 
  !udtTime = cTime(2012,5,20,0.,0.,0.) ! 14:00 HST 
  
  write(*,*) 'Starting time',udtTime
  write(*,*) 'Time step=',dtmin,'(min)  Max number of steps=',nsteps
  write(*,*) 'Calculations performed for latitude=',latitude
  write(*,*) 'Nx=',NSx,'Ny=',NSy,'File=',fileext
  write(*,*) 'Atmosphere',atmosphere
  
  call readdem(h)
  call difftopo(NSx,NSy,h,dx,dy,surfaceSlope,azFac)

  Qmax=0.; Qmeans=0.; Qrefm=0.
  nm=0   
  
  if (atmosphere) call getskysize(skysize) 

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
        Qn(:,:) = Qn(:,:)*I0 + S0*D0*skysize(:,:)/(2*pi)
        Qref = Qref*I0 + S0*D0
     endif

     if (edays>tmax-365) then
     !if (edays>tmax-1.) then
     !if (udtTime%iMonth==6) then
        Qmeans = Qmeans + Qn
        where (Qn>Qmax) Qmax=Qn
        Qrefm = Qrefm + Qref
        nm=nm+1
     endif

  enddo  ! end of time loop

  Qmeans=Qmeans/nm
  
  print *,'Qref=',Qrefm/nm  
  open(unit=21,file='qmean.dat',status='unknown',action='write')
  do i=2,NSx-1
     do j=2,NSy-1
        write(21,'(2(i4,1x),f9.2,2x,f6.4,2(1x,f6.1))') &
             & i,j,h(i,j),surfaceSlope(i,j),Qmeans(i,j),Qmax(i,j)
        !write(21,'(2(i4,1x),f9.2,2x,f6.4,5(1x,f6.1),1x,f5.1)') &  ! instanteneous values
        !     & i,j,h(i,j),surfaceSlope(i,j),Qn(i,j),Qmax(i,j)
     enddo
  enddo
  close(21)

end program insol_earth



subroutine getskysize(skysize)
!***********************************************************************
!   reads horizons file and calculates sky size
!***********************************************************************
  use filemanager, only : NSx,NSy,fileext
  use allinterfaces
  implicit none
  real(8), parameter :: pi=3.1415926535897932, d2r=pi/180.
  integer, parameter :: nres=360   ! # of azimuths
  real(8) smax(nres)
  integer i, j, ii, jj, ierr
  real(8), intent(OUT) :: skysize(NSx,NSy) 

  ! azimuth in degrees east of north, 0=north facing, 0...2*pi

  print *,'# azimuth rays = ',nres
  write(*,*) 'Nx=',NSx,'Ny=',NSy,'File=',fileext
  
  print *,'...reading horizons file ...'
  open(unit=21,file='horizons.'//fileext,status='old',action='read',iostat=ierr)
  if (ierr>0) stop 'skysize: Input file not found'
  
  do i=2,NSx-1
     do j=2,NSy-1
        read(21,*) ii,jj,smax(:)
        if (ii/=i .or. jj/=j) stop 'index mismatch'
        skysize(i,j) = sum(atan(smax))*2*pi/nres
        !print *,i,j,skysize(i,j),atan(maxval(smax))/d2r
     enddo
  enddo
  skysize = 2*pi - skysize

  close(21)
end subroutine getskysize
