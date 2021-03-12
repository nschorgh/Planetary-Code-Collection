program insol3d_earth
!***********************************************************************
!   insol_earth3d: program to calculate direct insolation with terrain 
!                  shadowing, Earth orbit, and Mauna Kea atmosphere
!
!   written by Norbert Schorghofer 2010-2016
!***********************************************************************
  use filemanager
  use allinterfaces
  use dateformat
  use newhorizons
  implicit none
  real(8), parameter :: pi=3.1415926535897932, d2r=pi/180.
  real(8), parameter :: zero=0.

  integer nsteps, n, i, j, nm
  real(8) tmax, edays, dtmin, latitude
  real(8) R, dZenithAngle, dAzimuth, longitude
  real(8) azSun, sinbeta, emax
  real(8), dimension(:,:), allocatable :: h, SlopeAngle, azFac
  real(8), dimension(:,:), allocatable :: Qn, Qsw   ! incoming
  real(8), dimension(:,:), allocatable :: Qmeans, Qmax, daytime
  real(8), dimension(NSx,NSy) :: ditime=0., longday=-999., shortday=1e32
  real(8) Qflat, Qflatm, alltime
  real(8) I0, D0, S0, unsd  ! atmosphere
  real(8) :: oldHours = -999.
  type(cTime) udtTime
  logical, parameter :: atmosphere = .true.

  allocate(h(NSx,NSy), SlopeAngle(NSx,NSy), azFac(NSx,NSy))
  allocate(Qn(NSx,NSy), Qsw(NSx,NSy))
  allocate(Qmeans(NSx,NSy), Qmax(NSx,NSy), daytime(NSx,NSy), source=0.d0)

  dtmin=10.
  ! dtmin=5.
  tmax = 365.+1
  
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
  call difftopo(NSx,NSy,h,dx,dy,SlopeAngle,azFac)

  nm=0; Qflatm=0.; alltime=0.
  
  print *,'...reading horizons file...'
  call readhorizons

  print *,'...calculating...'
  ! loop over time steps 
  do n=0,nsteps-1
     edays = (n+1)*dtmin/1440.
     print *,edays
     udtTime%dMinutes = udtTime%dMinutes + dtmin
     oldHours = udtTime%dHours
     call addtime(udtTime)

     call sunpos(udtTime, longitude, latitude, dZenithAngle, dAzimuth, R)
     sinbeta = cos(dZenithAngle*d2r)
     azSun = dAzimuth*d2r

     do i=2,NSx-1
        do j=2,NSy-1
           emax = getonehorizon(i,j,azSun)
           Qn(i,j) = flux_wshad(R,sinbeta,azSun,SlopeAngle(i,j),azFac(i,j),emax)
        end do
     end do
     Qflat = flux_wshad(R,sinbeta,azSun,zero,zero,zero)
     
     if (atmosphere) then
        unsd = mk_atmosphere(dZenithAngle*d2r,I0,D0)
        S0 = 1365./R**2  ! must be the same as in flux_wshad
        Qsw(:,:) = Qn(:,:)*I0 + S0*D0  ! do not use skysize
        Qflat = Qflat*I0 + S0*D0
     else
        Qsw = Qn
     endif
     
     if (edays>tmax-365.) then  ! one year
     !if (edays>tmax-1.) then  ! one day
     !if (udtTime%iMonth==6) then  ! June
        Qmeans = Qmeans + Qsw
        where (Qsw>Qmax) Qmax=Qsw
        Qflatm = Qflatm + Qflat
        nm = nm+1
        
        where (Qn>0.) daytime = daytime + dtmin
        alltime = alltime + dtmin
     endif

     ! length of day
     if (udtTime%dHours < oldHours) then ! new day
        if (edays>tmax-365.) then
           where (ditime > longday) longday=ditime
           where (ditime < shortday) shortday=ditime
        endif
        ditime = 0.
     endif
     where (Qn>0) ditime = ditime + dtmin
     
  enddo  ! end of time loop

  Qmeans = Qmeans/nm
  longday = longday/60.  ! min -> hr
  shortday = shortday/60.  ! min -> hr
  daytime = daytime/alltime*24
  
  print *,'Qflat=',Qflatm/nm
  open(unit=21,file='qmean.dat',status='unknown',action='write')
  do i=2,NSx-1
     do j=2,NSy-1
        write(21,'(2(i4,1x),f9.2,2x,f6.4,2(1x,f6.1),1x,f6.3)') &
             & i,j,h(i,j),SlopeAngle(i,j),Qmeans(i,j),Qmax(i,j),daytime(i,j)
        !write(21,'(2(i4,1x),f9.2,2x,f6.4,5(1x,f6.1),1x,f5.1)') &
        !     & i,j,h(i,j),SlopeAngle(i,j),Qn(i,j),Qmax(i,j)
        !write(22,'(2(i4,1x),f9.2,2(1x,f6.1),3(1x,f6.3))') i,j,h(i,j), &
        !     & -999.,-999.,daytime(i,j),shortday(i,j),longday(i,j)
     enddo
  enddo
  close(21)

end program insol3d_earth

