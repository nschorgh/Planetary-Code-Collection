program cratersQ_mars_shadowsonly
  ! map of direct solar irradiance
  use filemanager
  use allinterfaces
  use newhorizons
  implicit none
  real(8), parameter :: pi=3.1415926535897932, d2r=pi/180.
  real(8), parameter :: solsy = 668.60 ! solar days per Mars year
  real(8), parameter :: solarDay = 88775.244  ! Mars
  real(8) edays, marsR, marsLs, marsDec

  integer nsteps, n, i, j, nm
  real(8) tmax, dt, latitude, dtsec
  real(8) HA, sdays, azSun, emax, sinbeta
  real(8), allocatable, dimension(:,:) :: h, surfaceSlope, azFac
  real(8), allocatable, dimension(:,:) :: Qn, Qmean, Qmax
  real(8), allocatable, dimension(:,:) :: shadowtime, maxshadowtime

  integer julday, iyr, imm, iday
  real(8) jd, dt0_J2000
  real(8) jd_snap, longitude, buf, LTST
  
  integer, parameter :: Mx1=2, Mx2=NSx-1, My1=2, My2=NSy-1
  
  allocate(h(NSx,NSy), surfaceSlope(NSx,NSy), azFac(NSx,NSy))
  allocate(Qn(NSx,NSy), Qmean(NSx,NSy), Qmax(NSx,NSy))
  allocate(shadowtime(NSx,NSy), maxshadowtime(NSx,NSy))
  Qn=0.; Qmean=0; Qmax=0; shadowtime=0; maxshadowtime=0
  
  dt=0.02; 
  tmax = 2*solsy+1.
  latitude = -41.6

  ! set some constants
  nsteps=int(tmax/dt)       ! calculate total number of timesteps
  dtsec = dt*solarDay
  
  write(*,*) 'Time step=',dt,' Max number of steps=',nsteps
  write(*,*) 'Calculations performed for latitude=',latitude
  write(*,*) 'Nx=',NSx,'Ny=',NSy,'File=',fileext
  write(*,*) 'Region of interest: (',Mx1,',',My1,') x (',Mx2,',',My2,')'

  ! date at the beginning of run
  imm = 1; iday=1; iyr=2013
  jd=dble(julday(imm,iday,iyr))  !  JD for noon UTC on iyear/imm/iday
  call marsclock24(jd,dt0_J2000,marsLs,marsDec,marsR,0d0,LTST) ! calculate dt0_J2000
  !call marsorbit(dt0_j2000,0.d0,marsLs,marsDec,marsR)
  
  call readdem(h)
  call difftopo(NSx,NSy,h,dx,dy,surfaceSlope,azFac)

  latitude=latitude*d2r

  nm=0
  
  print *,'...reading horizons file...'
  call readhorizons

  ! image taken imm = 5; iday=15; iyr=2014 8:44 UTC  ESP_036561
  longitude = 360 - 202.3 ! west longitude
  imm = 5; iday=15; iyr=2014 
  jd_snap=dble(julday(imm,iday,iyr)) + (8.+44./60-12)/24.  ! noon UTC
  call marsclock24(jd_snap,buf,marsLs,marsDec,marsR,Longitude,LTST)

  open(unit=25,file='qsnap.dat',status='unknown',action='write') ! snapshot
  
  print *,'...calculating...'
  ! loop over time steps 
  do n=0,nsteps-1
     sdays = (n+1)*dtsec/solarDay
     edays = (n+1)*dtsec/86400.
     
     !call marsorbit(dt0_j2000,edays,marsLs,marsDec,marsR)
     !HA=2.*pi*mod(sdays,1.d0)   ! hour angle
     call marsclock24(jd+edays,buf,marsLs,marsDec,marsR,Longitude,LTST)
     HA=2.*pi*mod(LTST+12,24.d0)/24
     
     call equatorial2horizontal(marsDec,latitude,HA,sinbeta,azSun)
     
     if (mod(n,10)==0) print *,n,sdays,mod(marsLs/d2r,360.d0)
     
     do i=Mx1,Mx2
        do j=My1,My2
           if (h(i,j)<-32000.) cycle
           emax = getonehorizon(i,j,azSun)
           Qn(i,j)=flux_wshad(marsR,sinbeta,azSun,surfaceSlope(i,j),azFac(i,j),emax)
        enddo
     enddo

     ! only diagnostics below this line
     if (sdays > tmax-nint(solsy)) then
        Qmean(:,:) = Qmean(:,:) + Qn
        where (Qn>Qmax) Qmax=Qn
        nm=nm+1
     endif
     if (sdays > tmax-2*solsy) then  ! longest continuous period
        where (Qn==0.) 
           shadowtime=shadowtime+dt
        elsewhere
           shadowtime=0.
        end where
        where (shadowtime>maxshadowtime) maxshadowtime=shadowtime
     endif

     if (jd+edays > jd_snap-dt/2 .and. jd+edays <= jd_snap+dt/2) then
        print *,'writing snapshot'
        print *,'current hour angle',HA
        print *,'current altitude of sun',asin(sinbeta)/d2r
        print *,'current azimuth of sun',azSun/d2r
        do i=Mx1,Mx2
           do j=My1,My2
              write(25,'(2(i4,1x),f9.2,1x,f6.1)') i,j,h(i,j),Qn(i,j)
           enddo
        enddo
        close(25)
     endif
  
  enddo  ! end of time loop

  Qmean = Qmean/nm

  open(unit=21,file='qshadow.dat',status='unknown',action='write')
  do i=Mx1,Mx2
     do j=My1,My2
        write(21,'(2(i4,1x),f9.2,2(1x,f6.1),1x,f6.1)') &
             & i,j,h(i,j),Qmean(i,j),Qmax(i,j),maxshadowtime(i,j)
     enddo
  enddo
  close(21)
  
end program cratersQ_mars_shadowsonly
