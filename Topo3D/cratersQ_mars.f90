!***************************************************************************
! cratersQ_mars: Mars thermal model with direct insolation, subsurface
!                conduction, terrain shadowing, and approximate self-heating
!***************************************************************************


module miscparams
  ! parameters that never change
  real(8), parameter :: pi=3.1415926535897932, d2r=pi/180.
  real(8), parameter :: sigSB = 5.6704e-8
  real(8), parameter :: Lco2frost=6.0e5  ! [J/kg]
  real(8), parameter :: zero = 0.
  real(8), parameter :: earthDay = 86400. ! [s]
  real(8), parameter :: solsy = 668.60 ! solar days per Mars year
  real(8), parameter :: solarDay = 88775.244  ! Mars [s]

  ! thermal model parameters
  real(8), parameter :: Tco2frost=145. ! adjust according to elevation [K]
  real(8), parameter :: Tfrost = 200.  ! H2O frost point temperature, for diagnostics only [K]
  real(8), parameter :: fracIR=0.04, fracDust=0.02
  real(8), parameter :: emiss = 0.98
  real(8), parameter :: Fgeotherm = 0.0 ! [W/m^2]
  integer, parameter :: nz=70
  real(8), parameter :: thIn = 600.  ! Thermal inertia
end module miscparams


program cratersQ_mars
  use filemanager
  use allinterfaces
  use miscparams
  use newhorizons
  implicit none

  real(8), parameter :: albedo0=0.12, co2albedo=0.65
  real(8) :: latitude = -41.6  ! [degree] Palikir Crater
  real(8), parameter :: longitude = 360 - 202.3 ! west longitude [degree]

  integer nsteps, n, i, j, nm
  real(8) tmax, dt, dtsec, buf
  real(8) HA, sdays, azSun, emax, sinbeta
  real(8) edays, marsR, marsLs, marsDec
  real(8), allocatable, dimension(:,:) :: h, surfaceSlope, azFac
  real(8), allocatable, dimension(:,:) :: Qn   ! incoming
  real(8), allocatable, dimension(:,:) :: Tsurf, albedo, m
  real(8), allocatable, dimension(:,:) :: Qmean, Qmax, Tmean, Tmaxi, Qdirect
  real(8), allocatable, dimension(:,:) :: skyview, gterm
  real(8), allocatable, dimension(:,:) :: mmax, frosttime, maxfrosttime, Qnm1
  real(8), allocatable :: Fsurf(:,:), T(:,:,:), Tbottom(:,:), Tref(:)  ! subsurface
  real(8), allocatable, dimension(:,:) :: mmin, h2olast
  real(8) dE, Tsurfold, QIR, Qrefl, Qscat, Qlw
  logical, parameter :: subsurface=.true.  ! control panel
  real(8) jd, LTST, jd_end
  
  integer k
  real(8) jd_snap(3), jd_themis(2)  ! Julian dates of snapshots
  character(len=20) fns(3), fnt(2)  ! file names of snapshots

  integer, parameter :: Mx1=2, Mx2=NSx-1, My1=2, My2=NSy-1
  
  allocate(h(NSx,NSy), surfaceSlope(NSx,NSy), azFac(NSx,NSy))
  allocate(Qn(NSx,NSy), Qnm1(NSx,NSy), Qdirect(NSx,NSy))
  allocate(Tsurf(NSx,NSy), m(NSx,NSy))
  allocate(albedo(NSx,NSy)); albedo = albedo0
  allocate(skyview(NSx,NSy), gterm(NSx,NSy))
  allocate(Qmean(NSx,NSy), Qmax(NSx,NSy), Tmean(NSx,NSy), Tmaxi(NSx,NSy))
  allocate(frosttime(NSx,NSy), maxfrosttime(NSx,NSy))
  allocate(mmax(NSx,NSy), mmin(NSx,NSy), h2olast(NSx,NSy))
  Qn=0.; Qnm1=0.; m=0.; Qdirect=0.; skyview=0.
  Qmean=0.; Qmax=0.; Tmean=0.; Tmaxi=0.
  frosttime=0.; maxfrosttime=0.; mmax=0.; mmin=1e32
  h2olast=-9.
  
  dt=0.02
  tmax = 6.*solsy

  ! set some constants
  nsteps=int(tmax/dt)       ! calculate total number of timesteps
  dtsec = dt*solarDay
  
  write(*,*) 'Time step=',dt,' Max number of steps=',nsteps
  write(*,*) 'Calculations performed for latitude=',latitude
  write(*,*) 'fracIR=',fracIR,'fracDust=',fracDust
  write(*,*) 'Nx=',NSx,'Ny=',NSy,'File=',fileext
  write(*,*) 'Region of interest: (',Mx1,',',My1,') x (',Mx2,',',My2,')'
  write(*,*) 'Mean albedo=',sum(albedo)/size(albedo),'Emissivity=',emiss
  write(*,*) 'CO2 frost temperature=',Tco2frost
  write(*,*) 'Reflections:',.FALSE.,'Subsurface:',subsurface

  ! Set start date
  !jd=dble(julday(1,1,2009))  !  JD for noon UTC on imm,iday,iyear

  ! Alternatively set end date
  jd_end=dble(julday(2,28,2017)) 
  jd = nint(jd_end - tmax*solarDay/earthDay)
  
  write(*,*) 'Julian start and end date',jd,jd_end
  
  call readdem(h)
  call difftopo(NSx,NSy,h,dx,dy,surfaceSlope,azFac)

  latitude=latitude*d2r
  Tsurf=200.
  nm=0

  print *,'...reading horizons file...'
  call readhorizons
  do concurrent(i=2:NSx-1, j=2:Nsy-1)
     skyview(i,j) = getoneskysize_v2(i,j)/(2*pi)
     gterm(i,j) = getoneGterm(i,j,surfaceSlope(i,j),azFac(i,j))
  end do
  
  if (subsurface) then ! initialize subsurface component
     allocate(T(nz,Mx1:Mx2,My1:My2))
     allocate(Tref(nz))
     allocate(Fsurf(NSx,NSy))
     call subsurfaceconduction_mars(Tref(:),buf,dtsec,zero,zero,buf,buf,.true.)
     allocate(Tbottom(NSx,NSy))
     Tbottom(:,:)=-9
     Tsurf(:,:)=-9  ! max 3 digits
     Fsurf(:,:)=0.
  end if

  open(unit=22,file='timeseries_flat.dat',status='unknown',action='write')

  ! image taken imm = 5; iday=15; iyr=2014 8:44 UTC  ESP_036561
  jd_snap(1)=dble(julday(5,15,2014)) + (8.+44./60-12)/24.  !  JD for noon UTC on imm,iday,iyear
  !call marsclock24(jd_snap,buf,marsLs,marsDec,marsR,Longitude,LTST)
  fns(1)='qsnap_036561.dat' ! snapshot
  
  jd_snap(2)=dble(julday(8,18,2012)) + (8.+48./60-12)/24.; fns(2)='qsnap_028412.dat'
  jd_snap(3)=dble(julday(5,11,2012)) + (17.+31./60-12)/24.; fns(3)='qsnap_027146.dat'

  jd_themis(1)=dble(julday(10,30,2016)) + (1.+30./60-12)/24.; fnt(1)='qsnap_i65997002.dat'
  jd_themis(2)=dble(julday(2,26,2017)) + (18.+03./60-12)/24.; fnt(2)='qsnap_i67450002.dat' 
  
  print *,'...calculating...'
  ! loop over time steps 
  do n=0,nsteps-1
     sdays = (n+1)*dtsec/solarDay
     edays = (n+1)*dtsec/earthDay

     call marsclock24(jd+edays,buf,marsLs,marsDec,marsR,Longitude,LTST)
     HA=2.*pi*mod(LTST+12,24.d0)/24  ! hour angle
     
     call equatorial2horizontal(marsDec,latitude,HA,sinbeta,azSun)
     
     if (mod(n,10)==0) print *,n,sdays,marsLs/d2r

     call flux_mars2(marsR,marsDec,latitude,HA,fracIR,fracDust, &
          & zero,zero,zero,Qdirect(1,1),Qscat,Qlw)
     Qn(1,1) = (1-albedo(1,1))*(Qdirect(1,1)+Qscat) + emiss*Qlw
     do i=Mx1,Mx2
        do j=My1,My2
           ! incoming solar flux
           if (h(i,j)<-32000) cycle
           emax = getonehorizon(i,j,azSun)
           call flux_mars2(marsR,marsDec,latitude,HA,fracIR,fracDust, &
                & surfaceSlope(i,j),azFac(i,j),emax,Qdirect(i,j),Qscat,Qlw)
           ! absorbed direct insolation and contributions from atmosphere
           Qn(i,j) = (1-albedo(i,j))*(Qdirect(i,j)+Qscat*skyview(i,j)) &
                & + emiss*Qlw*skyview(i,j)
           ! contribution from land in field of view
           if (n>0) then 
              QIR = gterm(i,j)*emiss*sigSB*Tsurf(1,1)**4
              Qn(i,j) = Qn(i,j) + emiss*QIR
           endif
           Qrefl = gterm(i,j)*albedo(1,1)*(Qdirect(1,1)+Qscat)
           Qn(i,j) = Qn(i,j) + (1-albedo(i,j))*Qrefl
        enddo
     enddo
     if (n==0) Qnm1(:,:) = Qn(:,:)
     
     if (subsurface) then ! subsurface conduction
        do i=Mx1,Mx2
           do j=My1,My2
              if (h(i,j)<-32000) cycle
              call subsurfaceconduction_mars(T(:,i,j),Tsurf(i,j), &
                   & dtsec,Qnm1(i,j),Qn(i,j),m(i,j),Fsurf(i,j),.false.)
           enddo
        enddo
        call subsurfaceconduction_mars(Tref(:),Tsurf(1,1), &
             & dtsec,Qnm1(1,1),Qn(1,1),m(1,1),Fsurf(1,1),.false.)

     else  ! no subsurface conduction
        do i=Mx1,Mx2
           do j=My1,My2
              if (h(i,j)<-32000) cycle
              Tsurf(i,j) = (Qn(i,j)/emiss/sigSB)**0.25
              Tsurfold = (Qnm1(i,j)/emiss/sigSB)**0.25
              if (Tsurf(i,j)<Tco2frost.or.m(i,j)>0.) then   ! CO2 condensation
                 Tsurf(i,j)=Tco2frost
                 dE = - Qn(i,j) + emiss*sigSB*(Tsurf(i,j)**4 + Tsurfold**4)/2.
                 m(i,j) = m(i,j) + dtsec*dE/Lco2frost
              endif
           enddo
        enddo
        Tsurf(1,1) = (Qn(1,1)/emiss/sigSB)**0.25
        Tsurfold = (Qnm1(1,1)/emiss/sigSB)**0.25
        if (Tsurf(1,1)<Tco2frost.or.m(1,1)>0.) then   ! CO2 condensation
           Tsurf(1,1)=Tco2frost
           dE = - Qn(1,1) + emiss*sigSB*(Tsurf(1,1)**4 + Tsurfold**4)/2.
           m(1,1) = m(1,1) + dtsec*dE/Lco2frost
        endif
     endif
     Qnm1(:,:) = Qn(:,:)

     where (Tsurf>Tco2frost .or. m<=0.)
        albedo = albedo0
     elsewhere
        albedo = co2albedo
     end where

     ! only diagnostics below this line
     if (sdays > tmax-solsy) then
        Qmean(:,:) = Qmean(:,:) + Qn
        where (Qn>Qmax) Qmax=Qn
        Tmean = Tmean + Tsurf
        where (Tsurf>Tmaxi) Tmaxi=Tsurf
        where (m>mmax) mmax=m
        where (m<mmin) mmin=m ! >0 if there is CO2 growth all year
        if (subsurface) then
           Tbottom(Mx1:Mx2,My1:My2) = Tbottom(Mx1:Mx2,My1:My2)+T(nz,:,:)
        endif
        nm=nm+1

        write(22,'(f9.3,1x,f7.3,2x,f6.1,1x,f5.1,1x,f6.1)') &
             & sdays,mod(marsLs/d2r,360.d0),Qn(1,1),Tsurf(1,1),m(1,1) ! flat surf ref
     endif
     
     if (sdays > tmax-2*solsy) then  ! longest continuous period below H2O frost point
        where (Tsurf<Tfrost) 
           frosttime = frosttime+dt
        elsewhere
           frosttime = 0.
        end where
        where (frosttime>maxfrosttime)
           maxfrosttime = frosttime
           h2olast = marsLs  ! last time maxfrosttime increased
        end where
     endif

     do k=1,3
        if (jd+edays > jd_snap(k)-dt/2 .and. jd+edays <= jd_snap(k)+dt/2) then
           call writesnapshot(fns(k),h,Qdirect,m,Qn)
        endif
     end do
     do k=1,2
        if (jd+edays > jd_themis(k)-dt/2 .and. jd+edays <= jd_themis(k)+dt/2) then
           call writethemissnapshot(fnt(k),h,Tsurf)
        endif
     enddo
  
  enddo  ! end of time loop

  close(22)
  if (subsurface) deallocate(T)

  Qmean = Qmean/nm
  Tmean = Tmean/nm
  where (h2olast/=-9.) h2olast=h2olast/d2r

  open(unit=21,file='qmean.dat',status='unknown',action='write')
  do i=Mx1,Mx2
     do j=My1,My2
        write(21,'(2(i4,1x),f9.2,2x,f6.4,2(1x,f6.1),2(1x,f5.1),3(1x,f7.1),1x,f6.2)') &
             & i,j,h(i,j),surfaceSlope(i,j),Qmean(i,j),Qmax(i,j), &
             & Tmean(i,j),Tmaxi(i,j),mmax(i,j),maxfrosttime(i,j), &
             & mmin(i,j),h2olast(i,j)
     enddo
  enddo
  close(21)
  if (subsurface) then
     Tbottom=Tbottom/nm
     do i=Mx1,Mx2
        do j=My1,My2
           write(23,'(2(i5,1x),2(1x,f7.1))') i,j,Tmean(i,j),Tbottom(i,j)
        enddo
     enddo
  endif
  
end program cratersQ_mars



subroutine writethemissnapshot(fn,h,Tsurf)
  ! output surface temperature snapshot
  use filemanager, only : NSx, NSy
  implicit none
  character(len=*), intent(IN) :: fn
  real(8), intent(IN), dimension(NSx,NSy) :: h,Tsurf
  integer i,j

  print *,'entered writethemissnapshot'
  open(unit=27,file=fn,status='unknown',action='write')
  do i=2,NSx-1
     do j=2,NSy-1
        write(27,'(2(i4,1x),f9.2,1x,f5.1)') i,j,h(i,j),Tsurf(i,j)
     enddo
  enddo
  close(27)
end subroutine writethemissnapshot



subroutine writesnapshot(fn,h,Qdirect,m,Qn)
  ! output snapshot
  use filemanager, only : NSx, NSy
  implicit none
  character(len=*), intent(IN) :: fn
  real(8), intent(IN), dimension(NSx,NSy) :: h,Qdirect,m,Qn
  integer i,j

  print *,'entered writesnapshot'
  open(unit=27,file=fn,status='unknown',action='write')
  do i=2,NSx-1
     do j=2,NSy-1
        write(27,'(2(i4,1x),f9.2,1x,f6.1,1x,f7.1,1x,f6.1)') &
             & i,j,h(i,j),Qdirect(i,j),m(i,j),Qn(i,j)
     enddo
  enddo
  close(27)
end subroutine writesnapshot
