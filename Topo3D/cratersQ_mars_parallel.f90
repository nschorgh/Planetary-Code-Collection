!***************************************************************************
! cratersQ_mars_parallel:
!                Mars thermal model with direct insolation, subsurface
!                conduction, terrain shadowing, and approximate self-heating
!
! parallel version of cratersQ_mars
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
  integer, parameter :: nz=70
  real(8), parameter :: thIn = 300.  ! Thermal inertia
end module miscparams


program cratersQ_mars
  use filemanager
  use allinterfaces
  use miscparams
  use newhorizons
  implicit none

  real(8), parameter :: albedo0=0.12, co2albedo=0.65
  real(8) :: latitude = -41.6  ! Palikir Crater
  real(8), parameter :: longitude = 360 - 202.3 ! west longitude

  integer nsteps, n, i, j, nm
  integer Mx1, Mx2, My1, My2
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
  real(8), allocatable, dimension(:,:) :: mmin, co2last, co2first, h2olast
  real(8) Fsurf_flat, m_flat, Tsurf_flat, albedo_flat, Qn_flat, Qnm1_flat, Qdirect_flat
  real(8) dE, Tsurfold, Qlw, Qscat, QIR, Qrefl
  logical, parameter :: subsurface=.true.  ! control panel
  integer k, i0, j0
  real(8) jd, LTST, jd_end
  character(len=5) ext
  
  real(8) jd_snap(3), jd_themis(2)  ! Julian dates of snapshots
  character(len=20) fns(3), fnt(2)  ! file names of snapshots
  
  integer narg
  character(4) extc
  narg = iargc()
  if (narg==2) then ! parallel implementation
     call slicer(NSx,Mx1,Mx2,extc)
     My1=2; My2=NSy-1
     ext = '.'//extc
  else
     Mx1=2; Mx2=NSx-1; My1=2; My2=NSy-1 ! whole domain

     ! single point
     !Mx1=4369; My1=12850
     !Mx2=Mx1; My2=My1
     
     ! enlarged RSL region
     !My1=3500/4; My2=3660/4; Mx1=980/4; Mx2=1080/4 ! 16m
     !My1=3500/2; My2=3660/2; Mx1=980/2; Mx2=1080/2 ! 8m
     !My1=3500; My2=3660; Mx1=980; Mx2=1080 ! 4m
     !My1=3500*2; My2=3660*2; Mx1=980*2; Mx2=1080*2 ! 2m
     !My1=3500*4; My2=3660*4; Mx1=980*4; Mx2=1080*4 ! 1m
     
     i0=Mx1; j0=My1
     ext = '.dat'
  endif
  
  print *,'File extension ',ext
  
  allocate(h(NSx,NSy))
  allocate(surfaceSlope(Mx1:Mx2,My1:My2), azFac(Mx1:Mx2,My1:My2))
  allocate(Qn(Mx1:Mx2,My1:My2), source=zero); Qn_flat=0.
  allocate(Qnm1, Qdirect, source=Qn); Qnm1_flat=0.; Qdirect_flat=0.
  allocate(Tsurf(Mx1:Mx2,My1:My2))
  allocate(m(Mx1:Mx2,My1:My2), source=zero); m_flat=0.
  allocate(albedo(Mx1:Mx2,My1:My2), source=albedo0); albedo_flat=albedo0
  allocate(skyview(Mx1:Mx2,My1:My2)); skyview=1.
  allocate(gterm(Mx1:Mx2,My1:My2), source=zero)
  ! diagnostic variables
  allocate(Qmean, Qmax, source=Qn)
  allocate(Tmean(Mx1:Mx2,My1:My2), Tmaxi(Mx1:Mx2,My1:My2), source=zero)
  allocate(frosttime, maxfrosttime, mmax, mmin, source=m)
  allocate(co2last, co2first, h2olast, mold=m)
  co2last=-9.; co2first=-9.; h2olast=-9.; mmin=1e32
  
  dt=0.02
  tmax = 5.*solsy

  ! set some constants
  nsteps=int(tmax/dt)       ! calculate total number of timesteps
  dtsec = dt*solarDay
  
  write(*,*) 'Time step=',dt,' Max number of steps=',nsteps
  write(*,*) 'Calculations performed for latitude=',latitude
  write(*,*) 'fracIR=',fracIR,'fracDust=',fracDust
  write(*,*) 'Nx=',NSx,'Ny=',NSy,'File=',fileext
  write(*,*) 'Region of interest: (',Mx1,',',My1,') x (',Mx2,',',My2,')'
  write(*,*) 'Mean albedo=',sum(albedo)/size(albedo),'Emissivity=',emiss
  write(*,*) 'CO2 frost temperature=',Tco2frost,'CO2 albedo=',co2albedo
  write(*,*) 'Reflections:',.FALSE.,'Subsurface:',subsurface

  ! Set start date
  !jd=dble(julday(1,1,2009))  !  JD for noon UTC on imm,iday,iyear

  ! Alternatively set end date
  !jd_end=dble(julday(8,24,2014))  !  JD for noon UTC on imm,iday,iyear
  jd_end=dble(julday(2,28,2017))  ! includes last themis snapshot
  jd = nint(jd_end - tmax*solarDay/earthDay)

  write(*,*) 'Julian start and end date',jd,jd_end

  print *,'...reading topography...'
  call readdem(h)
  call difftopo2(h,surfaceSlope,azFac,Mx1,Mx2,My1,My2)

  latitude=latitude*d2r
  Tsurf=200.
  nm=0

  print *,'...reading horizons file...'
  call readhorizons(Mx1,Mx2,My1,My2)
  do concurrent(i=Mx1:Mx2, j=My1:My2)
     skyview(i,j) = getoneskysize_v2(i,j)/(2*pi)
     gterm(i,j) = getoneGterm(i,j,surfaceSlope(i,j),azFac(i,j))
  end do
  
  if (subsurface) then ! initialize subsurface component
     allocate(T(nz,Mx1:Mx2,My1:My2))
     allocate(Tref(nz))
     allocate(Fsurf(Mx1:Mx2,My1:My2))
     call subsurfaceconduction_mars(Tref(:),buf,dtsec,zero,zero,buf,buf,.true.,thIn=thIn)
     allocate(Tbottom(Mx1:Mx2,My1:My2))
     Tbottom(:,:)=-9.
     Tsurf(:,:)=-9.  ! max 3 digits
     Tsurf_flat=-9.
     Fsurf(:,:)=0.; Fsurf_flat = 0.
  end if

  open(unit=22,file='timeseries_flat'//ext,status='unknown',action='write')
  open(unit=25,file='timeseries_pnt'//ext,status='unknown',action='write')

  ! image taken imm = 5; iday=15; iyr=2014 8:44 UTC  ESP_036561
  jd_snap(1)=dble(julday(5,15,2014)) + (8.+44./60-12)/24.  !  JD for noon UTC on imm,iday,iyear
  fns(1)='qsnap_036561'//ext ! snapshot
  
  jd_snap(2)=dble(julday(8,18,2012)) + (8.+48./60-12)/24.; fns(2)='qsnap_028412'//ext
  jd_snap(3)=dble(julday(5,11,2012)) + (17.+31./60-12)/24.; fns(3)='qsnap_027146'//ext

  jd_themis(1)=dble(julday(10,30,2016)) + (1.+30./60-12)/24.; fnt(1)='qsnap_i65997002'//ext
  jd_themis(2)=dble(julday(2,26,2017)) + (18.+03./60-12)/24.; fnt(2)='qsnap_i67450002'//ext 


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
          & zero,zero,zero,Qdirect_flat,Qscat,Qlw)
     Qn_flat = (1-albedo_flat)*(Qdirect_flat+Qscat) + emiss*Qlw
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
              QIR = gterm(i,j)*emiss*sigSB*Tsurf_flat**4
              Qn(i,j) = Qn(i,j) + emiss*QIR
           endif
           Qrefl = gterm(i,j)*albedo_flat*(Qdirect_flat+Qscat)
           Qn(i,j) = Qn(i,j) + (1-albedo(i,j))*Qrefl
        enddo
     enddo
     if (n==0) then
        Qnm1(:,:) = Qn(:,:)
        Qnm1_flat = Qn_flat
     endif
     
     if (subsurface) then ! subsurface conduction
        do i=Mx1,Mx2
           do j=My1,My2
              if (h(i,j)<-32000) cycle
              call subsurfaceconduction_mars(T(:,i,j),Tsurf(i,j), &
                   & dtsec,Qnm1(i,j),Qn(i,j),m(i,j),Fsurf(i,j),.false., &
                   & Tco2frost=Tco2frost,emiss=emiss)
           enddo
        enddo
        call subsurfaceconduction_mars(Tref(:),Tsurf_flat, &
             & dtsec,Qnm1_flat,Qn_flat,m_flat,Fsurf_flat,.false.,Tco2frost=Tco2frost,emiss=emiss)

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
        Tsurf_flat = (Qn_flat/emiss/sigSB)**0.25
        Tsurfold = (Qnm1_flat/emiss/sigSB)**0.25
        if (Tsurf_flat<Tco2frost.or.m_flat>0.) then   ! CO2 condensation
           Tsurf_flat=Tco2frost
           dE = - Qn_flat + emiss*sigSB*(Tsurf_flat**4 + Tsurfold**4)/2.
           m_flat = m_flat + dtsec*dE/Lco2frost
        endif
     endif
     Qnm1(:,:) = Qn(:,:)
     Qnm1_flat = Qn_flat

     where (Tsurf>Tco2frost .or. m<=0.) ! all arrays must be of same size
        albedo = albedo0
     elsewhere
        albedo = co2albedo
     end where
     if (Tsurf_flat>Tco2frost .or. m_flat<=0.) then
        albedo_flat = albedo0
     else
        albedo_flat = co2albedo
     end if
     
     ! only diagnostics below this line
     if (sdays > tmax-solsy) then
        ! all arrays must be of same size
        Qmean = Qmean + Qn
        where (Qn>Qmax) Qmax=Qn
        Tmean = Tmean + Tsurf
        where (Tsurf>Tmaxi) Tmaxi=Tsurf
        where (m>mmax) mmax=m
        where (m<mmin) mmin=m ! >0 if there is CO2 growth all year
        where (m>0. .and. co2first<0.) co2first = marsLs
        where (m>0.) co2last = marsLs
        if (subsurface) Tbottom(:,:) = Tbottom(:,:)+T(nz,:,:)
        nm=nm+1

        write(22,'(f9.3,1x,f7.3,2x,f6.1,1x,f5.1,1x,f6.1)') &
             & sdays,mod(marsLs/d2r,360.d0),Qn_flat,Tsurf_flat,m_flat ! flat surf ref
        !if (i0>=Mx1 .and. i0<=Mx2 .and. j0>=My1 .and. j0<=My2) then
        !   write(25,'(f9.3,1x,f7.3,2(1x,i5),2x,f6.1,1x,f5.1,1x,f6.1)') &
        !        & sdays,mod(marsLs/d2r,360.d0),i0,j0,Qn(i0,j0),Tsurf(i0,j0),m(i0,j0)
        !end if
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
           call writesnapshot(fns(k),h,Qdirect,m,Qn,Mx1,Mx2,My1,My2)
        endif
     end do
     do k=1,2
        if (jd+edays > jd_themis(k)-dt/2 .and. jd+edays <= jd_themis(k)+dt/2) then
           call writethemissnapshot(fnt(k),h,Tsurf,Mx1,Mx2,My1,My2)
        endif
     enddo

  enddo  ! end of time loop

  close(22)
  close(25)
  if (subsurface) deallocate(T)

  Qmean = Qmean/nm
  Tmean = Tmean/nm
  where (co2first/=-9.) co2first=co2first/d2r
  where (co2last/=-9.) co2last=co2last/d2r
  where (h2olast/=-9.) h2olast=h2olast/d2r

  open(unit=21,file='qmean'//ext,status='unknown',action='write')
  do i=Mx1,Mx2
     do j=My1,My2
        write(21,'(2(i5,1x),f9.2,2x,f6.3,2(1x,f6.1),2(1x,f5.1),3(1x,f7.1),3(1x,f6.2),1x,f5.1,1x,f6.2)') &
             & i,j,h(i,j),surfaceSlope(i,j),Qmean(i,j),Qmax(i,j), &
             & Tmean(i,j),Tmaxi(i,j),mmax(i,j),maxfrosttime(i,j), &
             & mmin(i,j),co2first(i,j),co2last(i,j),h2olast(i,j)
     enddo
  enddo

  close(21)
  if (subsurface) then
     open(unit=23,file='tsurfbot'//ext,status='unknown',action='write')
     Tbottom=Tbottom/nm
     do i=Mx1,Mx2
        do j=My1,My2
           write(23,'(2(i5,1x),2(1x,f7.1))') i,j,Tmean(i,j),Tbottom(i,j)
        enddo
     enddo
     close(23)
  endif

end program cratersQ_mars



subroutine writethemissnapshot(fn,h,Tsurf,Mx1,Mx2,My1,My2)
  ! output surface temperature snapshot
  use filemanager, only : NSx, NSy
  implicit none
  integer, intent(IN) :: Mx1,Mx2,My1,My2
  character(len=*), intent(IN) :: fn
  real(8), intent(IN), dimension(NSx,NSy) :: h
  real(8), intent(IN), dimension(Mx1:Mx2,My1:My2) :: Tsurf
  integer i,j

  print *,'entered writethemissnapshot'
  open(unit=27,file=fn,status='unknown',action='write')
  do i=max(2,Mx1),min(NSx-1,Mx2)
     do j=max(2,My1),min(NSy-1,My2)
        write(27,'(2(i5,1x),f9.2,1x,f5.1)') i,j,h(i,j),Tsurf(i,j)
     enddo
  enddo
  close(27)
end subroutine writethemissnapshot



subroutine writesnapshot(fn,h,Qdirect,m,Qn,Mx1,Mx2,My1,My2)
  ! output snapshot
  use filemanager, only : NSx, NSy 
  implicit none
  integer, intent(IN) :: Mx1,Mx2,My1,My2
  character(len=*), intent(IN) :: fn
  real(8), intent(IN), dimension(NSx,NSy) :: h
  real(8), intent(IN), dimension(Mx1:Mx2,My1:My2) :: Qdirect,m,Qn
  integer i,j

  print *,'entered writesnapshot'
  open(unit=27,file=fn,status='unknown',action='write')
  do i=max(2,Mx1),min(NSx-1,Mx2)
     do j=max(2,My1),min(NSy-1,My2)
        write(27,'(2(i5,1x),f9.2,1x,f6.1,1x,f7.1,1x,f6.1)') &
             & i,j,h(i,j),Qdirect(i,j),m(i,j),Qn(i,j)
     enddo
  enddo
  close(27)
end subroutine writesnapshot
