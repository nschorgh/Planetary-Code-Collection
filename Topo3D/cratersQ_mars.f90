! Mars thermal model with horizons from 3D topography

module miscparams
  real(8), parameter :: pi=3.1415926535897932, d2r=pi/180.
  real(8), parameter :: sigSB = 5.6704e-8
  real(8), parameter :: Tco2frost=145., Lco2frost=6.0e5  ! Mars
  real(8), parameter :: zero = 0.
  real(8), parameter :: Tfrost = 200. ! H2O frost point temperature, for diagnostics only
  real(8), parameter :: earthDay = 86400.

  ! thermal model parameters
  real(8), parameter :: fracIR=0.04, fracDust=0.02
  real(8), parameter :: solsy = 668.60 ! solar days per Mars year
  real(8), parameter :: solarDay = 88775.244  ! Mars
  real(8), parameter :: Fgeotherm = 0.0
  integer, parameter :: nz=70
end module miscparams


program cratersQ_mars
  use filemanager
  use allinterfaces
  use miscparams
  use newhorizons
  implicit none

  real(8), parameter :: albedo0=0.12, co2albedo=0.65
  real(8), parameter :: emiss = 0.98
  real(8) :: latitude = -41.6
  real(8), parameter :: longitude = 360 - 202.3 ! west longitude

  integer nsteps, n, i, j, nm
  real(8) tmax, dt, dtsec, buf
  real(8) HA, sdays, azSun, emax, sinbeta
  real(8) edays, marsR, marsLs, marsDec
  real(8), allocatable, dimension(:,:) :: h, surfaceSlope, azFac
  real(8), allocatable, dimension(:,:) :: Qn   ! incoming
  real(8), allocatable, dimension(:,:) :: Tsurf, albedo, m
  real(8), allocatable, dimension(:,:) :: Qmean, Qmax, Tmean, Tmaxi, Qdirect
  real(8), allocatable, dimension(:,:) :: viewfactor, gterm
  real(8), allocatable, dimension(:,:) :: mmax, frosttime, maxfrosttime, Qnm1
  real(8), allocatable :: Fsurf(:,:), T(:,:,:), Tbottom(:,:), Tref(:)  ! subsurface
  real(8), allocatable, dimension(:,:) :: mmin, co2last, co2first, h2olast, h2olastt !, EH2Ocum
  real(8) dE, Tsurfold, QIR, Qscat, Qlw !, EH2O
  logical, parameter :: subsurface=.true.  ! control panel
  !integer, parameter :: NrMP=3   ! number of monitoring points
  integer k !, i0, j0, i00(NrMP), j00(NrMP)
  integer, external :: julday
  real(8) jd, jd_snap(3), LTST, jd_end, jd_themis(2)
  character(len=20) fnt(2), fns(3)  ! snapshot file names
  integer, parameter :: Mx1=2, Mx2=NSx-1, My1=2, My2=NSy-1
  
  allocate(h(NSx,NSy), surfaceSlope(NSx,NSy), azFac(NSx,NSy))
  allocate(Qn(NSx,NSy), Qnm1(NSx,NSy), Qdirect(NSx,NSy))
  allocate(Tsurf(NSx,NSy), m(NSx,NSy))
  allocate(albedo(NSx,NSy)); albedo = albedo0
  allocate(Fsurf(NSx,NSy))
  allocate(viewfactor(Mx1:Mx2,My1:My2), gterm(Mx1:Mx2,My1:My2))
  allocate(Qmean(NSx,NSy), Qmax(NSx,NSy), Tmean(NSx,NSy), Tmaxi(NSx,NSy))
  allocate(frosttime(NSx,NSy), maxfrosttime(NSx,NSy))
  allocate(mmax(NSx,NSy), mmin(NSx,NSy))
  allocate(co2last(NSx,NSy), co2first(NSx,NSy), h2olast(NSx,NSy), h2olastt(NSx,NSy))
  !allocate(EH2Ocum(Mx1:Mx2,My1:My2))
  Qn=0.; Qnm1=0.; Fsurf=0.; m=0.; Qdirect=0.; viewfactor=0.
  Qmean=0.; Qmax=0.; Tmean=0.; Tmaxi=0.
  frosttime=0.; maxfrosttime=0.; mmax=0.; mmin=0
  co2last=-9.; co2first=-9.; h2olast=-9.; h2olastt=-9.
  
  dt=0.02
  !tmax = 2*solsy+1.
  !tmax = solsy*10.5  ! should end at the beginning of spring for the respective hemisphere
  !tmax = 2.
  tmax = 5.*solsy

  ! set some constants
  nsteps=int(tmax/dt)       ! calculate total number of timesteps
  dtsec = dt*solarDay
  
  write(*,*) 'Time step=',dt,' Max number of steps=',nsteps
  write(*,*) 'Calculations performed for latitude=',latitude
  write(*,*) 'Nx=',NSx,'Ny=',NSy,'File=',fileext
  write(*,*) 'Region of interest: (',Mx1,',',My1,') x (',Mx2,',',My2,')'
  write(*,*) 'Mean albedo=',sum(albedo)/size(albedo),'Emissivity=',emiss
  write(*,*) 'Reflections:',.FALSE.,'Subsurface:',subsurface

  ! Set start date
  !jd=dble(julday(1,1,2009))  !  JD for noon UTC on imm,iday,iyear
  !call marsorbit(dt0_j2000,0.d0,marsLs,marsDec,marsR)
  !call marsclock24(jd,dt0_J2000,marsLs,marsDec,marsR,0d0,LTST) ! calculate dt0_J2000

  ! Alternatively set end date
  jd_end=dble(julday(2,28,2017)) 
  jd = nint(jd_end - tmax*solarDay/earthDay)
  
  call readdem(h)
  call difftopo(NSx,NSy,h,dx,dy,surfaceSlope,azFac)

  latitude=latitude*d2r
  Tsurf=200.
  nm=0
  !EH2Ocum = 0.

  !open(unit=31,file='monitorpoints.dat',action='read')
  !do k=1,NrMP
  !   read(31,*) i00(k),j00(k)
  !enddo
  !close(31)
  
  print *,'...reading horizons file...'
  call readhorizons('horizons.'//fileext)
  ! only use viewfactors with IR contribution
  !call getskysize(viewfactor)
  do concurrent(i=2:NSx-1, j=2:Nsy-1)
     viewfactor(i,j) = getoneskysize(i,j)
     gterm(i,j) = getoneGterm(i,j,surfaceSlope(i,j),azFac(i,j))
  end do
  viewfactor = viewfactor/(2*pi)
  
  if (subsurface) then
     allocate(T(Mx1:Mx2,My1:My2,nz))
     allocate(Tref(nz))
     call subsurfaceconduction_mars(Tref(:),buf,dtsec,zero,zero,zero,buf,buf,.true.)
     allocate(Tbottom(NSx,NSy))
     Tbottom(:,:)=-9
     Tsurf(:,:)=-9  ! max 3 digits
  end if

  open(unit=22,file='timeseries_flat.dat',status='unknown',action='write')

  ! image taken imm = 5; iday=15; iyr=2014 8:44 UTC  ESP_036561
  jd_snap(1)=dble(julday(5,15,2014)) + (8.+44./60-12)/24.  !  JD for noon UTC on imm,iday,iyear
  !call marsclock24(jd_snap,buf,marsLs,marsDec,marsR,Longitude,LTST)
  fns(1)='qsnap_036561.dat' ! snapshot
  
  jd_snap(2)=dble(julday(8,18,2012)) + (8.+48./60-12)/24.  
  fns(2)='qsnap_028412.dat'

  jd_snap(3)=dble(julday(5,11,2012)) + (17.+31./60-12)/24. 
  fns(3)='qsnap_027146.dat'

  jd_themis(1)=dble(julday(2,13,2017)) + (22.+05./60-12)/24. 
  fnt(1) = 'qsnap_i67294006.dat'

  jd_themis(2)=dble(julday(2,24,2017)) + (16.+36./60-12)/24. 
  fnt(2) = 'qsnap_i67425002.dat' 

  print *,'...calculating...'
  ! loop over time steps 
  do n=0,nsteps-1
     sdays = (n+1)*dtsec/solarDay
     edays = (n+1)*dtsec/earthDay

     !call marsorbit(dt0_j2000,edays,marsLs,marsDec,marsR)
     !call generalorbit(edays,a,ecc,omega,eps,marsLs,marsDec,marsR)
     !HA=2.*pi*mod(sdays,1.)   ! hour angle
     call marsclock24(jd+edays,buf,marsLs,marsDec,marsR,Longitude,LTST)
     HA=2.*pi*mod(LTST+12,24.d0)/24
     
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
           Qn(i,j) = (1-albedo(i,j))*(Qdirect(i,j)+Qscat*viewfactor(i,j)) &
                & + emiss*Qlw*viewfactor(i,j)
           ! contribution from land in field of view
           if (n>0) then 
              !QIR = (1-viewfactor(i,j))*emiss*sigSB*Tsurf(1,1)**4
              QIR = gterm(i,j)*emiss*sigSB*Tsurf(1,1)**4
              Qn(i,j) = Qn(i,j) + emiss*QIR
           endif
           !Qrefl = (1-viewfactor(i,j))*albedo(1,1)*(Qdirect(1,1)+Qscat)
           !Qrefl = gterm(i,j)*albedo(1,1)*(Qdirect(1,1)+Qscat)
           !Qn(i,j) = Qn(i,j) + (1-albedo(i,j))*Qrefl(i,j)
           ! Qdirect w/o atmosphere is for snapshot; different from Qdirect above
           !Qdirect(i,j)=flux_wshad(marsR,sinbeta,azSun,surfaceSlope(i,j),azFac(i,j),emax)
        enddo
     enddo
     if (n==0) Qnm1(:,:) = Qn(:,:)
     
     if (subsurface) then ! subsurface conduction
        do i=Mx1,Mx2
           do j=My1,My2
              if (h(i,j)<-32000) cycle
              call subsurfaceconduction_mars(T(i,j,:),Tsurf(i,j), &
                   & dtsec,Qnm1(i,j),Qn(i,j),emiss,m(i,j),Fsurf(i,j),.false.)
           enddo
        enddo
        call subsurfaceconduction_mars(Tref(:),Tsurf(1,1), &
             & dtsec,Qnm1(1,1),Qn(1,1),emiss,m(1,1),Fsurf(1,1),.false.)

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
        ! Tsurf(1,1) not implemented
     endif
     Qnm1(:,:) = Qn(:,:)

     where (Tsurf>Tco2frost.or.m<=0.)
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
        where (m>0. .and. co2first<0.) co2first = marsLs
        where (m>0.) co2last = marsLs        
        if (subsurface) then
           Tbottom(Mx1:Mx2,My1:My2) = Tbottom(Mx1:Mx2,My1:My2)+T(:,:,nz)
        endif
        nm=nm+1

        write(22,'(f9.3,1x,f7.3,2x,f6.1,1x,f5.1,1x,f6.1)') &
             & sdays,mod(marsLs/d2r,360.d0),Qn(1,1),Tsurf(1,1),m(1,1) ! flat surf ref
        !do k=1,size(i00)
           !i0=i00(k); j0=j00(k)
           !if (i0<Mx1 .or. i0>Mx2 .or. j0<My1 .or. j0>My2) cycle
           !EH2O = evap_ingersoll(Tsurf(i0,j0))
           !if (EH2O==EH2O) EH2Ocum = EH2Ocum + EH2O*dtsec   ! skips NaNs
           !EH2Ocum(i0,j0) = EH2Ocum(i0,j0) + EH2O*dtsec 
           !write(24,'(f9.3,1x,f7.3,2(1x,i4),2x,f6.1,1x,f5.1,2(1x,g12.4))') &
           !     & sdays,mod(marsLs/d2r,360.),i0,j0,Qn(i0,j0),Tsurf(i0,j0),EH2O,EH2Ocum(i0,j0)
        !enddo

        do i=2,NSx-1 ! second half
           do j=2,NSy-1
              if (Tsurf(i,j)>263. .and. sdays-h2olastt(i,j)<14.) then
                 write(40,*) i,j,h2olastt(i,j),sdays
              endif
           end do
        end do
     endif
     
     if (sdays > tmax-2*solsy) then  ! longest continuous period below H2O frost point (~200K)
        where (Tsurf<Tfrost) 
           frosttime=frosttime+dt
        elsewhere
           frosttime=0.
        end where
        where (frosttime>maxfrosttime)
           maxfrosttime=frosttime
           h2olast = marsLs  ! last time maxfrosttime increased
           h2olastt = sdays ! first half
           !EH2Ocum = 0.
        end where
     endif

     do k=1,3
        if (jd+edays > jd_snap(k)-dt/2 .and. jd+edays <= jd_snap(k)+dt/2) then
           call writesnapshot(fns(k),h,Qdirect,m,Qn,NSx,NSy)
        endif
     end do
     do k=1,2
        if (jd+edays > jd_themis(k)-dt/2 .and. jd+edays <= jd_themis(k)+dt/2) then
           call writethemissnapshot(fnt(k),h,Tsurf,NSx,NSy)
        endif
     enddo
  
  enddo  ! end of time loop

  close(22)
  if (subsurface) deallocate(T)

  Qmean = Qmean/nm
  Tmean = Tmean/nm
  where (co2first/=-9.) co2first=co2first/d2r
  where (co2last/=-9.) co2last=co2last/d2r
  where (h2olast/=-9.) h2olast=h2olast/d2r

  open(unit=21,file='qmean.dat',status='unknown',action='write')
  do i=Mx1,Mx2
     do j=My1,My2
        write(21,'(2(i4,1x),f9.2,2x,f6.4,2(1x,f6.1),2(1x,f5.1),3(1x,f7.1),3(1x,f6.2))') &
             & i,j,h(i,j),surfaceSlope(i,j),Qmean(i,j),Qmax(i,j), &
             & Tmean(i,j),Tmaxi(i,j),mmax(i,j),maxfrosttime(i,j), &
             & mmin(i,j),co2first(i,j),co2last(i,j),h2olast(i,j)
     enddo
  enddo
  close(21)
  if (subsurface) then
     Tbottom=Tbottom/nm
     do i=Mx1,Mx2
        do j=My1,My2
           write(23,'(2(i4,1x),2(1x,f7.1))') i,j,Tmean(i,j),Tbottom(i,j)
        enddo
     enddo
  endif
  
end program cratersQ_mars



subroutine subsurfaceconduction_mars(T,Tsurf,dtsec,Qn,Qnp1,emiss,m,Fsurf,init)
  use miscparams
  use conductionQ
  use conductionT
  implicit none
  real(8), intent(INOUT) :: T(nz), Tsurf, m, Fsurf
  real(8), intent(IN) :: dtsec,Qn,Qnp1,emiss
  logical, intent(IN) :: init
  integer i
  !real(8), parameter :: zmax=3., zfac=1.05d0  ! adjust
  real(8), parameter :: zmax=13., zfac=1.05d0  ! with rhoc=thIn*1000 (nz=70, 3x seasonal)
  real(8) Tinit, delta, thIn
  real(8) Fsurfold, dE, Tsurfold, Told(1:nz)
  real(8) z(nz), ti(nz), rhocv(nz)

  if (init) then ! initialize grid
     thIn = 400.
     ti(:) = thIn  ! adjust
     !rhocv(:) = 1200.*800.
     rhocv(:) = thIn*1000.  ! makes skin depth invariant
     
     delta = thIn/rhocv(1)*sqrt(solarDay/pi)  ! skin depth

     call setgrid(nz,z,zmax,zfac)
     if (z(6)>delta) then
        print *,'WARNING: fewer than 6 points within diurnal skin depth'
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
     write(*,*) '   Thermal inertia=',thIn,' rho*c=',rhocv(1)
     write(*,*) '   Diurnal and seasonal skin depths=',delta,delta*sqrt(solsy)
     write(*,*) '   Geothermal flux=',Fgeotherm

     call conductionT2_init(nz,z,dtsec,ti,rhocv,Fgeotherm)
     call conductionQ2_init(nz,z,dtsec,ti,rhocv,Fgeotherm)
     
     return
  endif
  
  if (Tsurf<=0.) then  ! initialize temperature profile
     Tinit=200.
     T(1:nz) = Tinit
     Tsurf = Tinit
  endif

  Tsurfold=Tsurf
  Fsurfold=Fsurf
  Told(1:nz)=T(1:nz)
  if (Tsurf>Tco2frost.or.m<=0.) then
     call conductionQ2(nz,Qn,Qnp1,T,emiss,Tsurf,Fsurf)
  endif
  if (Tsurf<Tco2frost.or.m>0.) then   ! CO2 condensation                                              
     T(1:nz)=Told
     call conductionT2(nz,Tsurfold,Tco2frost,T,Fsurf)
     Tsurf=Tco2frost
     dE = (- Qn - Qnp1 + Fsurfold + Fsurf + &
          &           emiss*sigSB*(Tsurfold**4+Tsurf**4))/2.
     m = m + dtsec*dE/Lco2frost
  endif

end subroutine subsurfaceconduction_mars



elemental function evap_ingersoll(T)
  ! Returns sublimation rate (kg/m^2/s)
  use allinterfaces, only : psv
  implicit none
  real(8) evap_ingersoll
  real(8), intent(IN) :: T
  real(8) p0,psat,D,Gbuf,rho,R,rhow,nu,drhooverrho,g

  p0=520.  ! atmospheric pressure
  psat=psv(T)
  R=8314.5
  D=20e-4 ! in m^2/s 
  g=3.7
  rhow = psat*18/(R*T)
  rho = p0*44/(R*T)
  !nu=7e-4  ! kinematic viscosity of CO2
  nu=16e-4*(T/200)**1.5
  
  !drhooverrho=(44-18)*psat/(44*p0-(44-18)*psat) ! Ingersoll (1970)
  drhooverrho=(44-18)*psat/(44*(p0-psat)) ! diverges at p0=psat
  Gbuf=(drhooverrho*g/nu**2)**(1./3.);
  !evap_ingersoll=0.17*D*rhow*Gbuf  ! Ingersoll (1970)
  evap_ingersoll=0.13*D*rhow*Gbuf
  
end function evap_ingersoll



subroutine writethemissnapshot(fn,h,Tsurf,NSx,NSy)
  ! output surface temperature snapshot
  implicit none
  integer, intent(IN) :: NSx,NSy
  character(len=*), intent(IN) :: fn
  real(8), intent(IN), dimension(NSx,NSy) :: h,Tsurf
  integer i,j
  
  open(unit=27,file=fn,status='unknown',action='write')
  do i=2,NSx-1
     do j=2,NSy-1
        write(27,'(2(i4,1x),f9.2,1x,f5.1)') i,j,h(i,j),Tsurf(i,j)
     enddo
  enddo
  close(27)
end subroutine writethemissnapshot



subroutine writesnapshot(fn,h,Qdirect,m,Qn,NSx,NSy)
  ! output snapshot
  implicit none
  integer, intent(IN) :: NSx,NSy
  character(len=*), intent(IN) :: fn
  real(8), intent(IN), dimension(NSx,NSy) :: h,Qdirect,m,Qn
  integer i,j

  open(unit=27,file=fn,status='unknown',action='write')
  do i=2,NSx-1
     do j=2,NSy-1
        write(27,'(2(i4,1x),f9.2,1x,f6.1,1x,f7.1,1x,f6.1)') &
             & i,j,h(i,j),Qdirect(i,j),m(i,j),Qn(i,j)
     enddo
  enddo
  close(27)
end subroutine writesnapshot
