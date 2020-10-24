!***************************************************************************
! cratersQ_mars_full: Mars thermal model with direct insolation, subsurface
!                     conduction, terrain shadowing, and full 3D reflection
!
! based on cratersQ_mars.f90 but with full 3D reflection
!***************************************************************************


module miscparams
  real(8), parameter :: pi=3.1415926535897932, d2r=pi/180.
  real(8), parameter :: sigSB = 5.6704e-8
  real(8), parameter :: Lco2frost = 6.0e5  ! [J/kg]
  real(8), parameter :: zero = 0.
  real(8), parameter :: solsy = 669. ! solar days per Mars year
  real(8), parameter :: solarDay = 88775.244  ! Mars [s]

  ! thermal model parameters
  real(8), parameter :: Tco2frost = 145. ! adjust according to elevation [K]
  real(8), parameter :: Tfrost = 200.  ! H2O frost point temperature, for diagnostics only
  real(8), parameter :: fracIR = 0.04, fracDust = 0.02
  real(8), parameter :: emiss = 0.98d0
  integer, parameter :: nz = 70
  real(8), parameter :: thIn = 400.
end module miscparams


PROGRAM cratersQ_mars
  use filemanager
  use allinterfaces
  use miscparams
  use newhorizons
  implicit none

  real(8), parameter :: albedo0=0.12d0, co2albedo=0.65d0
  real(8), parameter :: earthDay = 86400.

  integer nsteps, n, i, j, nm
  real(8) tmax, dt, latitude, dtsec, buf
  real(8) HA, sdays, azSun, emax, sinbeta
  real(8) edays, marsR, marsLs, marsDec
  real(8), dimension(NSx,NSy) :: h, surfaceSlope, azFac
  real(8), dimension(NSx,NSy) :: Qabs, Tsurf, albedo, m
  real(8), dimension(NSx,NSy) :: Qmean, Qmax, Tmean, Tmaxi
  real(8), dimension(NSx,NSy) :: skyview
  real(8), dimension(NSx,NSy) :: mmax, frosttime, maxfrosttime, Qnm1, mmin
  real(8), dimension(NSx,NSy) :: h2olast, meltfirst, h2olast_sols, meltfirst_sols, del
  real(8), allocatable :: Fsurf(:,:), T(:,:,:), Tbottom(:,:)  ! subsurface

  integer CCMAX, iii, jjj, k
  real(8), dimension(NSx,NSy) :: QIR, Qrefl, Qdirect, Qlw, Qscat
  real(8), dimension(NSx,NSy) :: Qvis, viewsize, QIRin, QIRre
  integer, dimension(NSx,NSy) :: cc
  integer(2), dimension(:,:,:), allocatable :: ii,jj
  real(4), dimension(:,:,:), allocatable :: VF
  
  real(8) jd, longitude, LTST, jd_end
  
  logical, parameter :: subsurface=.true.  ! control panel
  
  albedo = albedo0
  
  ! initializations
  Qabs=0.; Qnm1=0.; m=0.; skyview=1.
  Qmean=0.; Qmax=0.; Tmean=0.; Tmaxi=0.
  frosttime=0.; maxfrosttime=0.;
  mmax=0.; mmin=1e32
  h2olast=-9.; meltfirst=-9.
  h2olast_sols=-9; meltfirst_sols=-9.; del=-9.
  
  dt=0.02
  tmax = 6.*solsy
  latitude = -30.

  ! set some constants
  nsteps=int(tmax/dt)       ! calculate total number of timesteps
  dtsec = dt*solarDay
  
  write(*,*) 'Time step=',dt,' Max number of steps=',nsteps
  write(*,*) 'Calculations performed for latitude=',latitude
  write(*,*) 'fracIR=',fracIR,'fracDust=',fracDust
  write(*,*) 'Nx=',NSx,'Ny=',NSy,'File=',fileext
  write(*,*) 'Mean albedo=',sum(albedo)/size(albedo),'Emissivity=',emiss
  write(*,*) 'CO2 frost temperature=',Tco2frost,'CO2 albedo=',co2albedo
  write(*,*) 'Reflections:',.true.,'Subsurface:',subsurface

  ! Set start date
  !jd=dble(julday(1,1,2009))  !  JD for noon UTC on imm,iday,iyear
  !call marsorbit(dt0_j2000,0.d0,marsLs,marsDec,marsR)
  !call marsclock24(jd,dt0_J2000,marsLs,marsDec,marsR,0d0,LTST) ! calculate dt0_J2000

  ! Alternatively set end date
  !jd_end=dble(julday(8,24,2014))  !  JD for noon UTC on imm,iday,iyear
  jd_end = dble(julday(3,23,2019))  ! northern spring equinox (Ls=0)
  jd = jd_end - tmax*solarDay/earthDay
  
  call readdem(h)
  call difftopo(NSx,NSy,h,dx,dy,surfaceSlope,azFac)

  latitude=latitude*d2r
  Tsurf=200.
  nm=0

  print *,'...reading horizons file...'
  call readhorizons(2,NSx-1,2,NSy-1)
  ! only use skyviews with IR contribution
  do concurrent(i=2:NSx-1, j=2:NSy-1)
     !skyview(i,j) = getoneskysize_v2(i,j)/(2*pi)
     !gterm(i,j) = getoneGterm(i,j,surfaceSlope(i,j),azFac(i,j))
     skyview(i,j) = 1.-getoneGterm(i,j,surfaceSlope(i,j),azFac(i,j))
  end do
  print *,'max/min of skyview:',maxval(skyview(2:NSx-1,2:NSy-1)), &
       & minval(skyview(2:NSx-1,2:NSy-1))

  print *,'...reading huge fieldofviews file...'
  CCMAX = getmaxfieldsize(NSx,NSy,vfn)
  print *,'... max field of view size=',CCMAX
  allocate(ii(NSx,NSy,CCMAX), jj(NSx,NSy,CCMAX), VF(NSx,NSy,CCMAX))
  call getviewfactors(NSx,NSy,vfn,cc,ii,jj,VF,viewsize,CCMAX)
  
  Qrefl=0.; QIRre=0. 
  
  if (subsurface) then
     allocate(Fsurf(NSx,NSy), T(nz,NSx,NSy), Tbottom(NSx,NSy))
     Tbottom(:,:) = 0.
     Tsurf(:,:) = 200.
     Fsurf(:,:) = 0.
     ! initialize private variables in module
     call subsurfaceconduction_mars(T(:,1,1),buf,dtsec,zero,zero,buf,buf,.true.,thIn=thIn)
     T(:,:,:) = 200.
  end if

  longitude = 360 - 202.3 ! west longitude
  open(unit=25,file='timeseries_pnt.dat',status='unknown',action='write')

  print *,'...calculating...'
  ! loop over time steps 
  do n=0,nsteps-1
     sdays = (n+1)*dtsec/solarDay
     edays = (n+1)*dtsec/earthDay

     !call marsorbit(dt0_j2000,edays,marsLs,marsDec,marsR)
     !call generalorbit(edays,a,ecc,omega,eps,marsLs,marsDec,marsR)
     !HA=2.*pi*mod(sdays,1.d0)   ! hour angle
     call marsclock24(jd+edays,buf,marsLs,marsDec,marsR,Longitude,LTST)
     HA=2.*pi*mod(LTST+12,24.d0)/24
     
     call equatorial2horizontal(marsDec,latitude,HA,sinbeta,azSun)
     
     if (mod(n,10)==0) print *,n,sdays,marsLs/d2r

     do i=2,NSx-1
        do j=2,NSy-1
           ! incoming solar flux
           if (h(i,j)<-32000) cycle
           emax = getonehorizon(i,j,azSun)
           call flux_mars2(marsR,marsDec,latitude,HA,fracIR,fracDust, &
                & surfaceSlope(i,j),azFac(i,j),emax,Qdirect(i,j),Qscat(i,j),Qlw(i,j))
        enddo
     enddo

     Qvis = Qdirect + Qrefl + Qscat*skyview
     QIRin = QIR + QIRre + Qlw*skyview
     do i=2,NSx-1
        do j=2,NSy-1
           QIR(i,j)=0.; Qrefl(i,j)=0.; QIRre(i,j)=0.
           do k=1,cc(i,j)
              iii = ii(i,j,k); jjj = jj(i,j,k)
              Qrefl(i,j) = Qrefl(i,j) + VF(i,j,k)*albedo(iii,jjj)*Qvis(iii,jjj)
              QIR(i,j) = QIR(i,j) + VF(i,j,k)*emiss*sigSB*Tsurf(iii,jjj)**4
              QIRre(i,j) = QIRre(i,j) + VF(i,j,k)*(1-emiss)*QIRin(iii,jjj)
           enddo
        enddo
     enddo
     if (n==0) then
        Qnm1 = (1.-albedo)*(Qdirect+Qscat*skyview)+emiss*(QIR+Qlw*skyview)
     else
        Qnm1 = Qabs
     endif
     Qabs = (1.-albedo)*(Qdirect+Qrefl+Qscat*skyview)+emiss*(QIR+QIRre+Qlw*skyview)  ! Q absorbed


     do i=2,NSx-1
        do j=2,NSy-1
           if (h(i,j)<-32000) cycle
           if (subsurface) then ! subsurface conduction              
              call subsurfaceconduction_mars(T(:,i,j),Tsurf(i,j), &
                   & dtsec,Qnm1(i,j),Qabs(i,j),m(i,j),Fsurf(i,j),.false., &
                   & Tco2frost=Tco2frost,emiss=emiss)
           else  ! no subsurface conduction
              call equilibrT_mars(Tsurf(i,j),dtsec,Qnm1(i,j),Qabs(i,j),m(i,j),Tco2frost,emiss)
           endif
        enddo
     enddo

     where (Tsurf>Tco2frost .or. m<=0.) 
        albedo = albedo0
     elsewhere
        albedo = co2albedo
     end where
     
     ! only diagnostics below this line
     if (sdays > tmax-solsy) then
        Qmean(:,:) = Qmean(:,:) + Qabs
        where (Qabs>Qmax) Qmax=Qabs
        Tmean = Tmean + Tsurf
        where (Tsurf>Tmaxi) Tmaxi=Tsurf
        where (m>mmax) mmax=m
        where (m<mmin) mmin=m ! >0 if there is CO2 growth all year

        ! first time 273K is reached in southern spring
        where (meltfirst<0. .and. Tsurf>273.)
           meltfirst = marsLs
           meltfirst_sols = sdays
        end where
        ! this only works in the southern hemisphere

        if (subsurface) Tbottom(:,:) = Tbottom(:,:)+T(nz,:,:)
        nm=nm+1

        write(25,225) 2,2,sdays,marsLs/d2r,Qdirect(2,2),Qabs(2,2),Tsurf(2,2),m(2,2) ! flat
        i=41
        do j=46,49
           write(25,225) i,j,sdays,marsLs/d2r,Qdirect(i,j),Qabs(i,j),Tsurf(i,j),m(i,j)
        end do

     endif
225  format (2(i3,1x),f8.2,1x,f7.3,2(1x,f6.1),1x,f5.1,1x,f6.1)
     
     if (sdays > tmax-2*solsy) then  ! longest continuous period below H2O frost point
        where (Tsurf(:,:)<Tfrost)
           frosttime=frosttime+dt
        elsewhere
           frosttime=0.
        end where
        where (frosttime>maxfrosttime)
           maxfrosttime = frosttime
           h2olast = marsLs ! last time maxfrosttime increased
           h2olast_sols = sdays
        end where
     endif

  enddo  ! end of time loop

  close(25)
  if (subsurface) deallocate(T)

  Qmean = Qmean/nm
  Tmean = Tmean/nm
  where (h2olast/=-9.) h2olast=h2olast/d2r
  where (meltfirst/=-9.) meltfirst=meltfirst/d2r
  
  where (h2olast_sols < tmax-solsy) h2olast_sols = h2olast_sols+solsy
  where (h2olast_sols/=-9. .and. meltfirst_sols/=-9.) del = meltfirst_sols-h2olast_sols
  where (maxfrosttime<1) del=-9.

  open(unit=21,file='qmean.dat',status='unknown',action='write')
  do i=2,NSx-1
     do j=2,NSy-1
        write(21,'(2(i4,1x),f9.2,2x,f6.3,2(1x,f6.1),2(1x,f5.1),3(1x,f7.1),2(1x,f6.2),1x,f6.1)') &
             & i,j,h(i,j),surfaceSlope(i,j),Qmean(i,j),Qmax(i,j), &
             & Tmean(i,j),Tmaxi(i,j),mmax(i,j),maxfrosttime(i,j),mmin(i,j), &
             & h2olast(i,j),meltfirst(i,j),del(i,j)
     enddo
  enddo
  close(21)

  if (subsurface) then
     open(unit=23,file='tsurfbot.dat',status='unknown',action='write')
     Tbottom=Tbottom/nm
     do i=2,NSx-1
        do j=2,NSy-1
           write(23,'(2(i4,1x),2(1x,f7.1))') i,j,Tmean(i,j),Tbottom(i,j)
        enddo
     enddo
     close(23)
  endif
  
END PROGRAM cratersQ_mars
