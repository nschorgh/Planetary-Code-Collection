module thermalmodelparam_mars
  ! parameters for thermal model
  ! they are only used in the subroutines below
  real(8), parameter :: dt = 0.02  ! in units of Mars solar days
  !real(8), parameter :: Fgeotherm = 0.
  real(8), parameter :: Fgeotherm = 0.028  ! [W/m^2]
  real(8), parameter :: Lco2frost=6.0e5, co2albedo=0.60, co2emiss=1.
  real(8), parameter :: emiss0 = 1.  ! emissivity of dry surface
  integer, parameter :: EQUILTIME = 15 ! [Mars years]
end module thermalmodelparam_mars


!*****************************************************
! Subroutines for fast method
! written by Norbert Schorghofer 2007-2011
!*****************************************************


subroutine icelayer_mars(bigstep,nz,NP,thIn,rhoc,z,porosity,pfrost, &
     & Tb,zdepthF,zdepthE,porefill,Tmean1,Tmean3,zdepthG, &
     & latitude,albedo,p0,ecc,omega,eps,icefrac,zdepthT,avrho1, &
     & avrho1prescribed)
!*************************************************************************
! bigstep = time step [Earth years]
! latitude  [degree]
!*************************************************************************
  use miscparameters, only : d2r, NMAX, icedensity
  use allinterfaces, except_this_one => icelayer_mars
  !use omp_lib
  implicit none
  integer, intent(IN) :: nz, NP
  real(8), intent(IN) :: bigstep
  real(8), intent(IN) :: thIn(NP), rhoc(NP), z(NMAX), porosity, pfrost(NP)
  real(8), intent(INOUT) :: Tb(NP), porefill(nz,NP), zdepthF(NP), zdepthT(NP)
  real(8), intent(OUT), dimension(NP) :: zdepthE, Tmean1, Tmean3, zdepthG
  real(8), intent(IN), dimension(NP) :: latitude, albedo, p0
  real(8), intent(IN) :: ecc, omega, eps, icefrac
  real(8), intent(OUT) :: avrho1(NP)
  real(8), intent(IN), optional :: avrho1prescribed(NP)  ! <0 means absent

  integer k, typeF, typeG, typeT, j, jump, typeP
  real(8) fracIR, fracDust, ti(NMAX), rhocv(NMAX)
  real(8) Diff, ypp(nz), avdrho(NP), avdrhoP(NP), B, deltaz
  real(8), SAVE :: avdrho_old(100), zdepth_old(100)  ! NP<=100
  logical mode2

  !$omp parallel &
  !$omp    private(Diff,fracIR,fracDust,B,typeT,j,ti,rhocv,typeF,jump,typeG)
  !$omp do
  do k=1,NP   ! big loop

     Diff = 4e-4*600./p0(k)
     fracIR = 0.04*p0(k)/600.; fracDust = 0.02*p0(k)/600.
     B = Diff*bigstep*86400.*365.24/(porosity*icedensity)

     typeT = -9
     if (zdepthT(k)>=0. .and. zdepthT(k)<z(nz)) then
        do j=1,nz
           if (z(j)>zdepthT(k)) then ! ice
              typeT = j  ! first point in ice
              exit
           endif
        enddo
     endif

     call assignthermalproperties(nz,thIn(k),rhoc(k),ti,rhocv, &
          & typeT,icefrac,porosity,porefill(:,k))
     
     !----run thermal model
     call ajsub_mars(typeT, latitude(k)*d2r, albedo(k), pfrost(k), nz, z, & 
          &     ti, rhocv, fracIR, fracDust, p0(k), ecc, omega, eps, &
          &     avdrho(k), avdrhoP(k), avrho1(k), Tb(k), zdepthE(k), typeF, &
          &     zdepthF(k), ypp, porefill(:,k), Tmean1(k), Tmean3(k), B, &
          &     typeG, zdepthG(k), avrho1prescribed(k))

     if (typeF*zdepthF(k)<0.) stop 'error in icelayer_mars'
     ! diagnose
     if (zdepthT(k)>=0.) then
        jump = 0
        do j=1,nz
           if (zdepth_old(k)<z(j) .and. zdepthT(k)>z(j)) jump = jump+1
        enddo
     else
        jump = -9
     endif
     if (zdepthT(k)>=0. .and. avdrho(k)*avdrho_old(k)<0.) then 
        write(34,*) '# zdepth arrested'
        if (jump>1) write(34,*) '# previous step was too large',jump
     endif
     write(34,'(f8.3,1x,f6.2,1x,f11.5,2(1x,g11.4),1x,i3)') &
          &        bigstep,latitude(k),zdepthT(k),avdrho(k),avdrhoP(k),jump
     zdepth_old(k) = zdepthT(k)
     avdrho_old(k) = avdrho(k)

!----mode 2 growth  
     typeP = -9;  mode2 = .FALSE.
     do j=1,nz
        if (porefill(j,k)>0.) then
           typeP = j  ! first point with ice
           exit
        endif
     enddo
     if (typeT>0 .and. typeP>2 .and. zdepthE(k)>0.) then
        if (porefill(typeP,k)>=1. .and. porefill(typeP-1,k)==0. .and. &
             & zdepthE(k)<z(typeP) .and. &
             & z(typeP)-zdepthE(k)>2*(z(typeP)-z(typeP-1))) then  ! trick that avoids oscillations
           deltaz = -avdrhoP(k)/z(typeP)*18./8314.*B  ! conservation of mass 
           if (deltaz>z(typeP)-z(typeP-1)) then  ! also implies avdrhoP<0.
              mode2 = .TRUE.
           endif
        endif
     endif

     !call icechanges_poreonly(nz,z,typeF,typeG,avdrhoP(k),ypp,B,porefill(:,k))
     call icechanges(nz,z(:),typeF,avdrho(k),avdrhoP(k),ypp(:), &
          & Diff,porosity,icefrac,bigstep,zdepthT(k),porefill(:,k),typeG)

     if (mode2 .and. porefill(typeP,k)>=1. .and. porefill(typeP-1,k)==0.) then  ! nothing changed
        porefill(typeP-1,k)=1.  ! paste a layer
        write(34,*) '# mode 2 growth occurred',typeP,typeF,typeT
     endif

  enddo  ! end of big loop
  !$omp end do
  !$omp end parallel
end subroutine icelayer_mars



subroutine ajsub_mars(typeT, latitude, albedo0, pfrost, nz, z, ti, rhocv, &
     &     fracIR, fracDust, patm, ecc, omega, eps, avdrho, avdrhoP, avrho1, &
     &     Tb, zdepthE, typeF, zdepthF, ypp, porefill, Tmean1, Tmean3, &
     &     B, typeG, zdepthG, avrho1prescribed)
!***********************************************************************
!  A 1D thermal model that returns various averaged quantities
!
!  Tb<0 initializes temperatures
!  Tb>0 initializes temperatures with Tb 
!***********************************************************************
  use miscparameters
  use thermalmodelparam_mars
  use allinterfaces, except_this_one => ajsub_mars
  implicit none
  integer, intent(IN) :: nz, typeT
  real(8), intent(IN) :: latitude  ! in radians
  real(8), intent(IN) :: albedo0, pfrost, z(NMAX)
  real(8), intent(IN) :: ti(NMAX), rhocv(NMAX), fracIR, fracDust, patm
  real(8), intent(IN) :: ecc, omega, eps, porefill(nz)
  real(8), intent(OUT) :: avdrho, avdrhoP  ! difference in annual mean vapor density
  real(8), intent(OUT) :: avrho1  ! mean vapor density on surface
  real(8), intent(INOUT) :: Tb, Tmean1
  integer, intent(OUT) :: typeF  ! index of depth below which filling occurs
  real(8), intent(OUT) :: zdepthE, zdepthF 
  real(8), intent(OUT) :: ypp(nz) ! (d rho/dt)/Diff
  real(8), intent(OUT) :: Tmean3, zdepthG
  real(8), intent(IN) :: B  ! just passing through
  integer, intent(OUT) :: typeG
  real(8), intent(IN), optional :: avrho1prescribed
  real(8), parameter :: a = 1.52366 ! Mars semimajor axis in A.U.
  integer nsteps, n, i, nm, typeP
  real(8) tmax, time, Qn, Qnp1, tdays
  real(8) marsR, marsLs, marsDec, HA
  real(8) Tsurf, Tco2frost, albedo, Fsurf, m, dE, emiss, T(NMAX)
  real(8) Told(nz), Fsurfold, Tsurfold, Tmean0, avrho2
  real(8) rhosatav0, rhosatav(nz), rlow 
  real(8), external :: psv, tfrostco2
  
  tmax = EQUILTIME*solsperyear
  nsteps = int(tmax/dt)     ! calculate total number of timesteps

  Tco2frost = tfrostco2(patm) 

  if (Tb<=0.) then  ! initialize
     !Tmean0 = 210.15       ! black-body temperature of planet
     Tmean0 = (589.*(1.-albedo0)*cos(latitude)/(pi*emiss0*sigSB))**0.25 ! better estimate
     Tmean0 = Tmean0-5.
     write(34,*) '# initialized with temperature estimate at',latitude/d2r,'of',Tmean0,'K'
     T(1:nz) = Tmean0 
  else
     T(1:nz) = Tb
     ! not so good when Fgeotherm is on
  endif
  
  albedo = albedo0
  emiss = emiss0
  do i=1,nz
     if (T(i)<Tco2frost) T(i)=Tco2frost
  enddo
  Tsurf = T(1)
  m=0.; Fsurf=0.

  nm=0
  avrho1=0.; avrho2=0.
  Tmean1=0.; Tmean3=0.
  rhosatav0 = 0.
  rhosatav(:) = 0.

  time=0.
  call generalorbit(0.d0,a,ecc,omega,eps,marsLs,marsDec,marsR)
  HA = 2.*pi*time            ! hour angle
!  Qn=flux(marsR,marsDec,latitude,HA,albedo,fracir,fracdust,0.d0,0.d0)
  Qn = flux_mars77(marsR,marsDec,latitude,HA,albedo,fracir,fracdust)
  !----loop over time steps 
  do n=0,nsteps-1
     time = (n+1)*dt         !   time at n+1 
     tdays = time*(marsDay/earthDay) ! parenthesis may improve roundoff
     call generalorbit(tdays,a,ecc,omega,eps,marsLs,marsDec,marsR)
     HA = 2.*pi*mod(time,1.d0)  ! hour angle
!     Qnp1=flux(marsR,marsDec,latitude,HA,albedo,fracir,fracdust,0.d0,0.d0)
     Qnp1 = flux_mars77(marsR,marsDec,latitude,HA,albedo,fracir,fracdust)
     
     Tsurfold = Tsurf
     Fsurfold = Fsurf
     Told(1:nz) = T(1:nz)
     if (m<=0. .or. Tsurf>Tco2frost) then
        call conductionQ(nz,z,dt*marsDay,Qn,Qnp1,T,ti,rhocv,emiss, &
             &           Tsurf,Fgeotherm,Fsurf)
     endif
     if (Tsurf<Tco2frost .or. m>0.) then ! CO2 condensation
        T(1:nz) = Told(1:nz)
        call conductionT(nz,z,dt*marsDay,T,Tsurfold,Tco2frost,ti, &
             &              rhocv,Fgeotherm,Fsurf) 
        Tsurf = Tco2frost
        dE = (- Qn - Qnp1 + Fsurfold + Fsurf + &
             &           emiss*sigSB*(Tsurfold**4+Tsurf**4))/2.
        m = m + dt*marsDay*dE/Lco2frost
     endif
     if (Tsurf>Tco2frost .or. m<=0.) then
        albedo = albedo0
        emiss = emiss0
     else
        albedo = co2albedo
        emiss = co2emiss
     endif
     Qn=Qnp1
     
     if (time>=tmax-solsperyear) then
        Tmean1 = Tmean1 + Tsurf
        Tmean3 = Tmean3 + T(nz)
        avrho1 = avrho1 + min(psv(Tsurf),pfrost)/Tsurf
        rhosatav0 = rhosatav0 + psv(Tsurf)/Tsurf
        do i=1,nz
           rhosatav(i) = rhosatav(i) + psv(T(i))/T(i)
        enddo
        nm = nm+1
     endif

  enddo  ! end of time loop
  
  avrho1 = avrho1/nm
  if (present(avrho1prescribed)) then
     if (avrho1prescribed>=0.) avrho1=avrho1prescribed
  endif
  Tmean1 = Tmean1/nm; Tmean3 = Tmean3/nm
  rhosatav0 = rhosatav0/nm; rhosatav(:)=rhosatav(:)/nm
  if (typeT>0) then
     avrho2 = rhosatav(typeT)
  else
     avrho2 = rhosatav(nz) ! no ice
  endif
  avdrho = avrho2-avrho1
  typeP = -9 
  do i=1,nz
     if (porefill(i)>0.) then
        typeP = i  ! first point with ice
        exit
     endif
  enddo
  if (typeP>0) then
     avdrhoP = rhosatav(typeP) - avrho1
  else  
     avdrhoP = -9999.
  end if

  zdepthE = equildepth(nz, z, rhosatav, rhosatav0, avrho1)

  if (Fgeotherm>0.) then
     Tb = Tmean1 
     typeG = 1   ! will be overwritten by depths_avmeth
     rlow = 2*rhosatav(nz)-rhosatav(nz-1)
  else
     Tb = T(nz)
     typeG = -9
     rlow = rhosatav(nz-1)
  endif
  call depths_avmeth(typeT, nz, z, rhosatav(:), rhosatav0, rlow, avrho1,  &
       & porefill(:), typeF, zdepthF, B, ypp(:), typeG, zdepthG)

end subroutine ajsub_mars



subroutine outputmoduleparameters
  use miscparameters
  use thermalmodelparam_mars
  implicit none
  print *,'Parameters stored in modules'
  print *,'  Ice bulk density',icedensity,'kg/m^3'
  print *,'  dt=',dt,'Mars solar days'
  print *,'  Fgeotherm=',Fgeotherm,'W/m^2'
  write(6,'(2x,a27,1x,f5.3)') 'Emissivity of dry surface=',emiss0
  write(6,'(2x,a12,1x,f5.3,2x,a16,1x,f5.3)') 'CO2 albedo=',co2albedo,'CO2 emissivity=',co2emiss
  print *,'  Thermal model equilibration time',EQUILTIME,'Mars years'
end subroutine outputmoduleparameters



