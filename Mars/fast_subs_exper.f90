module thermalmodelparam_exper
  ! parameters for thermal model
  real(8), parameter :: dt = 0.4  ! in units of Mars solar days
  real(8), parameter :: Fgeotherm = 0.
  !real(8), parameter :: Fgeotherm = 0.028  ! W/m^2
  integer, parameter :: EQUILTIME = 10 ! (mars years)
end module thermalmodelparam_exper


!*****************************************************
! Subroutines for fast method
! applicable when surface temperature is prescribed
! written by Norbert Schorghofer 2007-2009
!*****************************************************


subroutine icelayer_exper(bigstep, nz, thIn, rhoc, z, porosity, &
     & pfrost, zdepthT, icefrac, zdepthF, zdepthE, porefill,  &
     & Tmean, Tampl, Diff, zdepthG)
!*************************************************************************
! bigstep = time step [Earth years]
!*************************************************************************
  use miscparameters, only : d2r, NMAX, icedensity
  use allinterfaces, except_this_one => icelayer_exper
  implicit none
  integer, intent(IN) :: nz
  real(8), intent(IN) :: bigstep
  real(8), intent(IN) :: thIn, rhoc, z(NMAX), porosity, pfrost
  real(8), intent(INOUT) :: zdepthT, zdepthF, porefill(nz)
  real(8), intent(OUT) :: zdepthE
  real(8), intent(IN) :: icefrac, Diff, Tmean, Tampl
  real(8), intent(OUT) :: zdepthG

  integer j, typeT, typeF, typeG
  real(8) ti(NMAX), rhocv(NMAX), ypp(nz), avdrho, avdrhoP, B
  real(8), parameter :: NULL = 0.
  ! real(8), external :: colint  ! defined in allinterfaces.mod

  B = Diff*bigstep*86400.*365.24/(porosity*icedensity)

  typeT = -9
  if (zdepthT>=0. .and. zdepthT<z(nz)) then
     do j=1,nz
        if (z(j)>zdepthT) then ! ice
           typeT = j  ! first point in ice
           exit
        endif
     enddo
  endif

  call assignthermalproperties(nz,thIn,rhoc,ti,rhocv, &
       & typeT,icefrac,porosity,porefill)
  
  !----run thermal model
  call ajsub_exper(typeT, nz, z, ti, rhocv, pfrost, Tmean, Tampl, &
       &     avdrho, avdrhoP, zdepthE, typeF, zdepthF, &
       &     ypp, porefill, B, typeG, zdepthG)
  
  if (typeF*zdepthF<0.) stop 'error in icelayer_poreonly'
  
  !call icechanges_poreonly(nz,z(:),typeF,typeG,avdrhoP,ypp(:),B,porefill(:))
  call icechanges(nz,z(:),typeF,avdrho,avdrhoP,ypp(:), &
          & Diff,porosity,icefrac,bigstep,zdepthT,porefill(:),typeG)

end subroutine icelayer_exper



subroutine ajsub_exper(typeT, nz, z, ti, rhocv, pfrost, Tmean, Tampl, &
     &     avdrho, avdrhoP, zdepthE, typeF, zdepthF, ypp, porefill, &
     &     B, typeG, zdepthG)
!***********************************************************************
!  A 1D thermal model that returns various averaged quantities
!  uses prescribed surface temperatures
!
!  a generalization of jsub2
!***********************************************************************
  use miscparameters, only : NMAX, solsperyear, marsDay
  use thermalmodelparam_exper
  use allinterfaces, except_this_one => ajsub_exper
  implicit none
  integer, intent(IN) :: nz, typeT
  real(8), intent(IN) :: z(NMAX), ti(NMAX), rhocv(NMAX), pfrost
  real(8), intent(IN) :: Tmean, Tampl
  real(8), intent(OUT) :: avdrho, avdrhoP  ! difference in annual mean vapor density
  real(8), intent(OUT) :: zdepthE
  integer, intent(OUT) :: typeF  ! index of depth below which filling occurs
  real(8), intent(INOUT) :: zdepthF 
  real(8), intent(OUT) :: ypp(nz) ! (d rho/dt)/Diff
  real(8), intent(IN) :: porefill(nz)
  real(8), intent(IN) :: B  ! just passing through
  integer, intent(OUT) :: typeG
  real(8), intent(OUT) :: zdepthG

  integer nsteps, n, i, nm, typeP
  real(8) tmax, time, Tsurf, Tsurfp1, Fsurf, T(NMAX)
  real(8) avrho1, avrho2, rhosatav0, rhosatav(nz), rlow
  real(8), external :: Tsurface, psv
  ! real(8), external :: equildepth   ! defined in allinterfaces.mod
  
  tmax = EQUILTIME*solsperyear
  nsteps=int(tmax/dt)       ! calculate total number of timesteps

  T(1:nz) = Tmean
  Tsurf = T(1)
  
  nm=0
  avrho1=0.; avrho2=0.
  rhosatav0 = 0.
  rhosatav(:) = 0.

  time=0.
  !-----loop over time steps 
  do n=0,nsteps-1
     time =(n+1)*dt         !   time at n+1 
     Tsurfp1 = Tsurface(time,Tmean,Tampl)
     call conductionT(nz,z,dt*marsDay,T,Tsurf,Tsurfp1,ti, &
             &              rhocv,Fgeotherm,Fsurf) 
     Tsurf=Tsurfp1
          
     if (time>=tmax-solsperyear) then
        avrho1 = avrho1+min(psv(Tsurf),pfrost)/Tsurf
        rhosatav0 = rhosatav0+psv(Tsurf)/Tsurf
        do i=1,nz
           rhosatav(i) = rhosatav(i)+psv(T(i))/T(i)
        enddo
        nm=nm+1
     endif

  enddo  ! end of time loop
  
  avrho1 = avrho1/nm
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
  endif

  zdepthE = equildepth(nz, z, rhosatav, rhosatav0, avrho1)

  if (Fgeotherm>0.) then
     typeG=1   ! will be overwritten by avmeth
     rlow=2*rhosatav(nz)-rhosatav(nz-1)
  else
     typeG=-9
     rlow=rhosatav(nz-1)
  endif
  call depths_avmeth(typeT, nz, z(:), rhosatav(:), rhosatav0, rlow, avrho1,  &
       & porefill(:), typeF, zdepthF, B, ypp(:), typeG, zdepthG)

end subroutine ajsub_exper




pure function Tsurface(time,Tmean,Tampl)
  ! sinusoidal surface temperature
  use miscparameters
  implicit none
  real(8) Tsurface
  real(8), intent(IN) :: time, Tmean, Tampl
  real(8), parameter :: omegaMars=2.*pi/(solsperyear)
  Tsurface=Tmean + Tampl*sin(omegaMars*time)
  return
end function Tsurface



subroutine outputmoduleparameters
  use miscparameters
  use thermalmodelparam_exper
  implicit none
  print *,'Parameters stored in modules'
  print *,'  Ice bulk density',icedensity,'kg/m^3'
  print *,'  dt=',dt,'Mars solar days'
  print *,'  Fgeotherm=',Fgeotherm,'W/m^2'
  print *,'  Thermal model equiltime',EQUILTIME,'Mars years'
end subroutine outputmoduleparameters







