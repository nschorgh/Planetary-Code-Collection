subroutine subsurfaceconduction_mars(T,Tsurf,dtsec,Qn,Qnp1,m,Fsurf,init, &
     & Tco2frost,thIn,emiss)
  use miscparams, only : pi, sigSB, solsy, Lco2frost, solarDay, nz
  use conductionQ
  use conductionT
  implicit none
  real(8), intent(INOUT) :: T(:)  ! in interface declared as T(:), not T(nz)
  real(8), intent(INOUT) :: Tsurf, m, Fsurf
  real(8), intent(IN) :: dtsec,Qn,Qnp1
  logical, intent(IN) :: init
  real(8), intent(IN), optional :: Tco2frost, thIn, emiss
  real(8), parameter :: Fgeotherm = 0.0 ! [W/m^2]
  integer i
  !real(8), parameter :: zmax=3., zfac=1.05d0  ! adjust
  real(8), parameter :: zmax=13., zfac=1.05d0  ! rhoc=thIn*1000 (nz=70, 3x seas)
  real(8) Tinit, delta
  real(8) Fsurfold, dE, Tsurfold, Told(nz)
  real(8) z(nz), ti(nz), rhocv(nz)

  if (init) then ! initialize grid
     if (.not.present(thIn)) then
        error stop 'missing argument'
     end if
     
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
     open(unit=30,file='z',action='write');
     write(30,*) (z(i),i=1,nz)
     close(30)

     write(*,*) 'Subsurface model parameters'
     write(*,*) '   nz=',nz,' zmax=',zmax,' zfac=',zfac
     write(*,*) '   Thermal inertia=',thIn,' rho*c=',rhocv(1)
     write(*,*) '   (Scaled) diurnal and seasonal skin depths=', &
          & delta,delta*sqrt(solsy)
     write(*,*) '   Geothermal flux=',Fgeotherm

     call conductionT2_init(nz,z,dtsec,ti,rhocv,Fgeotherm)
     call conductionQ2_init(nz,z,dtsec,ti,rhocv,Fgeotherm)
     
     return
  endif
  
  if (Tsurf<=0.) then  ! initialize temperature profile
     Tinit = 200.
     T(1:nz) = Tinit
     Tsurf = Tinit
  endif

  if (.not.present(Tco2frost) .or. .not.present(emiss)) then
     error stop 'missing argument'
  end if
  
  Tsurfold = Tsurf
  Fsurfold = Fsurf
  Told(1:nz) = T(1:nz)
  if (Tsurf>Tco2frost .or. m<=0.) then
     call conductionQ2(nz,Qn,Qnp1,T,emiss,Tsurf,Fsurf)
  endif
  if (Tsurf<Tco2frost .or. m>0.) then   ! CO2 condensation   
     T(1:nz) = Told
     call conductionT2(nz,Tsurfold,Tco2frost,T,Fsurf)
     Tsurf = Tco2frost
     dE = (- Qn - Qnp1 + Fsurfold + Fsurf + &
          &           emiss*sigSB*(Tsurfold**4+Tsurf**4))/2.
     m = m + dtsec*dE/Lco2frost
  endif

end subroutine subsurfaceconduction_mars



elemental subroutine equilibrT_mars(Tsurf, dtsec, Qnm1, Qn, m, Tco2frost, emiss)
  ! calculates equilibrium temperature and change in CO2 mass
  ! this is the analog to subsurfaceconduction_mars for thIn=0
  use miscparams, only : sigSB, Lco2frost
  implicit none
  real*8, intent(IN) :: Qn, Qnm1, Tco2frost, emiss, dtsec
  real*8, intent(OUT) :: Tsurf
  real*8, intent(INOUT) :: m
  !real(8), parameter :: sigSB = 5.6704e-8
  !real(8), parameter :: Lco2frost = 6.0e5 ! [J/kg]
  real*8 Tsurfold, dE
  
  Tsurf = (Qn/emiss/sigSB)**0.25
  Tsurfold = (Qnm1/emiss/sigSB)**0.25
  if (Tsurf<Tco2frost .or. m>0.) then   ! CO2 condensation
     Tsurf = Tco2frost
     dE = - Qn + emiss*sigSB*(Tsurf**4 + Tsurfold**4)/2.
     m = m + dtsec*dE/Lco2frost
  endif

end subroutine equilibrT_mars



pure function evap_ingersoll(T,p0)
  ! Returns sublimation rate of water ice into martian atmosphere [kg/m^2/s]
  ! Note: The latent heat of sublimation is 2.838 MJ/kg
  use allinterfaces, only : psv
  use, intrinsic :: ieee_arithmetic
  implicit none
  real(8) evap_ingersoll
  real(8), intent(IN) :: T
  real(8), intent(IN) :: p0  ! atmospheric pressure [Pa]
  real(8), parameter :: R=8314.5, g=3.7
  real(8) psat, Gbuf, rho, rhow, drhooverrho
  real(8) D   ! vapor diffusivity [m^2/s]
  real(8) nu  ! kinematic viscosity of CO2 [m^2/s]
  real(8) eta ! dynamic viscosity of CO2 [Pa.s]
  
  psat = psv(T)
  D = 0.1654*1e-4*1.013e5/p0*(T/273)**1.5  ! Schwertz & Brow (1951)
  rhow = psat*18/(R*T)
  rho = p0*44/(R*T)
  eta = 13.7e-6*(273+240)/(T+240)*(T/273)**1.5  ! Int. Crit. Tbl., vol 5
  nu = eta/rho
  
  !drhooverrho = (44-18)*psat/(44*p0-(44-18)*psat) ! Ingersoll (1970)
  drhooverrho = (44-18)*psat/(44*(p0-psat)) ! diverges at p0=psat
  if (drhooverrho>0.) then
     Gbuf = (drhooverrho*g/nu**2)**(1./3.)
     !evap_ingersoll = 0.17*D*rhow*Gbuf  ! Ingersoll (1970)
     evap_ingersoll = 0.12*D*rhow*Gbuf   ! Schorghofer (2020)
  else
     !evap_ingersoll = 1./0.  ! a smart compiler accepts this
     evap_ingersoll = ieee_value(1.d0,ieee_positive_inf)  ! for other compilers
  endif
  
end function evap_ingersoll
