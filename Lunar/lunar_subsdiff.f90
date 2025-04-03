program lunar_subsdiff
  ! subsurface diffusion of H2O on airless body
  ! used in upcoming publication, Schorghofer (2025)
  use params
  implicit none
  integer i, j
  real(8) time ! [s]
  real(8), dimension(0:NZ) :: T  ! [K]
  real(8) theta(0:NZ)   ! [#molecules/m^2]
  real(8) Sinfall, weath, Mtotal, outintrvl, theta0wo, Mprovided, S(0:NZ)
  real(8), external :: desorptionrate, desorptionrate_ice
  real(8), external :: weathering, colintmass, solver_wodiffusion

  call output_module_params

  open(22,file='Tprofiles.dat',action='write')
  open(23,file='aprofiles.dat',action='write')
  open(24,file='longseries.dat',action='write')
  open(26,file='TprofilesC.dat',action='write')  ! last cycle
  open(27,file='aprofilesC.dat',action='write')  ! last cycle
  open(28,file='Sprofiles.dat',action='write')

  theta(:) = 0.*thetaML   ! initial adsorption profile
  !theta = 0.5*thetaML
  Mprovided = 0.
  
  do i=1,ceiling(maxtime/dtsec)  ! begin time loop
     time = i*dtsec

     call temperatureprofile(NZ,time,T,Deltaz)

     !weath = 0.
     weath = weathering(theta(0),time)

     Sinfall = 0.
     !Sinfall = 1e19/(365.24*86400) *3
     Mprovided = Mprovided + dtsec*Sinfall/thetaML
     
     call solver(T,time,weath,Sinfall,theta)

     theta0wo = solver_wodiffusion(T(0),weath,Sinfall,theta(0))
     
     ! output time series
     if (mod(i,1000)==0) then
     !if (mod(i,STEPSPERSOL/2)==0) then
        Mtotal=colintmass(theta)
        write(*,('(f0.5,1x,f0.1,2(1x,g0.4),*(1x,f0.4))')) &
             & time/secyear,T(0),Mtotal,Mprovided,theta(0)/thetaML,theta0wo/thetaML
        write(24,('(f0.5,1x,g0.4,*(1x,f0.4))')) &
             & time/secyear,Mtotal,theta(0)/thetaML,theta(nz)/thetaML,theta0wo/thetaML
     endif

     ! output depth profiles
     outintrvl = 1000.*secyear  ! e.g. at 130K, static profile
     !outintrvl = 100.*secyear  ! e.g. 250+/-100K
     if (mod(time,outintrvl)<=dtsec/2. .or. mod(time,outintrvl)>outintrvl-dtsec/2.) then
        write(22,'(f0.4,*(1x,f0.3))') time/secyear, T(:)
        write(23,'(f0.4,*(1x,f0.5))') time/secyear, theta(:)/thetaML
        do j=0,nz
           S(j) = desorptionrate(T(j),theta(j))
        end do
        write(28,'(f0.4,*(1x,g0.5))') time/secyear, S(:)
     endif
     if (i> ceiling(maxtime/dtsec)-STEPSPERSOL) then
        write(26,'(f0.5,*(1x,f0.3))') time/secyear, T(:)
        write(27,'(f0.5,1x,g0.4,*(1x,f0.5))') time/secyear, colintmass(theta), theta(:)/thetaML
     endif
     
  enddo
  
  print *,'Total time',time,time/secyear
  close(22); close(23); close(24)
  close(26); close(27); close(28)

  open(25,file='lastprofile.dat',action='write')
  do j=0,NZ
     write(25,'(f5.3,1x,f6.2,1x,f0.2,1x,g0.4)') &
          & Deltaz*j,T(j),theta(j)/thetaML,desorptionrate(T(j),theta(j))
  end do
  write(25,*) Deltaz*(NZ+1),'Inf',desorptionrate_ice(T(NZ))
  close(25)
end program lunar_subsdiff


function weathering(theta0,time)
  ! ice loss due to space weathering
  use params, only : pi, lunation, thetaML, mu, Yrough, SSA
  implicit none
  real(8), intent(IN) :: theta0  ! adsorbate concentration on surface [#/m^2]
  real(8), intent(IN) :: time
  real(8) weathering, w
  real(8), parameter :: wsun = 1.16e14 ! space weathering (sunlit)
  real(8), parameter :: wLy = 7e10 ! Ly-alpha, Morgan & Shemansky (1991)
  real(8), parameter :: wI = 1e14 ! space weathering (impacts)
  real(8), parameter :: dLat = 85. ! latitude in degrees (for solar UV)

  w = 0.

  ! toggle the contributions on/off
  w = wsun/Yrough * min(theta0/thetaML,1.) ! solar photo-destruction, zenith
  ! incidence angle effect
  !w = w * max( cosd(dLat)*sin(-2*pi*time/lunation), 0.)  ! *cos(incidence angle)
  w = w * cosd(dLat) / pi  ! or use time-averaged value instead
  
  !w = w + wLy/Yrough * min(theta0/thetaML,1.)  ! all-sky Lyman-alpha
  !w = w + wI * min(SSA*mu*theta0,1.)  ! destruction by impacts
  weathering = w
end function weathering


function colintmass(theta)
  ! column integrated H2O mass
  use params, only : mu, lambda, Yrough, nz, Deltaz
  implicit none
  real(8) colintmass
  real(8), intent(IN) :: theta(0:nz)
  colintmass = 2*Yrough/lambda*mu * sum(theta) * Deltaz
end function colintmass


subroutine solver(T,time,w,Sinfall,theta)
  ! time step for H2O adsorbate density
  ! surface has depth-index 0
  use params, only : thetaML, lambda, Yrough, porosity, Deltaz, dtsec, nz
  implicit none
  real(8), intent(IN) :: T(0:nz), time
  real(8), intent(IN) :: w, Sinfall
  real(8), intent(INOUT) :: theta(0:nz)  ! [#molecules/m^2]
  real(8), parameter :: uconv = 86400*365.24/thetaML
  integer j
  real(8) alpha, SNp1, J0
  real(8), dimension(0:nz) :: S, Tnext
  real(8) S0next, theta0new
  real(8), external :: desorptionrate, desorptionrate_ice
  logical, parameter :: ICETABLE = .true.

  do j=0,nz
     S(j) = desorptionrate(T(j),theta(j))
  end do
  !write(*,'(f0.4,1x,f0.4,1x,f0.2,*(1x,g0.4))') theta(nz)/thetaML,T(nz),S(nz)
  
  ! Upper boundary condition
  ! Sinfall and w are input parameters
  J0 = - (porosity*lambda/Deltaz) * ( S(1)-S(0) ) ! upward flux is negative

  ! Population change on surface: simple Euler step
  theta0new = theta(0) + dtsec*(-S(0) + Sinfall - J0)/Yrough - dtsec*w

  !write(*,('(4(1x,g0.3))')) S(0)*uconv,w*uconv,J0*uconv
  
  if ((theta(0)>0.1 .and. (theta0new/theta(0)>1.2 .or. theta0new/theta(0)<0.8)) &
       & .or. theta0new<0.) then
     !print *,'Taking half step'
     call temperatureprofile(NZ,time+dtsec/2.,Tnext,Deltaz)
     S0next = desorptionrate(Tnext(0),theta(0))
     theta0new = theta(0) + dtsec*( -S0next + Sinfall - J0)/Yrough - dtsec*w
  endif

  theta(0) = theta0new
  if (theta(0)<0.) theta(0)=0.

  ! Interior
  alpha = dtsec * porosity/(2*Yrough) * (lambda/Deltaz)**2
  do j=1,nz-1
     theta(j) = theta(j) + alpha*(-2*S(j)+S(j+1)+S(j-1))
  end do

  ! Lower boundary condition
  if (ICETABLE) then
     SNp1 = desorptionrate_ice(T(nz))
  else
     SNp1 = S(nz-1)  ! impermeable (derivative vanishes)
  end if
  theta(nz) = theta(nz) + alpha*(-2*S(nz)+SNp1+S(nz-1))
  
  where (theta<0.) theta=0.
  
  !write(*,'(f0.2,1x,f0.4,3(1x,g0.4))') T(0),theta(0)/thetaML,Sinfall,w,S(0)
end subroutine solver


function solver_wodiffusion(T0,w,Sinfall,theta0)
  ! integrate H2O adsorbate density on surface without the subsurface
  use params, only : thetaML, Yrough, dtsec
  implicit none
  real(8), intent(IN) :: T0, w, Sinfall, theta0
  real(8) solver_wodiffusion, S0, theta0wo
  real(8), external :: desorptionrate

  S0 = desorptionrate(T0,theta0)
  
  ! Population change on surface: simple Euler step
  theta0wo = theta0 + dtsec*( -S0/Yrough - w + Sinfall/Yrough )

  if (theta0wo<0.) theta0wo=0.
  solver_wodiffusion = theta0wo
  
end function solver_wodiffusion

