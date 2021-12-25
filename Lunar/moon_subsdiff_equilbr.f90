program moon_subsdiff_equilibrium
  ! contiuum model for equilibrium subsurface adsorbate concentrations
  ! match subsurface desorption rates with time-averaged surface desorption rate
  implicit none
  integer, parameter :: MAXDEPTHIDX=500
  integer, parameter :: nz = 40
  real(8), parameter :: dz = 0.01  ! [m]
  integer i
  real(8) Esurf, Einfall
  real(8), dimension(0:nz) :: z, avEice, Tav, vvme, Tmin, Tmax
  real(8), external :: surface_av, bet_isotherm

  !Tav(:) = 0.
  !print *,'Domain depth = ',nz*dz
  !print *,'Mean free path and grain size = ',dz
  print *,'Number of vertical grid levels = ',nz

  Einfall = 100. ! kg/m^2/Gyr
  !Einfall = 1000.
  write(*,*) 'Einfall=',Einfall,'kg/m^2/Gyr'
  Esurf = surface_av(Einfall)
  write(*,'(a,2(1x,g10.4))') 'Esurf=', &
       & Esurf,Esurf*1.e9*365.2*86400 / 930. * 18*1.66e-27

  !call time_av_subl_fromdata('viper',nz,z,avEice,Tav)
  !call time_av_subl_fromdata('moonranger96',nz,z,avEice,Tav,Tmin,Tmax)

  open(22,file='profile.dat',action='write')
  
  do i=0,nz
     z(i) = i*dz
     call time_av_subl(z(i), avEice(i), Tav(i), Tmin(i), Tmax(i) )
     
     if (avEice(i) > Esurf) then
        vvme(i) = bet_isotherm( Esurf/avEice(i) )
     else
        vvme(i) = 999.
     end if
     write(22,'(f7.5,1x,f6.2,1x,g10.4,1x,g10.4,2(1x,f6.2))') &
          & z(i), Tav(i), Esurf/avEice(i), vvme(i), Tmin(i), Tmax(i)
  enddo
  close(22)
  print *,maxval(vvme)
end program moon_subsdiff_equilibrium



function surface_av(infallrate)
  ! time-averaged surface desorption rate = infall rate
  implicit none
  real(8) surface_av
  real(8), intent(IN) :: infallrate ! [kg/m^2/Ga]

  surface_av = infallrate/(1e9*365.2*86400)  ! Gyr -> s
  surface_av = surface_av /(18*1.66e-27)  ! -> #molecules/m^2/s
end function surface_av



subroutine time_av_subl(z,avE,Tav,Tmin,Tmax)
  ! time average of sublimation rate of ice using theoretical temperature profiles
  implicit none
  real(8), intent(IN) :: z
  real(8), intent(OUT) :: avE, Tav, Tmin, Tmax
  integer, parameter :: STEPSPERSOL = 360
  real(8), parameter :: pi = 3.1415926535, lunation = 29.52*86400
  real(8), parameter :: w = 2*pi/lunation
  real(8), parameter :: delta = 0.05  ! thermal skin depth
  !real(8), parameter :: g=0.  ! geothermal gradient [K/m]  
  real(8), parameter :: g=18.  ! geothermal gradient [K/m]  g = 0.018/0.001
  real(8), parameter :: Tm=113., Ta=40.
  integer j
  real(8) dt, time, T
  real(8), external :: sublrate1
  
  Tav = 0.
  avE = 0.
  dt = lunation/STEPSPERSOL
  Tmin = +999.; Tmax = -9.
  
  do j=1,STEPSPERSOL
     time = j*dt
     T = Tm + Ta*exp(-z/delta)*sin(z/delta-w*time) + g*z
     avE = avE + sublrate1(T)
     Tav = Tav + T
     if (T<Tmin) Tmin=T
     if (T>Tmax) Tmax=T
  enddo

  Tav = Tav / STEPSPERSOL
  avE = avE / STEPSPERSOL
end subroutine time_av_subl



subroutine time_av_subl_fromdata(filename,nz,z,avE,Tav,Tmax,Tmin)
  ! time average of sublimation rate of ice using temperature profiles from file
  implicit none
  character(len=*), intent(IN) :: filename
  integer, intent(IN) :: nz
  real(8), dimension(0:nz), intent(OUT) :: z, Tav, avE, Tmax, Tmin
  integer, parameter :: STEPSPERSOL = 96
  integer j, i, ierr
  real(8) T(0:nz)
  real(8), external :: sublrate1
  
  Tav(:) = 0.
  avE(:) = 0.
  z(:) = 0.
  Tmin(:) = +999.; Tmax(:) = -9.
  
  open(40,file='../DivinerMaps/z_Tprofiles.dat',action='read',iostat=ierr)
  open(41,file='../DivinerMaps/Tprofiles_'//filename//'.dat', &
       & action='read',iostat=ierr)
  if (ierr>0) stop 'file not found'
  
  read(40,*) z(0:nz)
  close(40)
  
  do j=1,STEPSPERSOL
     read(41,*) T
     do i=0,nz
        avE(i) = avE(i) + sublrate1(T(i))
        Tav(i) = Tav(i) + T(i)
        if (T(i)<Tmin(i)) Tmin(i)=T(i)
        if (T(i)>Tmax(i)) Tmax(i)=T(i)
     end do
  enddo
  close(41)
  
  Tav = Tav / STEPSPERSOL
  avE = avE / STEPSPERSOL
end subroutine time_av_subl_fromdata



real(8) function sublrate1(T)
  implicit none
  real(8), intent(IN) :: T   ! temperature [Kelvin]
  real(8), parameter :: pi=3.1415926535, kB=1.38065e-23, m=18*1.66e-27
  real(8) p, alphaT
  real(8), external :: psv, alpha
  
  p = psv(T)
  alphaT = alpha(T)
  sublrate1 = alphaT * p / sqrt(2*pi*m*kB*T)  ! [#molecules / m^2 / s]
  !sublrate2 = sublrate * padsr(vvm)  ! desorption rate
  !sublrate2 = sublrate2 * p_BET(vvm)
end function sublrate1

