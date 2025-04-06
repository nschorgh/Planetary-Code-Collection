module params
   implicit none
   
   ! Constants
   real(8), parameter :: thetaML = 1e19
   real(8), parameter :: mu = 18*1.66e-27 ! mass of H2O molecule
   real(8), parameter :: secyear = 365.24*86400, lunation = 29.52*86400.
   real(8), parameter :: pi = 3.1415926535897932d0

   ! Material parameters
   real(8), parameter :: lambda = 71e-6  ! grain spacing [m]
   real(8), parameter :: Yrough = 27  ! roughness factor
   real(8), parameter :: porosity = 0.4
   real(8), parameter :: SSA = 500. ! specific surface area [m^2/kg]
   ! Yrough = SSA*rhobar*lambda/2 must hold, rhobar~1500kg/m^3

   ! Numerical grid
   integer NZ      ! number of grid points
   real(8) Deltaz  ! grid spacing [m]
   real(8) maxtime ! run time [s]

   !real(8), parameter :: g=0.
   real(8), parameter :: g=9.  ! geothermal gradient [K/m]  g = 0.018/0.002
   real(8) Tm, Ta  ! mean temperature and temperature amplitude [K]

   ! Static situations
   parameter (Tm=130., Ta=0.)
   parameter (NZ = 200, Deltaz=0.5e-3, maxtime=(1e6+0.1)*secyear) ! 130K, static
   real(8), parameter :: dtsec = 10.*secyear

   ! Periodic situations
   integer STEPSPERSOL 
   !parameter (Tm=130., Ta=40.)
   !parameter (Tm=250., Ta=100.)
   !parameter (NZ = 200, Deltaz=0.5e-3, maxtime=(100e3+0.1)*secyear, STEPSPERSOL=24) ! 130+/-40K
   !parameter (NZ = 100, Deltaz=1e-5, maxtime=(1e3+0.1)*secyear, STEPSPERSOL=48) ! Sin3m 250+/-100K
   !real(8), parameter :: dtsec = lunation/STEPSPERSOL ! time step [s]
   
   ! adjust upper and lower boundary conditions manually
   ! adjust weathering rates manually
end module params


subroutine output_module_params()
  use params
  implicit none
  integer j
  open(20,file='z.dat',action='write')
  write(20,'(*(f0.7,1x))') (j*Deltaz, j=0,NZ)
  close(20)
  open(21,file='params.dat',action='write')
  write(21,*) 'Yrough=',Yrough,'SSA=',SSA
  write(21,*) 'lambda=',lambda,'porosity=',porosity
  write(21,*) 'Time step =',dtsec,'sec',dtsec/secyear,'years'
  write(21,*) 'Maximum run time =',maxtime/secyear,'years'
  write(21,*) 'Domain depth =',NZ*Deltaz
  write(21,*) 'Vertical resolution =',Deltaz
  write(21,*) 'Number of vertical grid levels =',NZ
  write(21,*) 'g=',g
  write(21,*) 'Tm=',Tm,'Ta=',Ta
  close(21)
end subroutine output_module_params


subroutine temperatureprofile(nz,time,T,dz)
  ! assign theoretical temperature profile
  use params, only : pi, lunation, Tm, Ta, g
  implicit none
  integer, intent(IN) :: nz
  real(8), intent(IN) :: time, dz
  real(8), intent(INOUT) :: T(0:nz)
  real(8), parameter :: delta = 0.05  ! thermal skin depth
  integer j
  real(8) z, phi

  phi = 2*pi*time/lunation
  do j=0,nz
     z = j*dz
     T(j) = Tm + Ta*exp(-z/delta)*sin(z/delta-phi) + g*z
  end do

end subroutine temperatureprofile
