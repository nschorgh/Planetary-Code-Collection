program mars_mapii
!***********************************************************************
! mars_mapii: program to calculate depth to equilibrium ice table for
!             a list of input parameters (can be the entire globe)
!***********************************************************************
  implicit none
  real*8, parameter :: pi=3.1415926535897932, d2r=pi/180., marsDay=88775.244

  integer nz, k
  real*8 dt, zmax, zfac, zdepth, icefrac, zacc
  real*8 latitude, thIn, albedo0, fracIR, fracDust, delta
  real*8 Fgeotherm, rhoc, lon, Tfrost, pfrost
  real*8 avdrho, junk, Tb, patm !, htopo
  real*8 zequil, zdepth_old, stretch
  real*8, external :: psv, rtbis, frostpoint, stretchfactor
      
  !-set input parameters
  dt=0.02
  nz=80; zfac=1.05;
  fracIR=0.04; fracDust=0.01
  Fgeotherm=0.0
  icefrac = 0.4  ! porosity for layertype 1
  
  zacc = 0.1  ! desired min. relative accuracy of ice table depth
  patm = 600.
  
  print *,'RUNNING MARS_MAP - EQUILIBRIUM ICE TABLE VIA ITERATION'
  write(*,*) 'Global model parameters'
  write(*,*) 'nz=',nz,' zfac=',zfac,' dt=',dt
  write(*,*) 'fracIR=',fracIR,' fracDust=',fracDust
  write(*,*) 'ice fraction=',icefrac
  write(*,*) 'Fgeotherm=',Fgeotherm
  write(*,*) 'Minimum ice depth accuracy dz/z=',zacc
  write(*,*) 'Atmospheric pressure=',patm
  
  open(unit=20,file='mapgrid.dat',status='old',action='read') ! the only input
  
  open(unit=30,file='z',action='write');
  open(unit=33,file='mapgrid2.dat',action='write')
  open(unit=34,file='means',action='write')
  
  do ! loop through sites
     read(20,*,end=80) lon,latitude,albedo0,thIn,Tfrost,zdepth
     !read(20,*,end=80) lon,latitude,albedo0,thIn,pfrost,patm
     
     !patm=520.*exp(-htopo/10800.)
     fracIR=0.04*patm/600.; fracDust=0.02*patm/600.  ! use if patm is available
     
     !Tfrost = frostpoint(pfrost)
     pfrost = psv(Tfrost)
     
     print *,lon,latitude,albedo0,thIn,Tfrost
     if (albedo0==-9999..or.thIn==-9999..or.Tfrost==-9999.) cycle
     ! Empirical relation from Mellon & Jakosky:
     rhoc = 800.*(150.+100.*sqrt(34.2+0.714*thIn))
     delta = thIn/rhoc*sqrt(marsDay/pi) ! diurnal skin depth
     zmax = 5.*26.*delta 
     Tb = -1.e32

     ! iterations of equilibrium ice table due to tifeedback
     if (zdepth<0.) zdepth = zmax
     zdepth_old = zdepth
     call jsub(zdepth, latitude*d2r, albedo0, thIn, pfrost, &
          &        nz/2, rhoc, fracIR, fracDust, patm, Fgeotherm, 2*dt, zfac, &
          &        icefrac, 1, Tb, junk, zequil)
     do k=1,12
        call jsub(zdepth, latitude*d2r, albedo0, thIn, pfrost, &
             &     nz,   rhoc, fracIR, fracDust, patm, Fgeotherm,   dt, zfac, &
             &     icefrac, 0, Tb, avdrho, zequil)

        zdepth = zequil
        if (zequil<0.) exit
        if (zdepth > zdepth_old) then ! if advancing downward
           ! move downward at reduced rate
           stretch = stretchfactor(thIn,rhoc,icefrac,1,0.d0) ! match jsub
           zdepth = zdepth_old + (zdepth-zdepth_old)/stretch
        end if
        !print *,k,'zdepth=',zdepth
        !if (zdepth > zdepth_old) exit ! appropriate if starting at zdepth=zmax
        if (abs(zdepth-zdepth_old) < zacc/2.*zdepth) exit
        zdepth_old = zdepth
     end do

     write(33,'(f7.2,1x,f7.3,2x,f5.3,2x,f5.1,1x,f7.1,2x,f10.4)') &
          &        lon,latitude,albedo0,thIn,Tfrost,zdepth
  enddo

80 continue
  close(20)
  close(30)
  close(33)
end program mars_mapii



function stretchfactor(thIn,rhoc,porosity,layertype,icefrac)
!***********************************************************************
! stretchfactor: returns ratio of sqrt(thermal diffusivities), a  
!                quantity relevant for the change of the temperature
!                gradient across the ice table
!    
!   INPUTS: 
!           rhoc = heat capacity per volume of ice-free regolith [J/m^3]
!           thIn = thermal inertia of ice-free regolith [SI-units]
!           porosity = void space / total volume
!           layertypes are explained below  
!           icefrac = fraction of ice in icelayer
!***********************************************************************
  implicit none
  real*8 stretchfactor
  integer, intent(IN) :: layertype
  real*8, intent(IN) :: thIn, rhoc, porosity, icefrac
  real*8 newrhoc, newti
  real*8, parameter :: NULL=0.
  
  select case (layertype)
  case (1)  ! interstitial ice
     call soilthprop(porosity,1.d0,rhoc,thIn,1,newrhoc,newti,NULL)
  case (2)  ! mostly ice (massive ice)
     call soilthprop(porosity,NULL,rhoc,thIn,2,newrhoc,newti,icefrac)
  case (3)  ! all ice (pure ice)
     call soilthprop(NULL,NULL,NULL,NULL,3,newrhoc,newti,NULL)
  case (4)  ! ice + rock + nothing else
     call soilthprop(porosity,NULL,rhoc,NULL,4,newrhoc,newti,NULL)
  case default
     error stop 'invalid layer type'
  end select

  ! thermal conductivity k = I^2/(rho*c)
  ! Thermal skin depth is proportional to sqrt(thermal diffusivity).
  ! thermal diffusivity kappa = k/(rho*c) = I^2/(rho*c)**2
  !stretchfactor = (newti/thIn)*(rhoc/newrhoc) ! sqrt(icy)/sqrt(ice-free)
  stretchfactor = (newti/thIn)**2 / (newrhoc/rhoc)  ! k(icy)/k(ice-free)
  
  !print *,'stretch=',stretchfactor
end function stretchfactor


