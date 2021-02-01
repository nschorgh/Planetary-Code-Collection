pure subroutine soilthprop(porosity,fill,rhocobs,tiobs,layertype, &
     &     newrhoc,newti,icefrac)
!***********************************************************************
! soilthprop: assign thermal properties of icy soil or dirty ice
!
!     porositiy = void space / total volume
!     rhof = density of free ice in space not occupied by regolith [kg/m^3]
!     fill = rhof/icedensity <=1 (only relevant for layertype 1)
!     rhocobs = heat capacity per volume of dry regolith [J/m^3]
!     tiobs = thermal inertia of dry regolith [SI-units]
!     layertype: 1=interstitial ice, 2=pure ice or ice with dirt
!                3=pure ice, 4=ice-cemented soil, 5=custom values
!     icefrac: fraction of ice in icelayer (only relevant for layertype 2)
!     output are newti and newrhoc
!***********************************************************************
  implicit none
  integer, intent(IN) :: layertype
  real*8, intent(IN) :: porosity, fill, rhocobs, tiobs
  real*8, intent(OUT) :: newti, newrhoc
  real*8, intent(IN) :: icefrac
  real*8 kobs, cice, icedensity, kice
  !parameter (cice=2000.d0, icedensity=926.d0, kice=2.4d0) ! unaffected by scaling
  parameter (cice=1540.d0, icedensity=927.d0, kice=3.2d0) ! at 198 Kelvin
  real*8 fA, ki0, ki, k
  real*8, parameter :: kw=3. ! Mellon et al., JGR 102, 19357 (1997)

  kobs = tiobs**2/rhocobs
  ! k, rhoc, and ti are defined in between grid points
  ! rhof and T are defined on grid points

  newrhoc = -9999.
  newti  = -9999.

  select case (layertype)
  case (1) ! interstitial ice
     newrhoc = rhocobs + porosity*fill*icedensity*cice
     if (fill>0.) then
        !--linear addition (option A)
        k = porosity*fill*kice + kobs
        !--Mellon et al. 1997 (option B)
        ki0 = porosity/(1/kobs-(1-porosity)/kw)
        fA = sqrt(fill)
        ki = (1-fA)*ki0 + fA*kice
        !k = kw*ki/((1-porosity)*ki+porosity*kw)
     else
        k = kobs
     endif
     newti = sqrt(newrhoc*k)
     
  case (2)  ! massive ice (pure or dirty ice)
     newrhoc = rhocobs*(1.-icefrac)/(1.-porosity) + icefrac*icedensity*cice
     k = icefrac*kice + (1.-icefrac)*kw
     newti = sqrt(newrhoc*k)
  
  case (3)  ! all ice, special case of layertype 2, which doesn't use porosity
     newrhoc = icedensity*cice
     k = kice 
     newti = sqrt(newrhoc*k)
  
  case (4)  ! pores completely filled with ice, special case of layertype 1
     newrhoc = rhocobs + porosity*icedensity*cice
     k = porosity*kice + kobs ! option A, end-member case of type 1, option A 
     !k = kw*kice/((1-porosity)*kice+porosity*kw) ! option B, harmonic average
     newti = sqrt(newrhoc*k)

  case (5)  ! custom values
     ! values from Mellon et al. (2004) for ice-cemented soil
     newrhoc = 2018.*1040.
     k = 2.5
     newti = sqrt(newrhoc*k)

  case default
     error stop 'invalid layer type'
     
  end select
  
end subroutine soilthprop


      
subroutine smartgrid(nz,z,zdepth,thIn,rhoc,porosity,ti,rhocv,layertype,icefrac)
!***********************************************************************
! smartgrid: returns intelligently spaced grid and appropriate 
!            values of thermal inertia ti and rho*c in icy layer
!                  
!     INPUTS: 
!             nz = number of grid points
!             z = grid spacing for dry regolith 
!                 (will be partially overwritten)
!             zdepth = depth where ice table starts
!                      negative values indicate no ice
!             rhoc = heat capacity per volume of dry regolith [J/m^3]
!             thIn = thermal inertia of dry regolith [SI-units]
!             porosity = void space / total volume
!             layertypes are explained below  
!             icefrac = fraction of ice in icelayer
!
!     OUTPUTS: z = grid coordinates
!              vectors ti and rhocv
!***********************************************************************
  implicit none
  integer, intent(IN) :: nz, layertype
  real*8, intent(IN) :: zdepth, thIn, rhoc, porosity, icefrac
  real*8, intent(INOUT) :: z(nz)
  real*8, intent(OUT) :: ti(nz), rhocv(nz)
  integer j, b
  real*8 stretch, newrhoc, newti
  real*8, parameter :: NULL=0.
  
  if (zdepth>0 .and. zdepth<z(nz)) then

     select case (layertype)
     case (1)  ! interstitial ice
        call soilthprop(porosity,1.d0,rhoc,thIn,1,newrhoc,newti,NULL)
     case (2)  ! mostly ice (massive ice)
        call soilthprop(porosity,NULL,rhoc,thIn,2,newrhoc,newti,icefrac)
     case (3)  ! all ice (pure ice)
        call soilthprop(NULL,NULL,NULL,NULL,3,newrhoc,newti,NULL)
     case (4)  ! ice + rock + nothing else (ice-cemented soil)
        call soilthprop(porosity,NULL,rhoc,NULL,4,newrhoc,newti,NULL)
     case default
        error stop 'invalid layer type'
     end select

     ! Thermal skin depth is proportional to sqrt(kappa)
     ! thermal diffusivity kappa = k/(rho*c) = I^2/(rho*c)**2
     stretch = (newti/thIn)*(rhoc/newrhoc) ! ratio of sqrt(thermal diffusivity)
     
     b = 1
     do j=1,nz
        if (z(j)<=zdepth) then 
           b = j+1
           rhocv(j) = rhoc
           ti(j) = thIn
        else
           rhocv(j) = newrhoc
           ti(j) = newti
        endif
        ! print *,j,z(j),ti(j),rhocv(j)
     enddo
     do j=b+1,nz
        z(j) = z(b) + stretch*(z(j)-z(b))
     enddo
     
     ! print *,'zmax=',z(nz),' b=',b,' stretch=',stretch
     ! print *,'depth at b-1, b ',z(b-1),z(b)
     ! print *,'ti1=',thIn,' ti2=',newti
     ! print *,'rhoc1=',rhoc,' rhoc2=',newrhoc 
  endif
  
end subroutine smartgrid

