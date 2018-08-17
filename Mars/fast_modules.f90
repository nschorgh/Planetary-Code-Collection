module miscparameters
  ! miscellaneous parameters that are very constant
  real(8), parameter :: pi=3.1415926535897932, d2r=pi/180.
  integer, parameter :: NMAX=1000
  real(8), parameter :: marsDay=88775.244, solsperyear=668.60
  real(8), parameter :: icedensity=927.
  real(8), parameter :: earthDay=86400.
  real(8), parameter :: sigSB=5.6704e-8
  ! for reference here are some parameters that are literally coded
  !   mass of H2O molecule = 18
  !   universal gas constant = 8314 
  !   length of Earth year in days = 365.24
end module miscparameters


module allinterfaces
  ! interfaces from Fortran 90 subroutines and functions
  ! comments have been stripped

  !begin fast_subs_mars.f90

  interface
     subroutine icelayer_mars(bigstep,nz,NP,thIn,rhoc,z,porosity,pfrost, &
          & Tb,zdepthF,zdepthE,porefill,Tmean1,Tmean3,zdepthG, &
          & latitude,albedo,p0,ecc,omega,eps,icefrac,zdepthT,avrho1, &
          & avrho1prescribed)
       use miscparameters, only : d2r, NMAX, icedensity
       implicit none
       integer, intent(IN) :: nz, NP
       real(8), intent(IN) :: bigstep
       real(8), intent(IN), dimension(NP) :: thIn, rhoc
       real(8), intent(IN) :: z(NMAX)
       real(8), intent(IN) :: porosity, pfrost(NP)
       real(8), intent(INOUT) :: Tb(NP), porefill(nz,NP), zdepthF(NP), zdepthT(NP)
       real(8), intent(OUT), dimension(NP) :: zdepthE, Tmean1, Tmean3, zdepthG
       real(8), intent(IN), dimension(NP) :: latitude, albedo, p0
       real(8), intent(IN) :: ecc, omega, eps, icefrac
       real(8), intent(OUT) :: avrho1(NP)
       real(8), intent(IN), optional :: avrho1prescribed(NP)
     end subroutine icelayer_mars
  end interface

  interface 
     subroutine ajsub_mars(typeT, latitude, albedo0, pfrost, nz, z, ti, rhocv, &
          &  fracIR, fracDust, p0, ecc, omega, eps, avdrho, avdrhoP, avrho1, &
          &  Tb, zdepthE, typeF, zdepthF, ypp, porefill, Tmean1, Tmean3, &
          &  B, typeG, zdepthG, avrho1prescribed)
       use miscparameters
       implicit none
       integer, intent(IN) :: nz, typeT
       real(8), intent(IN) :: latitude 
       real(8), intent(IN) :: albedo0, pfrost, z(NMAX)
       real(8), intent(IN) :: ti(NMAX), rhocv(NMAX), fracIR, fracDust, p0
       real(8), intent(IN) :: ecc, omega, eps, porefill(nz)
       real(8), intent(OUT) :: avdrho, avdrhoP 
       real(8), intent(OUT) :: avrho1
       real(8), intent(INOUT) :: Tb, Tmean1
       integer, intent(OUT) :: typeF 
       real(8), intent(OUT) :: zdepthE, zdepthF 
       real(8), intent(OUT) :: ypp(nz) 
       real(8), intent(OUT) :: Tmean3, zdepthG
       real(8), intent(IN) :: B
       integer, intent(OUT) :: typeG
       real(8), intent(IN), optional :: avrho1prescribed
     end subroutine ajsub_mars
  end interface
  
  !end of fast_subs_mars.f90
  !begin fast_subs_univ.f90

  interface
     pure function zint(y1,y2,z1,z2)
       implicit none
       real(8), intent(IN) :: y1,y2,z1,z2
       real(8) zint
     end function zint
  end interface

  interface
     pure function equildepth(nz, z, rhosatav, rhosatav0, avrho1)
       implicit none
       integer, intent(IN) :: nz
       real(8), intent(IN) :: z(nz), rhosatav(nz)
       real(8), intent(IN) :: rhosatav0, avrho1
       real(8) zint, equildepth
       external zint
     end function equildepth
  end interface

  interface
     subroutine depths_avmeth(typeT, nz, z, rhosatav, rhosatav0, rlow, avrho1,  &
          & porefill, typeF, zdepthF, B, ypp, typeG, zdepthG)
       use miscparameters, only : icedensity
       implicit none
       integer, intent(IN) :: nz, typeT
       real(8), intent(IN), dimension(nz) :: z, rhosatav, porefill
       real(8), intent(IN) :: rhosatav0, rlow, avrho1
       integer, intent(INOUT) :: typeF
       real(8), intent(INOUT) :: zdepthF
       real(8), intent(IN) :: B 
       real(8), intent(OUT) :: ypp(nz), zdepthG
       integer, intent(INOUT) :: typeG
       real(8), external :: zint, deriv1_onesided, colint
     end subroutine depths_avmeth
  end interface

  interface
     pure function constriction(porefill)
       implicit none
       real(8), intent(IN) :: porefill
       real(8) constriction
     end function constriction
  end interface
  
  interface
     subroutine assignthermalproperties(nz,thIn,rhoc,ti,rhocv, &
          &                      typeT,icefrac,porosity,porefill)
       implicit none
       integer, intent(IN) :: nz
       integer, intent(IN), optional :: typeT
       real(8), intent(IN), optional :: icefrac
       real(8), intent(IN) :: thIn, rhoc
       real(8), intent(IN), optional :: porosity, porefill(nz)
       real(8), intent(OUT) :: ti(nz), rhocv(nz)
     end subroutine assignthermalproperties
  end interface

  interface
     pure subroutine icechanges_poreonly(nz,z,typeF,typeG,avdrhoP,ypp,B,porefill)
       implicit none
       integer, intent(IN) :: nz, typeF, typeG
       real(8), intent(IN) :: z(nz), ypp(nz), avdrhoP, B
       real(8), intent(INOUT) :: porefill(nz)
     end subroutine icechanges_poreonly
  end interface

  interface
     pure subroutine icechanges(nz,z,typeF,avdrho,avdrhoP,ypp, &
          & Diff,porosity,icefrac,bigstep,zdepthT,porefill,typeG)
       implicit none
       integer, intent(IN) :: nz, typeF, typeG
       real(8), intent(IN) :: z(nz), ypp(nz), avdrho, avdrhoP
       real(8), intent(IN) :: Diff, porosity, icefrac, bigstep
       real(8), intent(INOUT) :: zdepthT, porefill(nz)
     end subroutine icechanges
  end interface

  interface
     subroutine compactoutput(unit,porefill,nz)
       implicit none
       integer, intent(IN) :: unit,nz
       real(8), intent(IN) :: porefill(nz)
     end subroutine compactoutput
  end interface

  !end of fast_subs_univ
  !begin derivs.f90 

  interface
     subroutine deriv1(z,nz,y,y0,yNp1,yp)
       implicit none
       integer, intent(IN) :: nz
       real*8, intent(IN) :: z(nz),y(nz),y0,yNp1
       real*8, intent(OUT) :: yp(nz)
     end subroutine deriv1
  end interface

  interface
     subroutine deriv2_full(z,nz,a,b,a0,b0,bNp1,yp2)
       implicit none
       integer, intent(IN) :: nz
       real*8, intent(IN) :: z(nz),a(nz),b(nz),a0,b0,bNp1
       real*8, intent(OUT) :: yp2(nz)
     end subroutine deriv2_full
  end interface

  interface
     subroutine deriv2(z,nz,a,b,a0,b0,bNp1,yp2,c)
       implicit none
       integer, intent(IN) :: nz
       real*8, intent(IN) :: z(nz),a(nz),b(nz),a0,b0,bNp1,c(3,nz)
       real*8, intent(OUT) :: yp2(nz)
     end subroutine deriv2
  end interface

  interface
     subroutine deriv2_simple(z,nz,y,y0,yNp1,yp2)
       implicit none
       integer, intent(IN) :: nz
       real*8, intent(IN) :: z(nz),y(nz),y0,yNp1
       real*8, intent(OUT) :: yp2(nz)
     end subroutine deriv2_simple
  end interface

  interface
     real(8) pure function deriv1_onesided(j,z,nz,y)
       implicit none
       integer, intent(IN) :: nz,j
       real(8), intent(IN) :: z(nz),y(nz)
     end function deriv1_onesided
  end interface

  !end of derivs.f90
  !begin fast_subs_exper.f90

  interface
     subroutine icelayer_exper(bigstep, nz, thIn, rhoc, z, porosity, pfrost, &
          & zdepthT, icefrac, zdepthF, zdepthE, porefill, Tmean, Tampl, Diff, zdepthG)
       use miscparameters, only : d2r, NMAX, icedensity
       implicit none
       integer, intent(IN) :: nz
       real(8), intent(IN) :: bigstep
       real(8), intent(IN) :: thIn, rhoc, z(NMAX), porosity, pfrost
       real(8), intent(INOUT) :: zdepthT, zdepthF, porefill(nz)
       real(8), intent(OUT) :: zdepthE
       real(8), intent(IN) :: icefrac, Diff, Tmean, Tampl
       real(8), intent(OUT) :: zdepthG
     end subroutine icelayer_exper
  end interface

  interface
     subroutine ajsub_exper(typeT, nz, z, ti, rhocv, pfrost, Tmean, Tampl, &
          &     avdrho, avdrhoP, zdepthE, typeF, zdepthF, ypp, porefill, & 
          &     B, typeG, zdepthG)
       use miscparameters, only : NMAX, solsperyear, marsDay
       implicit none
       integer, intent(IN) :: nz, typeT
       real(8), intent(IN) :: z(NMAX), ti(NMAX), rhocv(NMAX), pfrost
       real(8), intent(IN) :: Tmean, Tampl
       real(8), intent(OUT) :: avdrho, avdrhoP 
       real(8), intent(OUT) :: zdepthE
       integer, intent(OUT) :: typeF 
       real(8), intent(INOUT) :: zdepthF 
       real(8), intent(OUT) :: ypp(nz)
       real(8), intent(IN) :: porefill(nz)
       real(8), intent(IN) :: B 
       integer, intent(OUT) :: typeG
       real(8), intent(OUT) :: zdepthG
       real(8), external :: Tsurface, psv
       real(8), external :: equildepth 
     end subroutine ajsub_exper
  end interface
  
  !end fast_subs_exper

  ! Common/*.f
  interface
     pure function colint(y,z,nz,i1,i2)
       implicit none
       integer, intent(IN) :: nz, i1, i2
       real(8), intent(IN) :: y(nz),z(nz)
       real(8) colint
     end function colint
  end interface

end module allinterfaces
