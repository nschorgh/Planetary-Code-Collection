
module allinterfaces
  ! interfaces from subroutines and functions

  ! begin topos.f90
  interface
     subroutine readdem(h)
       use filemanager
       implicit none
       real(8), intent(OUT) :: h(NSx,NSy)
     end subroutine readdem
  end interface

  ! begin shadows_subs.f90
  interface
     pure subroutine findonehorizon(h,i0,j0,azRay,smax)
       use filemanager, only : NSx,NSy,RMAX
       integer, intent(IN) :: i0,j0
       real(8), intent(IN) :: h(NSx,NSy),azRay
       real(8), intent(OUT) :: smax
     end subroutine findonehorizon
  end interface
  
  interface
     subroutine findonehorizon_wsort(h,i0,j0,azRay,smax,visibility)
       ! finds horizon and determines visibility for one azimuth ray
       use filemanager
       integer, intent(IN) :: i0,j0
       real(8), intent(IN) :: h(NSx,NSy),azRay
       real(8), intent(OUT) :: smax
       logical, intent(INOUT) :: visibility(NSx,NSy)
     end subroutine findonehorizon_wsort
  end interface

  interface
     subroutine findonehorizon_allaz(h,i0,j0,naz,smax)
       use filemanager, only : NSx,NSy,RMAX
       implicit none
       integer, intent(IN) :: i0,j0,naz
       real(8), intent(IN) :: h(NSx,NSy)
       real(8), intent(OUT) :: smax(naz)
     end subroutine findonehorizon_allaz
  end interface

  interface
     elemental function diffangle(a1,a2)
       real(8) diffangle
       real(8), intent(IN) :: a1,a2
     end function diffangle
  end interface


  ! begin fieldofview_subs.f90
  interface
     subroutine find3dangle(h,i0,j0,unit,visibility)
       use filemanager
       implicit none
       integer, intent(IN) :: i0,j0,unit
       real(8), intent(IN) :: h(NSx,NSy)
       logical, intent(IN) :: visibility(NSx,NSy)
     end subroutine find3dangle
  end interface

  interface
     elemental subroutine xyz2thetaphi(x,y,z,theta,phi)
       implicit none
       real(8), intent(IN) :: x,y,z
       real(8), intent(OUT) :: theta,phi
     end subroutine xyz2thetaphi
  end interface

  interface
     pure function area_spherical_quadrangle(phi,theta)
       implicit none
       real(8), intent(IN) :: phi(4), theta(4)
       real(8) area_spherical_quadrangle
     end function area_spherical_quadrangle
  end interface

  interface
     elemental function distanceonsphere(phi1,theta1,phi2,theta2)
       implicit none
       real(8), intent(IN) :: phi1,phi2,theta1,theta2
       real(8) distanceonsphere
     end function distanceonsphere
  end interface

  interface
     subroutine refinevisibility(i0,j0,h,visibility)
       use filemanager, only : NSx,NSy
       implicit none
       integer, intent(IN) :: i0,j0
       logical, intent(INOUT) :: visibility(NSx,NSy)
       real(8), intent(IN) :: h(NSx,NSy)
     end subroutine refinevisibility
  end interface

  
  ! begin crater_common.f90
  interface
     elemental function horizontaldistance(i1,j1,i2,j2)
       implicit none
       real(8) horizontaldistance
       integer, intent(IN) :: i1,j1,i2,j2
     end function horizontaldistance
  end interface

  interface
     pure function azimuth(i1,j1,i2,j2)
       implicit none
       real(8) azimuth
       integer, intent(IN) :: i1,j1,i2,j2
     end function azimuth
  end interface

  interface
     real(8) function viewing_angle(i0,j0,i,j,h)
       use filemanager, only : NSx,NSy,dx,dy
       implicit none
       integer, intent(IN) :: i0,j0,i,j
       real(8), intent(IN) :: h(NSx,NSy)
     end function viewing_angle
  end interface

  interface
     pure subroutine difftopo1(i,j,h,surfaceSlope,az)
       use filemanager, only : NSx,NSy,dx,dy
       implicit none
       integer, intent(IN) :: i,j
       real(8), intent(IN) :: h(NSx,NSy)
       real(8), intent(OUT) :: surfaceSlope,az
     end subroutine difftopo1
  end interface

  ! begin model_subs.f90
  interface
     subroutine gethorizon(i0,j0,azSun,smax,first)
       use filemanager
       implicit none
       integer, intent(IN) :: i0,j0
       real(8), intent(IN) :: azSun
       real(8), intent(OUT) :: smax
       logical, intent(IN) :: first
     end subroutine gethorizon
  end interface

  interface
     pure subroutine difftopo(NSx,NSy,h,dx,dy,surfaceSlope,azFac)
       implicit none
       integer, intent(IN) :: NSx,NSy
       real(8), intent(IN) :: h(NSx,NSy),dx,dy
       real(8), intent(OUT), dimension(NSx,NSy) :: surfaceSlope,azFac
     end subroutine difftopo
  end interface

  interface
     elemental subroutine equatorial2horizontal(decl,latitude,HA,sinbeta,azimuth)
       real(8), intent(IN) :: decl,latitude,HA
       real(8), intent(OUT) :: sinbeta,azimuth
     end subroutine equatorial2horizontal
  end interface

  interface
     elemental real(8) function flux_wshad(R,sinbeta,azSun,surfaceSlope,azFac,smax)
       real(8), intent(IN) :: R,azSun,sinbeta,surfaceSlope,azFac,smax
     end function flux_wshad
  end interface

  interface
     subroutine getfieldofview(NSx,NSy,fileext,cc,ia,ja,dOh,skysize,CCMAX)
       integer, intent(IN) :: NSx, NSy
       character(len=*), intent(IN) :: fileext
       integer, intent(IN) :: CCMAX
       integer, intent(OUT) :: cc(NSx,NSy) ! number of cells in field of view
       integer(2), intent(OUT), dimension(NSx,NSy,CCMAX) :: ia, ja
       real(4), intent(OUT), dimension(NSx,NSy,CCMAX) :: dOh
       real(8), intent(OUT) :: skysize(NSx,NSy)
     end subroutine getfieldofview
  end interface

  interface
     subroutine getmaxfieldsize(NSx,NSy,fileext,maxsize,type)
       implicit none
       integer, intent(IN) :: NSx,NSy,type
       character(len=*), intent(IN) :: fileext
       integer, intent(OUT) :: maxsize
     end subroutine getmaxfieldsize
  end interface

  interface
     subroutine compactoutput(unit,value,nr)
       implicit none
       integer, intent(IN) :: unit,nr
       real(8), intent(IN) :: value(nr)
     end subroutine compactoutput
  end interface

  ! begin mk_atmosphere.f90
  interface
     real(8) function mk_atmosphere(Z,I0,D0)
       implicit none
       real(8), intent(IN) :: Z  ! solar zenith angle (radians)
       real(8), intent(OUT) :: I0  ! clear-sky direct irradiance (unitless fraction) = transmittance
       real(8), intent(OUT) :: D0  ! clear-sky diffuse irradiance (unitless fraction)
     end function mk_atmosphere
  end interface
  
  ! begin cratersQ_mars_p.f90
  interface
     elemental function flux_mars(R,decl,latitude,HA,albedo, &
          &   fracir,fracdust,surfaceSlope,azFac,smax)
       implicit none
       real(8) flux_mars
       real(8), intent(IN) :: R,decl,latitude,HA,albedo,fracIR,fracDust
       real(8), intent(IN) :: surfaceSlope,azFac,smax
     end function flux_mars
  end interface

  ! Fortran 77 programs
  interface
     subroutine hpsort(n,ra,ind)
       implicit none
       INTEGER, intent(IN) :: n
       REAL(8), intent(INOUT) :: ra(n)
       INTEGER, intent(OUT) :: ind(n)
     end subroutine hpsort
  end interface

  interface
     subroutine conductionQ(nz,z,dt,Qn,Qnp1,T,ti,rhoc,emiss,Tsurf,Fgeotherm,Fsurf)
       implicit none
       integer NMAX
       parameter (NMAX=1000)
       integer, intent(IN) :: nz
       real*8, intent(IN) :: z(NMAX), dt, Qn, Qnp1, ti(NMAX),rhoc(NMAX), emiss, Fgeotherm
       real*8, intent(INOUT) :: T(NMAX), Tsurf
       real*8, intent(OUT) :: Fsurf
     end subroutine conductionQ
  end interface

  interface
     subroutine conductionT(nz,z,dt,T,Tsurf,Tsurfp1,ti,rhoc,Fgeotherm,Fsurf)
      implicit none
      integer NMAX
      parameter (NMAX=1000)
      integer, intent(IN) :: nz
      real*8, intent(IN) :: z(NMAX), dt, T(NMAX), Tsurf, Tsurfp1, ti(NMAX), rhoc(NMAX)
      real*8, intent(IN) :: Fgeotherm
      real*8, intent(OUT) :: Fsurf
    end subroutine conductionT
 end interface

end module allinterfaces
