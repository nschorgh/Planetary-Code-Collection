MODULE allinterfaces
  ! interfaces from subroutines and functions

  ! begin topo3d_geometry.f90
  interface
     elemental function diffangle(a1,a2)
       real(8) diffangle
       real(8), intent(IN) :: a1,a2
     end function diffangle
  end interface
  interface
     subroutine compactoutput(unit,value,nr)
       implicit none
       integer, intent(IN) :: unit,nr
       real(8), intent(IN) :: value(nr)
     end subroutine compactoutput
  end interface
  interface
     elemental function horizontaldistance1(x1,y1,x2,y2)
       implicit none
       real(8) horizontaldistance1
       real(8), intent(IN) :: x1,y1,x2,y2
     end function horizontaldistance1
  end interface
  interface
     elemental function azimuth1(x1,y1,x2,y2)
       implicit none
       real(8) azimuth1
       real(8), intent(IN) :: x1,y1,x2,y2
     end function azimuth1
  end interface
  interface
     elemental function horizontaldistance_square(x1,y1,x2,y2)
       implicit none
       real(8) horizontaldistance_square
       real(8), intent(IN) :: x1,y1,x2,y2
     end function horizontaldistance_square
  end interface
  interface
     subroutine downsample(NSx,NSy,h,hhalf)
       implicit none
       integer, intent(IN) :: NSx,NSy
       real(8), intent(IN) :: h(NSx,NSy)
       real(8), intent(OUT) :: hhalf(NSx/2,NSy/2)
     end subroutine downsample
  end interface
  interface
     pure subroutine difftopo1(NSx,NSy,i,j,h,dx,dy,surfaceSlope,az)
       implicit none
       integer, intent(IN) :: NSx, NSy, i, j
       real(8), intent(IN) :: h(NSx,NSy), dx, dy
       real(8), intent(OUT) :: surfaceSlope, az
     end subroutine difftopo1
  end interface
  
  ! begin shadows_subs.f90
  interface
     pure subroutine findallhorizon1(h,i0,j0,naz,smax)
       use filemanager, only : NSx,NSy,dx,dy
       implicit none
       integer, intent(IN) :: i0,j0,naz
       real(8), intent(IN) :: h(NSx,NSy)
       real(8), intent(OUT) :: smax(naz)
     end subroutine findallhorizon1
  end interface

  ! begin fieldofview_subs.f90
  interface
     function spherical_area(h,i,j,i0,j0)
       use filemanager, only : NSx, NSy, dx, dy
       implicit none
       real(8) spherical_area
       real(8), intent(IN) :: h(NSx,NSy)
       integer, intent(IN) :: i, j, i0, j0
     end function spherical_area
  end interface
  interface
     elemental subroutine xyz2thetaphi(x,y,z,theta,phi)
       implicit none
       real(8), intent(IN) :: x,y,z
       real(8), intent(OUT) :: theta,phi
     end subroutine xyz2thetaphi
  end interface
  interface
     pure subroutine h4tonormal(hip1,him1,hjp1,hjm1,nx,ny,nz)
       implicit none
       real(8), intent(IN) :: hip1, him1, hjp1, hjm1
       real(8), intent(OUT) :: nx, ny, nz
     end subroutine h4tonormal
  end interface
  interface
     elemental function area_cart(dx, dy, x20, y20, z20, n2x, n2y, n2z)
       implicit none
       real(8) area_cart
       real(8), intent(IN) :: dx, dy, x20, y20, z20
       real(8), intent(IN) :: n2x, n2y, n2z
     end function area_cart
  end interface
  interface
     elemental function viewfactor_cart(dx, dy, x20, y20, z20, &
          & n1x, n1y, n1z, n2x, n2y, n2z)
       implicit none
       real(8) viewfactor_cart
       real(8), intent(IN) :: dx, dy, x20, y20, z20
       real(8), intent(IN) :: n1x, n1y, n1z, n2x, n2y, n2z
     end function viewfactor_cart
  end interface

  ! begin topo3d_common.f90
  interface
     subroutine readdem(h)
       use filemanager
       implicit none
       real(8), intent(OUT) :: h(NSx,NSy)
     end subroutine readdem
  end interface
  interface
     pure function area_spherical_quadrangle(phi,theta)
       implicit none
       real(8), intent(IN) :: phi(4), theta(4)
       real(8) area_spherical_quadrangle
     end function area_spherical_quadrangle
  end interface
  interface
     pure function area_spherical_triangle(phi,theta)
       real(8), intent(IN) :: phi(3), theta(3)
       real(8) area_spherical_triangle
     end function area_spherical_triangle
  end interface
  interface
     subroutine slicer(NSx,ilower,iupper,extc)
       implicit none
       integer, intent(IN) :: NSx
       integer, intent(OUT) :: ilower, iupper
       character(5), intent(OUT) :: extc
     end subroutine slicer
  end interface
  
  ! begin topo3d_subs.f90
  interface
     pure subroutine difftopo(NSx,NSy,h,dx,dy,SlopeAngle,azFac)
       implicit none
       integer, intent(IN) :: NSx,NSy
       real(8), intent(IN) :: h(NSx,NSy),dx,dy
       real(8), intent(OUT), dimension(NSx,NSy) :: SlopeAngle,azFac
     end subroutine difftopo
  end interface
  interface
     pure subroutine difftopo2(h,SlopeAngle,azFac,Mx1,Mx2,My1,My2)
       use filemanager, only : NSx,NSy,dx,dy
       implicit none
       integer, intent(IN) :: Mx1,Mx2,My1,My2
       real(8), intent(IN) :: h(NSx,NSy)
       real(8), intent(OUT), dimension(Mx1:Mx2,My1:My2) :: SlopeAngle,azFac
     end subroutine difftopo2
  end interface
  interface
     elemental subroutine equatorial2horizontal(decl,latitude,HA, &
          & sinbeta,azimuth)
       real(8), intent(IN) :: decl,latitude,HA
       real(8), intent(OUT) :: sinbeta,azimuth
     end subroutine equatorial2horizontal
  end interface
  interface
     subroutine readviewfactors(NSx,NSy,vfn,CCMAX,cc,ia,ja,VF,viewsize)
       implicit none
       integer, intent(IN) :: NSx, NSy
       character(len=*), intent(IN) :: vfn
       integer, intent(IN) :: CCMAX
       integer, intent(OUT) :: cc(NSx,NSy)
       integer(2), intent(OUT), dimension(NSx,NSy,CCMAX) :: ia, ja
       real(4), intent(OUT), dimension(NSx,NSy,CCMAX) :: VF
       real(8), intent(OUT) :: viewsize(NSx,NSy)
     end subroutine readviewfactors
  end interface
  interface
     integer function getmaxfieldsize(NSx,NSy,ffn)
       implicit none
       integer, intent(IN) :: NSx,NSy
       character(len=*), intent(IN) :: ffn
     end function getmaxfieldsize
  end interface
  interface
     subroutine readtruncSVD(U,w,V,Nflat,T,pth,ext)
       implicit none
       integer, intent(IN) :: Nflat, T
       real(8), intent(OUT) :: U(Nflat,T), w(Nflat), V(Nflat,T)
       character(len=*), intent(IN) :: pth, ext
     end subroutine readtruncSVD
  end interface
  interface
     integer function countcolumns()
     end function countcolumns
  end interface
  interface
     subroutine update_terrain_irradiance(VF,cc,ii,jj,&
          & Qn,albedo,emiss,Qrefl,QIRin,Tsurf)
       use filemanager, only : NSx, NSy
       implicit none
       real(4), intent(IN) :: VF(:,:,:)
       integer, intent(IN) :: cc(NSx,NSy)
       integer(2), intent(IN) :: ii(:,:,:), jj(:,:,:)
       real(8), intent(IN), dimension(NSx,NSy) :: Qn, albedo, Tsurf
       real(8), intent(IN) :: emiss
       real(8), intent(INOUT), dimension(NSx,NSy) :: Qrefl, QIRin
     end subroutine update_terrain_irradiance
  end interface
  interface
     subroutine update_terrain_irradiance_full(VF, &
          & Qn,albedo,emiss,Qrefl,QIRin,Tsurf)
       use filemanager, only : NSx, NSy
       implicit none
       integer, parameter :: Nflat = (NSx-2)*(NSy-2)
       real(4), intent(IN) :: VF(NSx,NSy,Nflat)
       real(8), intent(IN), dimension(NSx,NSy) :: Qn, albedo, Tsurf
       real(8), intent(IN) :: emiss
       real(8), intent(INOUT), dimension(NSx,NSy) :: Qrefl, QIRin
     end subroutine update_terrain_irradiance_full
  end interface
  interface
     subroutine update_terrain_irradiance_SVD(T,U,w,V, &
          & Qn,albedo,emiss,Qrefl,QIRin,Tsurf)
       use filemanager, only : NSx, NSy
       implicit none
       integer, parameter :: Nflat = (NSx-2)*(NSy-2)
       integer, intent(IN) :: T
       real(8), intent(IN) :: U(:,:), w(T), V(:,:)
       real(8), intent(IN), dimension(NSx,NSy) :: Qn, albedo, Tsurf
       real(8), intent(IN) :: emiss
       real(8), intent(INOUT), dimension(NSx,NSy) :: Qrefl, QIRin
     end subroutine update_terrain_irradiance_SVD
  end interface
  
  ! mk_atmosphere.f90
  interface
     real(8) function mk_atmosphere(Z,I0,D0)
       implicit none
       real(8), intent(IN) :: Z
       real(8), intent(OUT) :: I0, D0
     end function mk_atmosphere
  end interface

  ! routines in Mars/
  interface
     pure subroutine flux_mars2(R,decl,latitude,HA,fracIR,fracDust, &
          &   SlopeAngle,azFac,emax,Qdir,Qscat,Qlw)
       implicit none
       real(8), intent(IN) :: R, decl, latitude, HA, SlopeAngle, azFac, emax
       real(8), intent(IN) :: fracIR, fracDust
       real(8), intent(OUT) :: Qdir, Qscat, Qlw
     end subroutine flux_mars2
  end interface
  interface
     subroutine marsorbit(dt0,tj,Ls,dec,r)
       implicit none
       real*8, intent(IN) :: dt0, tj
       real*8, intent(OUT) :: Ls, dec, r
     end subroutine marsorbit
  end interface
  interface
     subroutine marsclock24(JDUT,Deltat_J2000,Ls,dec,RM,Longitude_W,LTST)
       implicit none
       real*8, intent(IN) :: JDUT
       real*8, intent(OUT) :: Deltat_J2000
       real*8, intent(OUT) :: Ls
       real*8, intent(OUT) :: dec, RM
       real*8, intent(IN) :: Longitude_W
       real*8, intent(OUT) :: LTST
     end subroutine marsclock24
  end interface
  
  ! topod3d_subs_mars.f90
  interface
     subroutine subsurfaceconduction_mars(T,Tsurf,dtsec,Qn,Qnp1,m,Fsurf,init, &
          & Tco2frost,thIn,emiss)
       implicit none
       real(8), intent(INOUT) :: T(:), Tsurf, m, Fsurf
       real(8), intent(IN) :: dtsec, Qn, Qnp1
       logical, intent(IN) :: init
       real(8), intent(IN), optional :: Tco2frost, thIn, emiss
     end subroutine subsurfaceconduction_mars
  end interface
  interface
     elemental subroutine equilibrT_mars(Tsurf,dtsec,Qnm1,Qn,m,Tco2frost,emiss)
       implicit none
       real(8), intent(IN) :: Qn, Qnm1, Tco2frost, emiss, dtsec
       real(8), intent(OUT) :: Tsurf
       real(8), intent(INOUT) :: m
     end subroutine equilibrT_mars
  end interface
  interface
     pure function evap_ingersoll(T,p0)
       implicit none
       real(8) evap_ingersoll
       real(8), intent(IN) :: T,p0
     end function evap_ingersoll
  end interface

  ! f90 routines in Common/
  interface
     pure function flux_wshad(R,sinbeta,azSun,SlopeAngle,azFac,emax)
       real(8) flux_wshad
       real(8), intent(IN) :: R,azSun,sinbeta,SlopeAngle,azFac,emax
     end function flux_wshad
  end interface
  interface
     subroutine conductionQ(nz,z,dt,Qn,Qnp1,T,ti,rhoc,emiss,Tsurf,Fgeotherm, &
          & Fsurf)
       implicit none
       integer, intent(IN) :: nz
       real*8, intent(IN) :: z(nz), dt, Qn, Qnp1
       real*8, intent(INOUT) :: T(nz), Tsurf
       real*8, intent(IN) :: ti(nz), rhoc(nz), emiss, Fgeotherm
       real*8, intent(OUT) :: Fsurf
     end subroutine conductionQ
  end interface
  interface
     subroutine conductionT(nz,z,dt,T,Tsurf,Tsurfp1,ti,rhoc,Fgeotherm,Fsurf)
       implicit none
       integer, intent(IN) :: nz
       real*8, intent(IN) :: z(nz), dt, Tsurf, Tsurfp1, ti(nz), rhoc(nz)
       real*8, intent(IN) :: Fgeotherm
       real*8, intent(INOUT) :: T(nz)
       real*8, intent(OUT) :: Fsurf
     end subroutine conductionT
  end interface
  
  interface
     pure function psv(T)
       implicit none
       real(8), intent(IN) :: T
       real(8) psv
     end function psv
  end interface

END MODULE allinterfaces
