program cratersQ_equilbr
!***********************************************************************
!   cratersQ_equilbr: program to calculate surface energy balance
!             for constant sun position; equilibrium solution
!***********************************************************************
  use filemanager
  use allinterfaces
  use newhorizons
  implicit none
  real(8), parameter :: pi=3.1415926535897932, d2r=pi/180.
  real(8), parameter :: sigSB = 5.6704e-8  

  integer n, i, j, CCMAX
  real(8), dimension(NSx,NSy) :: h, SlopeAngle, azFac
  real(8), dimension(NSx,NSy) :: Qn, QIRin, Qrefl   ! incoming
  real(8) R, azSun, smax, emiss, betaSun, Qshadow1, Qshadow2
  real(8), dimension(NSx,NSy) :: Tsurf, viewsize, Qabs, albedo
  integer, dimension(NSx,NSy) :: cc ! for method 1
  integer(2), dimension(:,:,:), allocatable :: ii,jj ! for method 1
  real(4), dimension(:,:,:), allocatable :: VF  ! for methods 1 & 2
  integer, parameter :: Nflat = (NSx-2) * (NSy-2), T=200  ! for method 3
  real(8) w(Nflat)
  real(8), allocatable :: u(:,:), v(:,:)
  
  integer, parameter :: method = 1  ! choose input format 1...3
  
  ! azimuth in degrees east of north, 0=north facing
  albedo(:,:) = 0.12
  emiss = 0.95

  R=1.
  azSun = 180.*d2r
  betaSun = 10.*d2r

  call readdem(h)
  call difftopo(NSx,NSy,h,dx,dy,SlopeAngle,azFac)

  print *,'...reading horizons file...'
  call readhorizons

  print *,'...reading huge viewfactors file...'

  select case (method)
  case(1) ! fild with non-zero view factors
     CCMAX = getmaxfieldsize(NSx,NSy,vfn)
     print *,'... max field of view size=',CCMAX
     allocate(ii(NSx,NSy,CCMAX), jj(NSx,NSy,CCMAX), VF(NSx,NSy,CCMAX))
     !call getfieldofview(NSx,NSy,ffn,cc,ii,jj,dO12,landsize,CCMAX)
     call getviewfactors(NSx,NSy,vfn,cc,ii,jj,VF,viewsize,CCMAX)
  case(2) ! full square view factor matrix
     CCMAX = (NSx-2) * (NSy-2)
     allocate(VF(NSx,NSy,CCMAX))
     call getviewfactors_full(NSx,NSy,vfn,viewsize,VF)
  case(3) ! truncated or full SVD (outputs of xsvdcmp)
     print *,'...reading SVD files...'
     allocate(u(Nflat,T), v(Nflat,T))
     call readtruncSVD(U,w,V,Nflat,T,'./','.dat')
  case default
     stop 'no valid method'
  end select
     
  print *,'...calculating...'  
  print *,'beta=',betaSun/d2r,'azSun=',azSun/d2r

  do i=2,NSx-1
     do j=2,NSy-1
        smax = getonehorizon(i,j,azSun)
        !smax = 0.
        Qn(i,j) = flux_wshad(R,sin(betaSun),azSun,SlopeAngle(i,j),azFac(i,j),smax)
     enddo
  enddo

  Tsurf(:,:) = ((1-albedo)*Qn/sigSB)**0.25
  Qrefl=0.; QIRin=0.

  do n=1,7
     !print *,'Iteration',n

     select case(method)
     case(1)
        call update_terrain_irradiance(VF,cc,ii,jj,Qn,albedo,emiss,Qrefl,QIRin,Tsurf)
     case(2)
        call update_terrain_irradiance_full(VF,Qn,albedo,emiss,Qrefl,QIRin,Tsurf)
     case(3)
        call update_terrain_irradiance_SVD(T,U,w(1:T),V,Qn,albedo,emiss,Qrefl,QIRin,Tsurf)
     end select
     
     Qshadow1 = 0.; Qshadow2 = 0.
     do i=2,NSx-1
        Qshadow1 = Qshadow1 + sum(Qrefl(i,:),1,Qn(i,:)==0)
        Qshadow2 = Qshadow2 + sum(QIRin(i,:),1,Qn(i,:)==0)
     enddo
     
     ! Qabs is Q absorbed
     Qabs(:,:) = (1.-albedo(:,:))*(Qn+Qrefl) + emiss*QIRin
     
     Tsurf(:,:) = (Qabs/sigSB/emiss)**0.25

     print *,'Iteration',n,Qshadow1,Qshadow2
  enddo

  open(unit=21,file='qinst.dat',action='write')
  do i=2,NSx-1
     do j=2,NSy-1
        write(21,'(2(i4,1x),f9.2,2x,f6.4,4(1x,f6.1),1x,f5.1)') &
             & i,j,h(i,j),SlopeAngle(i,j), &
             & Qn(i,j),Qirin(i,j),Qrefl(i,j),Qabs(i,j),Tsurf(i,j)
     enddo
  enddo
  close(21)

end program cratersQ_equilbr


