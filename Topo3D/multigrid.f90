! subroutines for multigrid acceleration of horizons calculation


MODULE multigrid
  use filemanager, only : NSx, NSy, dx, dy
  real(8), allocatable, private, target :: h2(:,:), h3(:,:), h4(:,:), h5(:,:)
  real(8), allocatable, private, target :: h6(:,:), h7(:,:), h8(:,:), h9(:,:)
  real(8), allocatable, private, target :: h10(:,:)

  interface
     pure subroutine horizon_core(x0,y0,h00,smax,i,j,h,P)
       use filemanager, only : NSx, NSy, dx, dy
       use azRays
       implicit none
       real(8), intent(IN) :: x0, y0, h00 
       real(8), intent(INOUT) :: smax(naz)
       integer, intent(IN) :: i, j, P 
       real(8), intent(IN) :: h(NSx/P,NSy/P)
     end subroutine horizon_core
  end interface

  
contains

  subroutine L2coarse(L,hcoarse)
    ! returns pointer to coarse topography array, levels 2...10
    implicit none
    integer, intent(IN) :: L
    real(8), intent(OUT), pointer :: hcoarse(:,:)
        
    select case (L)
    case (10)
       hcoarse => h10
    case (9)
       hcoarse => h9
    case (8)
       hcoarse => h8
    case (7)
       hcoarse => h7
    case (6)
       hcoarse => h6
    case (5)
       hcoarse => h5
    case (4)
       hcoarse => h4
    case (3)
       hcoarse => h3
    case (2)
       hcoarse => h2
    case default
       error stop 'L2coarse: level out of range'
    end select
  end subroutine L2coarse


  subroutine downsample_all(h,LMAX,LACT)
    use allinterfaces, only: downsample
    implicit none
    integer, intent(IN) :: LMAX
    real(8), intent(IN) :: h(NSx,NSy)
    integer, intent(OUT) :: LACT  ! levels actually allocated
    integer, parameter :: MINPWIDTH = 5
    integer NS

    NS = min(NSx,NSy)
    
    LACT = 1
    if (NS/2 > MINPWIDTH .and. LMAX>=2) then
       allocate(h2(NSx/2,NSy/2))
       call downsample(NSx,NSy,h,h2)
       LACT = 2
       if (NS/4 > MINPWIDTH .and. LMAX>=3) then
          allocate(h3(NSx/4,NSy/4))
          call downsample(NSx/2,NSy/2,h2,h3)
          LACT = 3
          if (NS/8 > MINPWIDTH .and. LMAX>=4) then
             allocate(h4(NSx/8,NSy/8))
             call downsample(NSx/4,NSy/4,h3,h4)
             LACT = 4
             if (NS/16 > MINPWIDTH .and. LMAX>=5) then
                allocate(h5(NSx/16,NSy/16));
                call downsample(NSx/8,NSy/8,h4,h5)
                LACT = 5
                if (NS/32 > MINPWIDTH .and. LMAX>=6) then
                   allocate(h6(NSx/32,NSy/32))
                   call downsample(NSx/16,NSy/16,h5,h6)
                   LACT = 6
                   if (NS/64 > MINPWIDTH .and. LMAX>=7) then
                      allocate(h7(NSx/64,NSy/64))
                      call downsample(NSx/32,NSy/32,h6,h7)
                      LACT = 7
                      if (NS/128 > MINPWIDTH .and. LMAX>=8) then
                         allocate(h8(NSx/128,NSy/128))
                         call downsample(NSx/64,NSy/64,h7,h8)
                         LACT = 8
                         if (NS/256 > MINPWIDTH .and. LMAX>=9) then
                            allocate(h9(NSx/256,NSy/256))
                            call downsample(NSx/128,NSy/128,h8,h9)
                            LACT = 9
                            if (NS/512 > MINPWIDTH .and. LMAX>=10) then
                               allocate(h10(NSx/512,NSy/512))
                               call downsample(NSx/256,NSy/256,h9,h10)
                               LACT = 10
                            endif
                         endif
                      endif
                   endif
                endif
             endif
          endif
       endif
    endif
  end subroutine downsample_all
  
  
  subroutine findallhorizon_MGR(h,i0,j0,naz,smax,RMG,L)
    ! based on shadow_subs.f90, but for multigrid method
    ! recursive implementation
    use allinterfaces, only: horizontaldistance1
    implicit none
    integer, intent(IN) :: i0, j0, naz, L
    real(8), intent(IN) :: h(NSx,NSy)
    real(8), intent(IN) :: RMG
    real(8), intent(OUT) :: smax(naz)
    integer P, ii, jj
    real(8) r, x0, y0, h00
    real(8), pointer :: hcoarse(:,:)

    smax(:) = 0.
    
    x0 = i0*dx; y0 = j0*dy; h00 = h(i0,j0)
    
    P = 2**(L-1)
    if (L>10 .or. L<2) error stop 'findallhorizon_MGR: invalid grid level'
    ! The top loop for the coarsest grid is different from all others
    do ii=1,NSx/P + 1; do jj=1,NSy/P+1  ! +1 so a last odd one gets included too
       r = horizontaldistance1(P*ii*dx,P*jj*dy,x0,y0)
       if (r>=P/2*RMG) then  ! do 1 coarse grid cell (root node)
          call L2coarse(L,hcoarse)
          call horizon_core(x0,y0,h00,smax,ii,jj,hcoarse,P)
          !print *,'Level',L,P*ii,P*jj,smax(90)
       else
          call findallhorizon_recursive(ii,jj,h,x0,y0,h00,naz,smax,RMG,L-1)
       endif
    enddo; enddo
  end subroutine findallhorizon_MGR


  recursive subroutine findallhorizon_recursive(i,j,h,x0,y0,h00,naz,smax,RMG,L)
    use allinterfaces, only: horizontaldistance1
    implicit none
    integer, intent(IN) :: naz, L, i, j
    real(8), intent(IN) :: h(NSx,NSy)
    real(8), intent(IN) :: x0, y0, h00, RMG
    real(8), intent(INOUT) :: smax(naz)
    integer ii, jj, P
    real(8) r
    real(8), pointer :: hcoarse(:,:)
    
    ! start with the coarsest grid
    ! if distance is too close, switch to finer grid
    
    P = 2**(L-1)
    if (L>=10) error stop &
         & 'findallhorizon_recursive: exceeds maximum number of levels'
    if (L<1) error stop 'findallhorizon_recursive: level is zero or negative'
    do ii=2*i-1,2*i; do jj=2*j-1,2*j ! 4 cells
       r = horizontaldistance1(P*ii*dx,P*jj*dy,x0,y0)
       if (P==1) r = 0.  ! should be redundant
       if (r>=P/2*RMG) then ! do 1 coarse grid cell
          if (L>1) then
             call L2coarse(L,hcoarse)
             call horizon_core(x0,y0,h00,smax,ii,jj,hcoarse,P)
          else ! L=1, reached finest grid
             if (ii<=1 .or. jj<=1 .or. ii>=NSx .or. jj>=NSy) cycle ! strip margin
             call horizon_core(x0,y0,h00,smax,ii,jj,h,P)
          end if
       else ! do 4 finer cells (descend on quad-tree)
          call findallhorizon_recursive(ii,jj,h,x0,y0,h00,naz,smax,RMG,L-1)
       endif
    enddo; enddo
  end subroutine findallhorizon_recursive

END MODULE multigrid

