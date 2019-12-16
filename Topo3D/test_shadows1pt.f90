program test_shadows1pt
  ! calculate horizons for one point
  use filemanager, only : NSx, NSy, fileext, dx, dy, RMAX
  use allinterfaces
  use azRays, naz_h => naz
  implicit none
  real(8), parameter :: r2d = 180./pi
  integer i, j, k
  real(8) h(NSx,NSy)

  ! azimuth in degrees east of north, 0=north facing, 0...2*pi
  
  call readdem(h)
  print *,'...finished reading topography... ',fileext
  print *,'# domain size =',NSx*dx,'x',NSy*dy
  print *,'# azimuth rays = ',naz_h
  print *,'# fully sampled radius =',min(dx,dy)*naz_h/(2*pi)

  i=14; j=29
  print *,i,j
  
  block ! plain
    real(8) smax(naz_h)
    call findallhorizon1(h,i,j,naz_h,smax(:))
  
    do k=1,naz_h
       write(31,'(f5.1,2x,f7.3)') azRay(k)*r2d,atan(smax(k))*r2d
    end do

  end block


  block ! with multigrid
    use newmultigrid
    integer LMAX, LACT
    real(8) RMG, smax2(naz_h)
    
    print *,'multigrid parameters'
    RMG = naz_h*max(dx,dy)/(2*pi)
    LMAX = floor( log(sqrt( (NSx*dx)**2+(NSy*dy)**2 ) /RMG) / log(2.) )
    print *,'# log2(domain size/RMG) =',LMAX
    LMAX = min(10,LMAX)
    call downsample_all(h,LMAX,LACT)
    LMAX = min(LACT,LMAX)
    print *,'# levels allocated = ',LMAX
    
    if (LMAX>1) then
       call findallhorizon_MGR(h,i,j,naz_h,smax2,RMG,LMAX)
       
       do k=1,naz_h
          write(32,'(f5.1,2x,f7.3)') azRay(k)*r2d,atan(smax2(k))*r2d
       end do
    end if
  end block


  block ! with sort and visibility
    use findvisibletopo, naz_f => naz
    real(8) smax3(naz_f), azSun
    logical :: visibility(NSx,NSy) = .false.

    call findallhorizon_wsort_v3(h,i,j,smax3,visibility)
    
    do k=1,naz_f
       azSun = 360./real(naz_f)*(k-1)
       write(33,'(f5.1,2x,f7.3)') azSun,atan(smax3(k))*r2d
    end do
    
  end block
  
end program test_shadows1pt
