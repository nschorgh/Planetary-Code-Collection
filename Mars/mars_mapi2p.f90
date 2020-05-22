program mars_mapi2p
!***********************************************************************
!   mars_mapi2p: program to calculate depth to ice table for a list of
!            input parameters (can be the entire globe, can also be
!            a list of slopes); incorporates terrain irradiance for
!            planar slopes from a flat floor at a different temperature  
!
!     Can be launched with mars_mapi2p_go.cmd
!     1/20/05: Added new feature to bracket root before loop  -oded
!     5/31/05: vectorized, converted to Fortran 90  -norbert
!     5/16/20: modernized  -norbert
!***********************************************************************

  implicit none
  integer, parameter :: NS=10  ! number of slopes for each site
  real*8, parameter :: pi=3.1415926535897932, d2r=pi/180., zero=0.
  real*8, parameter :: marsDay=88775.244
  
  integer nz, iargc, job_nr, j, line_nr, k
  real*8 dt, zmax, zfac, zdepth0, icefrac, zacc
  real*8 latitude, thIn, albedo0, fracIR, fracDust, delta
  real*8 Fgeotherm, rhoc, lon, Tfrost, pfrost, slpd, azFacd, patm
  real*8 slp(NS), azFac(NS), zdepth(NS), avdrho(NS), Tb(NS), zz(NS)
  real*8 psv, rtbis
  external psv, rtbis
  character(40) infile, outfile
  character(10) line_nr_string, job_nr_string

!-set global input parameters
  dt = 0.02
  nz=80; zfac=1.05
  fracIR=0.04; fracDust=0.02
  !Fgeotherm = 0.028
  Fgeotherm = 0.0
  icefrac = 0.4
  !icefrac = 0.0
  zacc = 0.1  ! desired min. relative accuracy of ice table depth
  
  if (iargc() /= 2) then
     print *,'USAGE: mars_mapi2p file.in job_nr'
     print *,"For example: 'a.out mapgrid.slp 2' uses 2nd line of file mapgrid.slp as input"
     stop
  endif

  call getarg( 1, infile)
  call getarg( 2, job_nr_string)
  read(job_nr_string,'(i4)') job_nr  ! string->integer
  
  write (*,*) 'infile:  ',infile
  write (*,*) 'nextrun: ',job_nr
  
  print *,'RUNNING MARS_MAP-ICE TABLE'
  write(*,*) 'Global model parameters'
  write(*,*) 'nz=',nz,' zfac=',zfac,' dt=',dt
  write(*,*) 'fracIR=',fracIR,' fracDust=',fracDust
  write(*,*) 'ice fraction=',icefrac
  write(*,*) 'Fgeotherm=',Fgeotherm
  write(*,*) 'Number of sites=',NS
  write(*,*) 'Minimum ice depth accuracy dz/z=',zacc

  open(unit=20,file=infile,status='old',action='read')  ! the only input
  do j=1,job_nr
     read(20,*,end=90) line_nr,lon,latitude,albedo0,thIn,Tfrost, &
          &        slpd,azFacd,zdepth0
  enddo
  close(20)
  print *,line_nr,lon,latitude,albedo0,thIn,Tfrost,slpd,azFacd,zdepth0

  write (line_nr_string, "(I0)") line_nr  ! integer->string
  outfile = 'mapgrid2.'//line_nr_string

  if (latitude>=0.) then    ! northern hemisphere
     patm = 700. ! Tco2frost = 148 K 
  else                      ! southern hemisphere
     patm = 400. ! Tco2frost = 145 K
  endif
  write(*,*) 'Atmospheric pressure=',patm
  
  azFac=azFacd*d2r
  if (NS==1)  slp(1) = 0.
  if (NS==2)  slp(1:2) = (/ zero, slpd /)*d2r
  if (NS>2) then
     slp=-9999
     ! first slope MUST be zero; length of the list must be NS
     slp = (/ 0., 10., 20., 30., 40., 50., 60., 70., 80., 90. /)
     if (minval(slp)<0.) stop 'inappropriate slopes or wrong NS'
     slp = slp*d2r
  endif
  if (slp(1)/=0.) stop 'first slope must be zero'
  
  if (albedo0==-9999. .or. thIn==-9999. .or. Tfrost==-9999.) then
     zdepth = zdepth0   ! otherwise undefined on output
     goto 80
  endif

  ! zdepth0 input is ignored
  ! Empirical relation from Mellon & Jakosky:
  rhoc = 800.*(150.+100.*sqrt(34.2+0.714*thIn))
  delta = thIn/rhoc*sqrt(marsDay/pi) ! diurnal skin depth
  zmax = 5.*26.*delta 
  pfrost = psv(Tfrost)
  Tb(1) = -1.e32

  print *,'ice depth on flat slope ...'
  call jsubv(1, zmax, latitude*d2r, albedo0, thIn, pfrost, &
       &     nz/2, rhoc, fracIR, fracDust, Fgeotherm, 2.*dt, zfac, &
       &     icefrac, zero, zero, 1, avdrho(1), Tb(1), patm)
  call jsubv(1, zmax, latitude*d2r, albedo0, thIn, pfrost, &
       &     nz,   rhoc, fracIR, fracDust, Fgeotherm,    dt, zfac, &
       &     icefrac, zero, zero, 0, avdrho(1), Tb(1), patm)
  print *, 'ice depth: ','  rho_ice-rho_surf'
  print *, zmax,'#',avdrho(1)
  if (avdrho(1)>=0.) then   ! no ice
     zdepth = -9999.
  else  
     zdepth0 = rtbis(delta/4.,zmax, zacc,avdrho(1), &
          &     latitude*d2r,albedo0,thIn,pfrost,nz,rhoc,fracIR, &
          &     fracDust,Fgeotherm,dt,zfac,icefrac,zero,zero,patm)
  endif
  print *,'Equilibrium ice table depth on flat slope = ',zdepth0

  print *,'ice depth on all other slopes ...'
  Tb(:) = -1.e32
  zz = spread(zmax,1,NS)
  if (zdepth0>0.) zz(1)=zdepth0  ! otherwise zz(1)=zmax
  call jsubv(NS, zz, latitude*d2r, albedo0, thIn, pfrost, &
       &     nz/2, rhoc, fracIR, fracDust, Fgeotherm, 2.*dt, zfac, &
       &     icefrac, slp, azFac, 1, avdrho, Tb, patm)
  call jsubv(NS, zz, latitude*d2r, albedo0, thIn, pfrost, &
       &     nz,   rhoc, fracIR, fracDust, Fgeotherm,    dt, zfac, &
       &     icefrac, slp, azFac, 0, avdrho, Tb, patm)
  print *, 'ice depth: ','  rho_ice-rho_surf'
  print *, zz,'#',avdrho
  where (avdrho>=0.) zdepth= -9999.  ! no ice
  if (minval(avdrho(2:NS))<0.) then  ! stable somewhere
     zz(1)=zdepth0;  zz(2:NS) = spread(delta/4.,1,NS-1)
     call rtbisv(NS,zz,(/ zdepth0, spread(zmax,1,NS-1) /), &
          &     zacc,avdrho(:), &
          &     latitude*d2r,albedo0,thIn,pfrost,nz,rhoc,fracIR, &
          &     fracDust,Fgeotherm,dt,zfac,icefrac,slp(:),azFac(:),zdepth(:),patm)
  endif
  zdepth(1) = zdepth0
  print *,'Equilibrium ice table depths= ',zdepth

80 continue
  open(unit=33,file=outfile,status='unknown',action='write')
  do k=1,NS
     write(33,'(i5,1x,f7.2,1x,f7.3,2x,f0.3,1x,2(f7.1,2x),f5.2,2x,f6.1,2x,f10.4)') &
          &     line_nr,lon,latitude,albedo0,thIn,Tfrost,slp(k)/d2r,azFac(k)/d2r,zdepth(k)
  enddo
  close(33)
  
90 continue
end program mars_mapi2p




function rtbis(x1,x2,xacc,fmid, &
     &     latitude,albedo0,thIn,pfrost,nz,rhoc,fracIR,fracDust, &
     &     Fgeotherm,dt,zfac,icefrac,surfaceSlope,azFac,patm)
  ! finds root with bisection method a la Numerical Recipes (C)
  implicit none
  REAL*8 rtbis,x1,x2,xacc
  INTEGER, PARAMETER :: JMAX=40
  INTEGER j,nz
  REAL*8 dx,f,fmid,xmid,Tb,rhoc
  real*8 latitude,albedo0,thIn,pfrost,fracIR,fracDust,Fgeotherm,dt
  real*8 zfac,icefrac,surfaceSlope,azFac,patm
  real*8 xlower,xupper,fupper,flower
  
  Tb = -1.e32
  call jsubv(1,x1, &
       &     latitude,albedo0,thIn,pfrost,nz/2,rhoc,fracIR,fracDust, &
       &     Fgeotherm,2.*dt,zfac,icefrac,surfaceSlope,azFac,1,f,Tb,patm)
  call jsubv(1,x1, &
       &     latitude,albedo0,thIn,pfrost,nz,rhoc,fracIR,fracDust, &
       &     Fgeotherm,dt,zfac,icefrac,surfaceSlope,azFac,0,f,Tb,patm)
  print *,x1,f
  if (f<0.) then  ! ice stable at the uppermost location
     rtbis = x1   ! equilibrium ice table is less than x1
     fmid = f
     return
  endif
  if(f*fmid >= 0.) stop 'root must be bracketed in rtbis'
  rtbis = x2
  dx = x1-x2
  xupper=x1; fupper=f
  xlower=x2; flower=fmid
  do j=1,JMAX
     dx = dx*.5
     xmid = rtbis+dx
     Tb = -1.e32
     call jsubv(1,xmid, &
          &     latitude,albedo0,thIn,pfrost,nz/2,rhoc,fracIR,fracDust, &
          &     Fgeotherm,2.*dt,zfac,icefrac,surfaceSlope,azFac,1,fmid,Tb,patm)
     call jsubv(1,xmid, &
          &     latitude,albedo0,thIn,pfrost,nz  ,rhoc,fracIR,fracDust, &
          &     Fgeotherm,   dt,zfac,icefrac,surfaceSlope,azFac,0,fmid,Tb,patm)
     print *,xmid,fmid
     if(fmid <= 0.) then
        rtbis = xmid
        xlower=xmid; flower=fmid
     else
        xupper=xmid; fupper=fmid
     endif
     if(abs(dx/xmid)<xacc .or. fmid==0.) then
        
!-------do linear interpolation at last
        rtbis = (fupper*xlower - flower*xupper)/(fupper-flower)

!-------report last stable ice table instead
!        rtbis = xlower
!        fmid = flower
        return
     endif
  enddo
  print *,'too many bisections in rtbis'
end function rtbis



      
subroutine rtbisv(NS,x1,x2,xacc,fmid, &
     &     latitude,albedo0,thIn,pfrost,nz,rhoc,fracIR,fracDust, &
     &     Fgeotherm,dt,zfac,icefrac,surfaceSlope,azFac,zdepth,patm)
!***********************************************************************
! finds roots with bisection method where a root exists
!
! first slope is special
!                      ice depth must be contained in x1(1) and x2(1)
!                      if -9999 set xmid=x2(2) = zmax, then leave constant 
!                               x1(1) cannot be overwritten, that's okay
!                               return zdepth(1)=-9999
!                      if finite set xmid, then leave constant
!                               return zdepth(1)=x1(1)
!                      not clear what would happen if called with NS=1
!                               will stop
!
! fupper>0 & flower>0  case I: ice unstable top and bottom
!                      set xmid = zmax = x2(2), then leave constant
!                      x1(1) cannot be overwritten, that's okay
!                      return zdepth = -9999
!
! fupper>0 & flower<0  case II: root bracketed
!                      determine root, set xmid and fmid, 
!                      return zdepth
!
! fupper<0 & flower<0  case III: ice stable top and bottom
!                      set zdepth to x1, set xmid to x1, then leave constant
!                      return x1
!
! fupper<0 & flower>0  physically impossible case: program will stop
!***********************************************************************
  implicit none
  integer, intent(IN) :: NS, nz
  real*8, intent(IN) :: x1(NS),x2(NS),xacc,latitude,albedo0,thIn
  real*8, intent(INOUT) :: fmid(NS)
  real*8, intent(OUT) :: zdepth(NS)
  real*8, intent(IN) :: pfrost,rhoc,fracIR,fracDust,Fgeotherm,dt
  real*8, intent(IN) :: zfac,icefrac,surfaceSlope(NS),azFac(NS),patm
  INTEGER, PARAMETER :: JMAX=40
  INTEGER j, k
  REAL*8, dimension(NS) :: dx,f,xmid,Tb,xlower,xupper,fupper,flower
  real*8 zmax
  zmax = x2(2)

  ! flat slope (index 1)
  if (x1(1)/=x2(1)) stop 'rtbisv: incorrect input'
  if (NS<=1) stop 'rtbisv: NS needs to be larger than 1'
  if (x1(1)<0.) then
     xmid = zmax
     dx(1) = 0.
     zdepth(1) = -9999.
  else
     xmid = x1(1)
     dx(1) = 0.
     zdepth(1) = x1(1)
  endif

  ! all other slopes (index 2:NS)
  dx = 0.  ! satisfies accuracy requirement
  do k=2,NS
     if(fmid(k)>0.) then ! case I: ice unstable at bottom
        zdepth(k) = -9999.
        xmid(k) = zmax
     end if
  enddo
  if (minval(fmid)>0.) return  ! nothing left to do
  Tb(:) = -1.e32
  call jsubv(NS,x1, &
       &     latitude,albedo0,thIn,pfrost,nz/2,rhoc,fracIR,fracDust, &
       &     Fgeotherm,2.*dt,zfac,icefrac,surfaceSlope,azFac,1,f,Tb,patm)
  call jsubv(NS,x1, &
       &     latitude,albedo0,thIn,pfrost,nz  ,rhoc,fracIR,fracDust, &
       &     Fgeotherm,   dt,zfac,icefrac,surfaceSlope,azFac,0,f,Tb,patm)
  print *,x1,'#',f
  do k=2,NS
     if(f(k)<0.) then  ! case III: ice stable at the uppermost location
        zdepth(k) = x1(k)
        fmid(k) = f(k)
        xmid(k) = x1(k)
        dx(k) = 0.
     endif
  enddo
  if (maxval(f)<0.) return  ! ice stable at the uppermost location everywhere
  if(minval(f*fmid) >= 0.) stop 'root must be bracketed in rtbis'
  xupper=x1; fupper=f
  xlower=x2; flower=fmid
  do k=2,NS
     if (fupper(k)*flower(k)<0.) then  ! case II
        zdepth(k) = x2(k)
        dx(k) = x1(k)-x2(k)
     end if
     if (fmid(k)>0. .and. f(k)<0.) stop 'rtbisv: impossible case'
  enddo
  do j=1,JMAX
     do k=2,NS
        if (fupper(k)*flower(k)<0.) then ! case II
           dx(k) = dx(k)*.5
           xmid(k) = zdepth(k)+dx(k)
        endif
     end do
     Tb(:) = -1.e32
     call jsubv(NS,xmid, &
          &     latitude,albedo0,thIn,pfrost,nz/2,rhoc,fracIR,fracDust, &
          &     Fgeotherm,2.*dt,zfac,icefrac,surfaceSlope,azFac,1,fmid,Tb,patm)
     call jsubv(NS,xmid, &
          &     latitude,albedo0,thIn,pfrost,nz  ,rhoc,fracIR,fracDust, &
          &     Fgeotherm,   dt,zfac,icefrac,surfaceSlope,azFac,0,fmid,Tb,patm)
     print *,xmid,'#',fmid
     do k=2,NS
        if (fupper(k)*flower(k)<0.) then 
           if(fmid(k)<=0.) then   ! ice stable at xmid
              zdepth(k)=xmid(k)
              xlower(k)=xmid(k); flower(k)=fmid(k)
           else
              xupper(k)=xmid(k); fupper(k)=fmid(k)
           endif
        end if
     enddo
     if(maxval(abs(dx/xmid))<xacc) then ! .or. fmid.eq.0.) then
        do k=2,NS
           if (fupper(k)*flower(k)<0.) then
              ! do linear interpolation at last
              zdepth(k) = (fupper(k)*xlower(k) - flower(k)*xupper(k))/  &
                   &             (fupper(k)-flower(k))
              
              ! report last stable ice table instead
              ! zdepth(k) = xlower(k)
              ! fmid(k) = flower(k)
           endif
        enddo
        return
     endif
  enddo
  print *,'too many bisections in rtbisv'
end subroutine rtbisv


