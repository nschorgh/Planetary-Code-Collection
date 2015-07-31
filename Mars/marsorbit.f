C================================================
C Subroutine to return various orbital parameters
C
C INPUTS: 
C dt0 is the reference dt_J2000
C tj is the time (in days) relative to the reference day
C
C How to call
C     jd=dble(julday(imm,iday,iyr))  !  JD for noon UTC on iyear/imm/iday
C     temp1 = (jd-2451545.d0)/36525.d0
C     dcor = (64.184d0 + 95.*temp1 + 35.*temp1**2) ! correction in sec
C     dt0_j2000 = jd + dcor/earthDay - 2451545.d0
C     call marsorbit(dt0_j2000,0.d0,Ls,dec,r)
C
C OUTPUTS:
C Ls = areocentric longitude (radians)
C dec = planetocentric solar declination (radians)
C r = Mars heliocentric distance (AU)
C
C Data from Allison & McEwen, 2000
C Written by Samar Khatiwala, 2001
C================================================

      SUBROUTINE marsorbit(dt0,tj,Ls,dec,r)
      implicit none
      real*8 dt0,tj,Ls,dec,r

C     Some constants
      real*8 pi,d2r
      parameter (pi=3.1415926535897932,d2r=pi/180.d0)
      integer i
      real*8 alpha, PBS, M, dj, eps
      real*8 A(7),tau(7),phi(7),c1,c2,c3,c4,c5
    
      data A/0.007d0,0.006d0,0.004d0,0.004d0,0.002d0,0.002d0,0.002d0/
      data tau/2.2353d0,2.7543d0,1.1177d0,15.7866d0,2.1354d0,2.4694d0,
     &         32.8493d0/
C      data phi/d2r*49.409d0,d2r*168.173d0,d2r*191.837d0,d2r*21.736d0,
C     &         d2r*15.704d0,d2r*95.528d0,d2r*49.095d0/
      data phi/0.8623497301178783d0,2.9351725629564238d0,
     &         3.3481872771483618d0,0.3793647662134875d0,
     &         0.2740865057331895d0,1.6672781278451432d0,
     &         0.8568693962666161d0/
      data c1/0.0172024188932616d0/   !  d2r*0.985626d0
      data c2/0.3383669820841407d0/   !  d2r*19.3870d0
      data c3/0.0091458874362701d0/   !  d2r*0.52402075d0
      data c4/0.4396833451624115d0/   !  d2r*25.192d0
      data c5/6.021385919380437d-09/  !  d2r*3.45d-7

      dj = dt0 + tj ! actual dt_j2000


      PBS=0.d0
      do i=1,7
         PBS = PBS + A(i)*cos(c1*dj/tau(i) + phi(i))
      enddo

      M=c2 + c3*dj  ! mean anomaly in radians
      alpha=270.3863d0 + 0.52403840d0*dj ! fictitious mean sun angle in degrees

C     Ls in degrees      
      Ls = alpha + (10.691d0 + 3.d-7*dj)*sin(M) + 0.623d0*sin(2*M) + 
     &     0.05d0*sin(3*M) + 0.005d0*sin(4*M) + 0.0005d0*sin(5*M) + PBS
      Ls=mod(Ls,360.d0)*d2r  !  Ls in radians
      eps = c4 + c5*dj  ! obliquity in radians
c      eps = c4 ! CONSTANT obliquity in radians
      dec = asin(sin(eps)*sin(Ls))  !  declination in radians
      r = 1.5236d0*(1.00436d0 - 0.09309d0*cos(M) - 0.00436d0*cos(2*M) - 
     &    0.00031d0*cos(3*M))  !  distance in AU


      END

