      real*8 function psv(T)
C     equilibrium pressure of H2O ice (Pascal)
C     input is temperature (Kelvin)
      implicit none
      real*8 T,A,B
      parameter (A=-6143.7, B=28.9074)
c     eq. (2) in Murphy & Koop, Q. J. R. Meteor. Soc. 131, 1539 (2005) 
      psv = exp(A/T+B)  ! Clapeyron
c     differs from psv.f by only 0.1%
      end


      real*8 function frostpoint(p)
C     inverse of psv
      implicit none
      real*8 p,A,B
      parameter (A=-6143.7, B=28.9074)
      frostpoint = A / (log(p) - B)
      end
