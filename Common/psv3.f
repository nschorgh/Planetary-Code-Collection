      real*8 function psv(T)
C     equilibrium pressure of H2O ice (Pascal)
C     input is temperature (Kelvin)
      implicit none
      real*8 T
C     eq. (7) in Murphy & Koop, Q. J. R. Meteor. Soc. 131, 1539 (2005)
      psv = exp(9.550426 - 5723.265/T + 3.53068*log(T) - 0.00728332*T)
      end


      real*8 function frostpoint(p)
C     inverse of psv
      implicit none
      real*8 p
C     eq. (8) in Murphy & Koop (2005)
      frostpoint = (1.814625*log(p) + 6190.134)/(29.120 - log(p))
      end
