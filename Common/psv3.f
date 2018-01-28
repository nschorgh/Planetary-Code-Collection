C     updated parametrization for saturation vapor pressure of H2O ice
C     results differ by only 0.1% from those in psv.f,
C     but the equation has a reference


      real*8 function psv(T)
C     equilibrium pressure of H2O ice (Pascal)
C     input is temperature (Kelvin)
C     Murphy & Koop, Q. J. R. Meteor. Soc. 131, 1539 (2005)
      implicit none
      real*8 T,A,B
      parameter (A=-6143.7, B=28.9074)
      psv = exp(A/T+B);
      end


      real*8 function frostpoint(p)
C     T=frostpoint(T); inverse of psv
      implicit none
      real*8 p,A,B
      parameter (A=-6143.7, B=28.9074)
      frostpoint = A / (log(p) - B)
      end
