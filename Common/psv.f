      real*8 function psv(T)
C     saturation vapor pressure of H2O ice [Pascal]
C     input is temperature [Kelvin]
      implicit none
      real*8 T

C-----parametrization 1
c      real*8 DHmelt,DHvap,DHsub,R,pt,Tt,C
c      parameter (DHmelt=6008.,DHvap=45050.)
c      parameter (DHsub=DHmelt+DHvap) ! sublimation enthalpy [J/mol]
c      parameter (R=8.314,pt=6.11e2,Tt=273.16)
c      C = (DHsub/R)*(1./T - 1./Tt)
c      psv = pt*exp(-C)

C-----parametrization 2
C     eq. (2) in Murphy & Koop, Q. J. R. Meteor. Soc. 131, 1539 (2005)
C     differs from parametrization 1 by only 0.1%
      real*8 A,B
      parameter (A=-6143.7, B=28.9074)
      psv = exp(A/T+B)  ! Clapeyron

C-----parametrization 3      
C     eq. (7) in Murphy & Koop, Q. J. R. Meteor. Soc. 131, 1539 (2005)
c     psv = exp(9.550426 - 5723.265/T + 3.53068*log(T) - 0.00728332*T)
      
      end


      
      real*8 function frostpoint(p)
C     inverse of psv
C     input is partial pressure [Pascal]
C     output is temperature [Kelvin]
      implicit none
      real*8 p
      
C-----inverse of parametrization 1
c      real*8 DHmelt,DHvap,DHsub,R,pt,Tt
c      parameter (DHmelt=6008.,DHvap=45050.)
c      parameter (DHsub=DHmelt+DHvap)
c      parameter (R=8.314,pt=6.11e2,Tt=273.16)
c      frostpoint = 1./(1./Tt-R/DHsub*log(p/pt))
      
C-----inverse of parametrization 2
C     inverse of eq. (2) in Murphy & Koop (2005)
      real*8 A,B
      parameter (A=-6143.7, B=28.9074)
      frostpoint = A / (log(p) - B)

C-----approximate inverse of parametrization 3
C     eq. (8) in Murphy & Koop (2005)
c      frostpoint = (1.814625*log(p) + 6190.134)/(29.120 - log(p))
      
      end
