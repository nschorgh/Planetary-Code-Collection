      real*8 function psv(T)
C     saturation vapor pressure along solid-vapor line of H20
      implicit none
      real*8 T,DHmelt,DHvap,DHsub,R,pt,Tt,C
      parameter (DHmelt=6008.,DHvap=45050.)
      parameter (DHsub=DHmelt+DHvap) ! sublimation enthalpy (J/mol)
      parameter (R=8.314,pt=6.11e2,Tt=273.16)
      C=(DHsub/R)*(1./T - 1./Tt)
      psv=pt*exp(-C)
      end


      real*8 function frostpoint(p)
C     T=frostpoint(p); inverse of psv
      real*8 p,DHmelt,DHvap,DHsub,R,pt,Tt
      parameter (DHmelt=6008.,DHvap=45050.)
      parameter (DHsub=DHmelt+DHvap)
      parameter (R=8.314,pt=6.11e2,Tt=273.16)
      frostpoint=1./(1./Tt-R/DHsub*log(p/pt))
      end
