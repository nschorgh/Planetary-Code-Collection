      real*8 function psvco2(T)
C     solid-vapor transition for CO2
C     returns saturation pressure in Pascal, input is temperature in Kelvin
      implicit none
      real*8 T

      psvco2 = exp(23.3494 - 3182.48/T)*100.
      end


      real*8 function tfrostco2(p)
C     the inverse of function psvco2
C     input is pressure in Pascal, output is temperature in Kelvin
      implicit none
      real*8 p
c      tfrostco2 = 1301.679/(6.81228-log10(p*1.e-5))+3.494
c      tfrostco2 = -1352.6/(log10(p*0.01)-9.8318)
      tfrostco2 = 3182.48/(23.3494+log(100./p))
      end




C-----------------------------------------------------------------------
C Antoine equation parameters from NIST, 154-186K 
C based on Giauque and Egan, 1937
c A=6.81228; B=1301.679; C=-3.494;
c p=10.^(A-(B./(T+C)))*1.e5;
     

C Expressions from Int. Crit. Tabl. 3.207, based on many references
c psvco2 = 10**(-1352.6/T + 9.8318)*100.

C -135C to -56.7C (138K - 216K)
c p=10^(-0.05223*26179.3/T + 9.9082)*100.;

C -183C to -135C (90K - 138K)
c p=10^(-1275.62/T + 0.006833*T + 8.3071)*100;

C -110 to -80C (163K - 193K)
c p=10^(- 1279.11/T + 1.75*log10(T) - 0.0020757*T + 5.85242)*100;
c p=10^(- 1352.6/T + 9.8318)*100;

C Mars book values are systematically lower than Int. Crit. Tabl. values

C Mars book, p959
c p=exp(23.3494 - 3182.48/T)*100.;   % 120-160K
c p=exp(25.2194 - 3311.57/T - 6.71e-3*T)*100;  % 100-170K

C Mars book p960, based on Miller & Smythe (1970)
c p=exp(26.1228 - 3385.26/T - 9.4461e-3*T)*100;  % 123-173K
C-------------------------------------------------------------------------
