real(8) function psvco2(T)
  ! solid-vapor transition for CO2
  ! returns saturation pressure in Pascal, input is temperature in Kelvin
  implicit none
  real(8), intent(IN) :: T

  psvco2 = exp(23.3494 - 3182.48/T)*100.
end function psvco2


real(8) function tfrostco2(p)
  ! the inverse of function psvco2
  ! input is pressure in Pascal, output is temperature in Kelvin
  implicit none
  real(8), intent(IN) :: p
  tfrostco2 = 3182.48/(23.3494+log(100./p))
end function tfrostco2




!------------------------------------------------------------------------------
! Antoine equation parameters from NIST, 154K-196K 
! based on Giauque and Egan (1937)
! A=6.81228, B=1301.679, C=-3.494
! p = 10**(A-(B/(T+C)))*1.e5

! Expressions from Int. Crit. Tabl. 3.207, based on many references
! mm2Pa = 133.32
! -135C to -56.7C (138K-216K)
! p = 10**(-0.05223*26179.3/T + 9.9082)*mm2Pa
! -183C to -135C (90K-138K)
! p = 10**(-1275.62/T + 0.006833*T + 8.3071)*mm2Pa
! Expressions from Int. Crit. Tabl. 3.208, based on Henning
! -110 to -80C (163K-193K)
! p = 10**(- 1279.11/T + 1.75*log10(T) - 0.0020757*T + 5.85242)*mm2Pa
! p = 10**(- 1352.6/T + 9.8318)*mm2Pa
      
! Mars book (1992), p959, fit by chapter authors
! p = exp(23.3494 - 3182.48/T)*100.   ! 120-160K
! p = exp(25.2194 - 3311.57/T - 6.71e-3*T)*100  ! 100-170K
! Mars book (1992), p960, based on Miller & Smythe (1970)
! p = exp(26.1228 - 3385.26/T - 9.4461e-3*T)*100 ! 123-173K

! Fray & Schmitt, PSS 57, 2053 (2009)
! A0=1.476e1, A1=-2.571e3, A2=-7.781e4, A3=4.325e6, A4=-1.207e8, A5=1.350e9
! p = exp(A0+A1/T+A2/T**2+A3/T**3+A4/T**4+A5/T**5)*1e5 ! 40K-194.7K
! A0=1.861e1, A1=-4.154e3, A2=1.041e5
! p = exp(A0 + A1/T + A2/T**2)*1e5  ! 194.7K-216.58K
!------------------------------------------------------------------------------
