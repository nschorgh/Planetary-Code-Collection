
subroutine deriv1(z,nz,y,y0,yNp1,yp)
  ! first derivative of a function y(z) on irregular grid
  ! upper b.c.: y(0)=y0
  ! lower b.c.: yp =0
  implicit none
  integer, intent(IN) :: nz
  real(8), intent(IN) :: z(nz),y(nz),y0,yNp1
  real(8), intent(OUT) :: yp(nz)
  integer j
  real(8) hm,hp,c1,c2,c3

  hp = z(2)-z(1)
  !hm = z(1)
  c1 = z(1)/(hp*z(2))
  c2 = 1/z(1) - 1/(z(2)-z(1))
  c3 = -hp/(z(1)*z(2))
  yp(1) = c1*y(2) + c2*y(1) + c3*y0
  do j=2,nz-1
     hp = z(j+1)-z(j)
     hm = z(j)-z(j-1)
     c1 = +hm/(hp*(z(j+1)-z(j-1)))
     c2 = 1/hm - 1/hp
     c3 = -hp/(hm*(z(j+1)-z(j-1)))
     yp(j) = c1*y(j+1) + c2*y(j) + c3*y(j-1)
  enddo
  yp(nz) = (yNp1 - y(nz-1))/(2.*(z(nz)-z(nz-1)))
end subroutine deriv1



subroutine deriv2_full(z,nz,a,b,a0,b0,bNp1,yp2)
  ! second derivative (a*b_z)_z on irregular grid
  ! upper b.c.: a(0)=a0, b(0)=b0
  ! 2nd order, without cross-terms
  implicit none
  integer, intent(IN) :: nz
  real(8), intent(IN) :: z(nz),a(nz),b(nz),a0,b0,bNp1
  real(8), intent(OUT) :: yp2(nz)
  integer j
  real(8) hm,hp,c1,c2,c3
  
  c1 = 1./(z(1)*z(2))
  c2 = 1./((z(2)-z(1))*z(1))
  c3 = 1./((z(2)-z(1))*z(2))
  yp2(1) = -c2*a(1)*b(1) &   
       &  -c1*a0*b(1) + c1*(a(1)+a0)*b0 &
       &  -c3*a(2)*b(1) + c3*(a(1)+a(2))*b(2)
  do j=2,nz-1
     hp = z(j+1)-z(j)
     hm = z(j)-z(j-1)
     c1 = 1./(hm*(z(j+1)-z(j-1)))
     c2 = 1./(hp*hm)
     c3 = 1./(hp*(z(j+1)-z(j-1)))
     yp2(j) = -c2*a(j)*b(j) &   
          &  -c1*a(j-1)*b(j) + c1*(a(j)+a(j-1))*b(j-1) &
          &  -c3*a(j+1)*b(j) + c3*(a(j)+a(j+1))*b(j+1)
  enddo

  ! lower bc: b_z = 0
  !yp(nz)= 2.*a(nz)*(b(nz-1)-b(nz))/(z(nz)-z(nz-1))**2
  !more general lower bc
  yp2(nz) = a(nz)*(bNp1 - 2*b(nz) + b(nz-1))/(z(nz)-z(nz-1))**2
end subroutine deriv2_full



subroutine deriv2_simple(z,nz,y,y0,yNp1,yp2)
  ! second derivative y_zz on irregular grid
  ! boundary conditions: y(0)=y0, y(nz+1)=yNp1
  implicit none
  integer, intent(IN) :: nz
  real(8), intent(IN) :: z(nz),y(nz),y0,yNp1
  real(8), intent(OUT) :: yp2(nz)
  integer j
  real(8) hm,hp,c1,c2,c3

  c1 = +2./((z(2)-z(1))*z(2))
  c2 = -2./((z(2)-z(1))*z(1))
  c3 = +2./(z(1)*z(2))
  yp2(1) = c1*y(2) + c2*y(1) + c3*y0

  do j=2,nz-1
     hp = z(j+1)-z(j)
     hm = z(j)-z(j-1)
     c1 = +2./(hp*(z(j+1)-z(j-1)))
     c2 = -2./(hp*hm)
     c3 = +2./(hm*(z(j+1)-z(j-1)))
     yp2(j) = c1*y(j+1) + c2*y(j) + c3*y(j-1)
  enddo

  yp2(nz) = (yNp1 - 2*y(nz) + y(nz-1))/(z(nz)-z(nz-1))**2
end subroutine deriv2_simple



real(8) function deriv1_onesided(j,z,nz,y)
  ! first derivative of function y(z) at z(j)
  ! one-sided derivative on irregular grid
  implicit none
  integer, intent(IN) :: nz,j
  real(8), intent(IN) :: z(nz),y(nz)
  real(8) h1,h2,c1,c2,c3
  if (j<1 .or. j>nz-2) then
     deriv1_onesided = -9999.
  else
     h1 = z(j+1)-z(j)
     h2 = z(j+2)-z(j+1)
     c1 = -(2*h1+h2)/(h1*(h1+h2))
     c2 =  (h1+h2)/(h1*h2)
     c3 = -h1/(h2*(h1+h2))
     deriv1_onesided = c1*y(j) + c2*y(j+1) + c3*y(j+2)
  endif
end function deriv1_onesided

