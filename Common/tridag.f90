PURE SUBROUTINE tridag(a,b,c,r,x,N)
  ! solves Ax = r, where A is a tridiagonal matrix consisting of vectors a, b, c
  ! N ... number of equations
  ! a ... subdiagonal, index 1 is not used
  ! b ... main diagonal, indexed from 1,...,N
  ! c ... superdiagonal, index N is not used
  ! Method: A=LU, first solve Ls=r and then Ux=s
  ! stable when |b_i| > |a_i|+|c_i|
  implicit none
  INTEGER, intent(IN) :: N
  REAL*8, intent(IN) :: a(N),b(N),c(N),r(N)
  REAL*8, intent(OUT) :: x(N)
  INTEGER j
  REAL*8 bet,scratch(N)

  bet = b(1)
  x(1) = r(1)/b(1)
  do j = 2,N  ! downward sweep
     scratch(j) = c(j-1)/bet
     bet = b(j)-a(j)*scratch(j)  ! must never be 0
     x(j) = ( r(j)-a(j)*x(j-1) ) / bet
  end do
  do j = N-1,1,-1  ! upward sweep
     x(j) = x(j)-scratch(j+1)*x(j+1)
  end do
end SUBROUTINE tridag
