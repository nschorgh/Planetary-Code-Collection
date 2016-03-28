C================================================
C Tridiagonal solver
C================================================
      SUBROUTINE tridag(a,b,c,r,u,n)
      INTEGER n,NMAX
      REAL*8 a(n),b(n),c(n),r(n),u(n)
      PARAMETER (NMAX=1000)
      INTEGER j
      REAL*8 bet,gam(NMAX)
      if(b(1).eq.0.) then
         print *, 'tridag: rewrite equations'
         stop
      endif 
      if(n.gt.NMAX) then
         print *, 'tridag: too many points, set NMAX>',n
     &        ,' and recompile'
         stop
      endif 
      bet=b(1)
      u(1)=r(1)/bet
      do 11 j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
c        if(bet.eq.0.)pause 'tridag failed'
        u(j)=(r(j)-a(j)*u(j-1))/bet
11    continue
      do 12 j=n-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
12    continue
      return
      END

