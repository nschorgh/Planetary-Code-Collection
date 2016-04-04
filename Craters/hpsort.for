      SUBROUTINE hpsort(n,ra,ind)
      implicit none
      INTEGER, intent(IN) :: n
      REAL(8), intent(INOUT) :: ra(n)
      INTEGER, intent(OUT) :: ind(n)
      INTEGER rind
      INTEGER i,ir,j,l
      REAL(8) rra
      if (n==1) ind(1)=1
      if (n.lt.2) return
      l=n/2+1
      ir=n
      do i=1,n
         ind(i)=i
      enddo
10    continue
        if(l.gt.1)then
          l=l-1
          rra=ra(l); rind=ind(l)
        else
          rra=ra(ir); rind=ind(ir)
          ra(ir)=ra(1); ind(ir)=ind(1)
          ir=ir-1
          if(ir.eq.1)then
            ra(1)=rra; ind(1)=rind
            return
          endif
        endif
        i=l
        j=l+l
20      if(j.le.ir)then
          if(j.lt.ir)then
            if(ra(j).lt.ra(j+1))j=j+1
          endif
          if(rra.lt.ra(j))then
            ra(i)=ra(j); ind(i)=ind(j)
            i=j
            j=j+j
          else
            j=ir+1
          endif
        goto 20
        endif
        ra(i)=rra; ind(i)=rind
      goto 10
      END
C  (C) Copr. 1986-92 Numerical Recipes Software 0(9p#31&#5(+.
