      integer function select (n, sa, sb, a, b, order)
c
c     Coded by Diederik Fokkema
c     Modified Martin van Gijzen, test on division by zero
c
c     .. Parameters ..
c
      implicit none
      integer n, order
      double complex sa, sb, a(*), b(*)
c
c     .. Local ..
c
      integer i, j
      double precision dtmp
c
c     .. Executable statements ..
c
      j = 1
      if (order.eq.0) then
         do i=1,n
            if ( b(i) .ne. 0.d0 ) then
               dtmp = abs(a(i)/b(i)-sa/sb)
               if (dtmp.lt.abs(a(j)/b(j)-sa/sb)) j=i
            end if
         enddo
      elseif (order.eq.-1) then
         do i=1,n
            if ( b(i) .eq. 0.d0 ) then
               if ( dble(a(i)) .lt. 0.d0 ) j=i
            else
               dtmp = dble(a(i)/b(i))
               if (dtmp.lt.dble(a(j)/b(j))) j=i
            end if
         enddo
      elseif (order.eq.1) then
         do i=1,n
            if ( b(i) .eq. 0.d0 ) then
               if ( dble(a(i)) .gt. 0.d0 ) j=i
            else
               dtmp = dble(a(i)/b(i))
               if (dtmp.gt.dble(a(j)/b(j))) j=i
            end if
         enddo
      elseif (order.eq.-2) then
         do i=1,n
            if ( b(i) .eq. 0.d0 ) then
               if ( imag(a(i)) .lt. 0.d0 ) j=i
            else
               dtmp = imag(a(i)/b(i))
               if (dtmp.lt.imag(a(j)/b(j))) j=i
            end if
         enddo
      elseif (order.eq.2) then
         do i=1,n
            if ( b(i) .eq. 0.d0 ) then
               if ( imag(a(i)) .gt. 0.d0 ) j=i
            else
               dtmp = imag(a(i)/b(i))
               if (dtmp.gt.imag(a(j)/b(j))) j=i
            end if
         enddo
      else
         call error ('unknown order in select')
      endif
      select = j
c
c     --- Return
c
      end
