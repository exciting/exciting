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
      double precision dtmp, optval
c
c     .. Executable statements ..
c
      j = 1
      if ( order .le. 0 ) then
         optval =  1.d99
      else
         optval = -1.d99
      end if
      if (order.eq.0) then
c
c...        Nearest to target
c
         do i=1,n
            if ( b(i) .ne. 0.d0 ) then
               dtmp = abs(a(i)/b(i)-sa/sb)
               if (dtmp.lt.optval) then
                  j=i
                  optval = dtmp
               end if
            end if
         enddo
      elseif (order.eq.-1) then
c
c...        Smallest real part
c
         do i=1,n
            if ( b(i) .eq. 0.d0 ) then
               if ( dble(a(i)) .lt. 0.d0 ) then
                  j=i
                  optval = -1.d99
               end if
            else
               dtmp = dble(a(i)/b(i))
               if (dtmp.lt.optval) then
                  j=i
                  optval = dtmp
               end if
            end if
         enddo
      elseif (order.eq.1) then
c
c...        Largest real part
c
         do i=1,n
            if ( b(i) .eq. 0.d0 ) then
               if ( dble(a(i)) .gt. 0.d0 ) then
                  j=i
                  optval = 1.d99
               end if
            else
               dtmp = dble(a(i)/b(i))
               if (dtmp.gt.optval) then
                  j=i
                  optval = dtmp
               end if
            end if
         enddo
      elseif (order.eq.-2) then
c
c...        Smallest imaginari part
c
         do i=1,n
            if ( b(i) .eq. 0.d0 ) then
               if ( imag(a(i)) .lt. 0.d0 ) then
                  j=i
                  optval = -1.d99
               end if
            else
               dtmp = imag(a(i)/b(i))
               if (dtmp.lt.optval) then
                  j=i
                  optval = dtmp
               end if
            end if
         enddo
      elseif (order.eq.2) then
c
c...        Largest imaginari part
c
         do i=1,n
            if ( b(i) .eq. 0.d0 ) then
               if ( imag(a(i)) .gt. 0.d0 ) then
                  j=i
                  optval = 1.d99
               end if
            else
               dtmp = imag(a(i)/b(i))
               if (dtmp.gt.optval) then
                  j=i
                  optval = dtmp
               end if
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
