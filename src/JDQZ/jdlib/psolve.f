      subroutine psolve (n, x, nq, q, kz, invqkz, ldqkz, ipiv, f)
c
c     Coded by Diederik Fokkema
c
c     .. Parameters ..
c
      implicit none
      integer n, nq, ldqkz, ipiv(*)
      double complex x(*), q(n,*), kz(n,*),
     $     invqkz(ldqkz,*), f(*)
c
c     .. local .. 
c
      integer info
      double complex zero, one
      parameter (zero=(0.0d0,0.0d0), one=(1.0d0,0.0d0))
c
c     .. Executable Statements ..
c
      call precon( n, x )
      call zgemv ('c', n, nq, one, q, n, x, 1, zero, f, 1)
      call zgetrs('n', nq, 1, invqkz, ldqkz, ipiv, f, ldqkz, info)
      call zgemv ('n', n, nq, -one, kz, n, f, 1, one, x, 1)
c
c     --- Return
c
      end


