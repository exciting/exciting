!BOP
!
! !ROUTINE: gaussj
!
! !INTERFACE:
      subroutine gaussj(a,n,np,b,m,mp)
      
! !DESCRIPTION:
!      
! Linear equation solution by Gauss-Jordan elimination, equation
!\ref{gaussj-1} below.
!
!\begin{equation}\label{gaussj-1}
!\begin{aligned}
!  \begin{bmatrix}
!    a_{11}&a_{12}&a_{13}&a_{14}\\
!    a_{21}&a_{22}&a_{23}&a_{24}\\
!    a_{31}&a_{32}&a_{33}&a_{34}\\
!    a_{41}&a_{42}&a_{43}&a_{44}
!  \end{bmatrix}\cdot &
!  \begin{bmatrix}
!    \begin{pmatrix}
!      x_{11}\\
!      x_{21}\\
!      x_{31}\\
!      x_{41}
!    \end{pmatrix}\sqcup
!    \begin{pmatrix}
!      x_{12}\\
!      x_{22}\\
!      x_{32}\\
!      x_{42}
!    \end{pmatrix}\sqcup
!    \begin{pmatrix}
!      x_{13}\\
!      x_{23}\\
!      x_{33}\\
!      x_{43}
!    \end{pmatrix}\sqcup
!    \begin{pmatrix}
!      y_{11}&y_{12}&y_{13}&y_{14}\\
!      y_{21}&y_{22}&y_{23}&y_{24}\\
!      y_{31}&y_{32}&y_{33}&y_{34}\\
!      y_{41}&y_{42}&y_{43}&y_{44}
!    \end{pmatrix}
!  \end{bmatrix}\\
!= & \begin{bmatrix}
!    \begin{pmatrix}
!      b_{11}\\
!      b_{21}\\
!      b_{31}\\
!      b_{41}
!    \end{pmatrix}\sqcup
!    \begin{pmatrix}
!      b_{12}\\
!      b_{22}\\
!      b_{32}\\
!      b_{42}
!    \end{pmatrix}\sqcup
!    \begin{pmatrix}
!      b_{13}\\
!      b_{23}\\
!      b_{33}\\
!      b_{43}
!    \end{pmatrix}\sqcup
!    \begin{pmatrix}
!      1 & 0 & 0 & 0 \\
!      0 & 1 & 0 & 0 \\
!      0 & 0 & 1 & 0 \\
!      0 & 0 & 0 & 1 \\
!    \end{pmatrix}
!  \end{bmatrix}
!\end{aligned}
!\end{equation}
! \texttt{a(1:n,1:n)}
! is an imput matrix stored in an array of physical dimensions \texttt{np}
! by \texttt{np}.\texttt{b(1:n,1:m)} is an input matrix containing the
! \texttt{m} right-hand side vectors, stored in an array of physical
! dimensions \texttt{np} by \texttt{mp}. On output, \texttt{a(1:n,1:n)} is
! replaced by its matrix inverse, and \texttt{b(1:n,1:m)} is replaced by
! the corresponding set of solution vectors.
!
! !INPUT PARAMETERS:
      
      implicit none
 
      integer(4), intent(in) :: m
      integer(4), intent(in) :: mp
      integer(4), intent(in) :: n
      integer(4), intent(in) :: np
      
! !INPUT/OUTPUT PARAMETERS:      

      real(8), intent(inout) :: a(1:np,1:np)
      real(8), intent(inout) :: b(1:np,1:mp)

! !LOCAL VARIABLES:
      
      integer(4) :: i 
      integer(4) :: icol 
      integer(4) :: irow 
      integer(4) :: j 
      integer(4) :: k 
      integer(4) :: l 
      integer(4) :: ll 
      integer(4) :: indxc(1:np) ! Used for bookkeeping on the pivoting
      integer(4) :: indxr(1:np) ! Used for bookkeeping on the pivoting
      integer(4) :: ipiv(1:np)  ! Used for bookkeeping on the pivoting
      
      real(8) :: big
      real(8) :: dum
      real(8) :: pivinv

! !REVISION HISTORY:
!
! Original subroutine: gaussj.for (c) copr. 1986-92 numerical recipes
! software &30i..
! Last modified: 7th. Jul. 2005 by RGA
!
!EOP
!BOC 
      do j=1,n
        ipiv(j)=0
      enddo ! j
!      
!     This is the main loop over the columns to be reduced
!      
      do i=1,n
        big=0.0d0
!        
!        This is the outer loop of the search for a pivot element
!        
        do j=1,n
          if(ipiv(j).ne.1)then
            do  k=1,n
              if (ipiv(k).eq.0) then
                if (abs(a(j,k)).ge.big)then
                  big=abs(a(j,k))
                  irow=j
                  icol=k
                endif
              else if (ipiv(k).gt.1) then
                stop 'singular matrix in gaussj'
              endif
            enddo ! k
          endif
        enddo ! j
        ipiv(icol)=ipiv(icol)+1
!        
!       We now have the pivot element, so we interchange rows, if needed,
!       to put the pivot element on the diagonal. The columns are not
!       physically interchanged, only relabeled: indxc(i), the column of
!       the ith pivot element, is the ith column that is reduced, while
!       indxr(i) is the row in which that pivot   element was originally
!       located. If indxr(i) \neq indxc(i) there is an implied column
!       interchange. With this form of bookkeeping, the solution b's will
!       end up in the correct order, and the inverse matrix will be
!       scrambled by columns.
!                 
        if (irow.ne.icol) then
          do l=1,n
            dum=a(irow,l)
            a(irow,l)=a(icol,l)
            a(icol,l)=dum
          enddo ! l
          do l=1,m
            dum=b(irow,l)
            b(irow,l)=b(icol,l)
            b(icol,l)=dum
          enddo ! l
        endif
        indxr(i)=irow
        indxc(i)=icol
!        
!       We are now ready to divide the pivot row by the pivot element,
!       located at irow and icol.
!        
        if (a(icol,icol).eq.0.0d0) stop 'singular matrix in gaussj'
        pivinv=1./a(icol,icol)
        a(icol,icol)=1.
        do l=1,n
          a(icol,l)=a(icol,l)*pivinv
        enddo ! l
        do l=1,m
          b(icol,l)=b(icol,l)*pivinv
        enddo ! l
!        
!       Next, we reduce the rows...
!        
        do ll=1,n
          if(ll.ne.icol)then   ! ... except for the pivor one, of course.
            dum=a(ll,icol)
            a(ll,icol)=0.0d0
            do l=1,n
              a(ll,l)=a(ll,l)-a(icol,l)*dum
            enddo ! l
            do l=1,m
              b(ll,l)=b(ll,l)-b(icol,l)*dum
            enddo ! l
          endif
        enddo ! ll
      enddo ! i
!      
!     This is the end of the main loop over columns of the reduction. It
!     only remains to unscramble the solution in view of the column
!     interchanges. We do this by interchanging pairs of columns in the
!     reverse order that the permutation was bilt up.
!         
      do l=n,1,-1
        if(indxr(l).ne.indxc(l))then
          do k=1,n
            dum=a(k,indxr(l))
            a(k,indxr(l))=a(k,indxc(l))
            a(k,indxc(l))=dum
          enddo ! k
        endif
      enddo ! l

      return

      end subroutine gaussj
!EOC
