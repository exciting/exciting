!
!
!
!
Module m_genkerr
  Implicit None
Contains
  !
  !
  Subroutine genkerr (w, eps, kerr)
    Use mod_constants, Only: pi, zi
    Use modxs
    Implicit None
    ! arguments
    Real (8), Intent (In) :: w (:)
    Complex (8), Intent (In) :: eps (:,:,:)
    Complex (8), Intent (Out) :: kerr (:)
    ! local variables
    Character (*), Parameter :: thisnam = 'genkerr'
    Complex(8), Allocatable :: sqeps11(:)
    Integer :: n1 (1), n, iw, sig
    n1 = shape (w)
    n = n1 (1)
         
    allocate(sqeps11(n))

    !select lower branch of sqrt
    Do iw = 2, n
       sqeps11(iw) = sqrt(eps(1,1,iw))
       If (aimag(sqeps11(iw)).gt.0.0d0) Then
          sqeps11(iw) = -sqeps11(iw)
       End If
    End Do
    sqeps11(1) = cmplx(1.0d0,0.0d0,8)
    
    kerr (:) = -eps(1,2,:)/((eps(1,1,:)-1.0d0) * sqeps11(:))

    deallocate(sqeps11)
  End Subroutine genkerr
  !
End Module m_genkerr
