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
    ! optical conductivity

    kerr (:) = -eps(1,2,:)/((eps(1,1,:)-1.0d0) * Sqrt(eps(1,1,:)))

  End Subroutine genkerr
  !
End Module m_genkerr
