!BOP
! !ROUTINE: sorteval
! !INTERFACE:
!
Module m_sorteval
      Implicit None
Contains
Subroutine sorteval (vpl, ra, ca)
! !USES:
      Use modmain
! !DESCRIPTION:
!   Rearrange real array {\tt ra} or complex array {\tt ca} through
!   the same permutations that would sort
!   the energy eigenvalues into ascending order for a given k point {\tt vpl}.
!   No sorting of the eigenvalues itself is performed.
!   Uses the sortidx subroutine which in turn uses the Heapsort algorithm.
! !REVISION HISTORY:
!   Created Jul 2013 SR
!EOP
!BOC
       ! arguments
    Real (8), Intent (In) :: vpl (3)
    Real (8), Optional, Intent (InOut) :: ra (nstsv)
    Complex (8), Optional, Intent (InOut) :: ca (nstsv)
  ! local variables
    Real (8), Allocatable :: evalsv_t (:)
    Real (8), Allocatable :: ra_t (:)
    Complex (8), Allocatable :: ca_t (:)
    Integer,  Allocatable :: idx (:)
    Integer :: ist

    Allocate (evalsv_t(nstsv))
    If (present(ra)) Allocate (ra_t(nstsv))
    If (present(ca)) Allocate (ca_t(nstsv))
    Allocate (idx(nstsv))

    Call getevalsv (vpl, evalsv_t)

    Call sortidx (nstsv, evalsv_t, idx)

    If (present(ra)) Then
        ra_t = ra
        Do ist=1, nstsv
            ra(ist) = ra_t(idx(ist))
        End Do
    Else If (present(ca)) Then
        ca_t = ca
        Do ist=1, nstsv
            ca(ist) = ca_t(idx(ist))
        End Do
    End If

    Deallocate(evalsv_t, idx)
    If (present(ra)) Deallocate (ra_t)
    If (present(ca)) Deallocate (ca_t)
End Subroutine
!EOC
End Module m_sorteval
