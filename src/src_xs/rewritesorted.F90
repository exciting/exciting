!BOP
! !ROUTINE: rewritesorted
! !INTERFACE:
Subroutine rewritesorted
! !USES:
      Use m_sorteval
      Use modmain
      Use modxs
! !DESCRIPTION:
!   In a collinear spin calculation, the second variational hamiltonian is
!   block diagonalized in subroutine seceqnsv. This leads to block-sorted
!   eigenvalues, eigenvectors and occupations files into increasing
!   energy for each spin chanel. The XS part needs fully sorted arrays
!   intro incresing energy in order to correctly find bandlimits.
!   Accordingly, this subroutine rewrites the corresponding files
!   sorted into increasing energy when the calculation is collinear.
!
! !REVISION HISTORY:
!   Created July 2013 SR
!EOP
!BOC
      Implicit None
  ! local variables
      Integer :: ik, ist
      Logical, External :: iscollinear
      If ( iscollinear() ) Then
        Allocate (evecsv(nstsv, nstsv))
        Do ik = 1, nkpt
     ! read files
            Call getevecsv (vkl(1, ik), evecsv)
            Call getevalsv (vkl(1, ik), evalsv(1, ik))
            Call getoccsv (vkl(1, ik), occsv(1, ik))
     ! sort arrays
            Call sorteval (vkl(1, ik), ra=evalsv(1:nstsv, ik))
            Do ist=1, nstsv
                Call sorteval (vkl(1, ik), ca=evecsv(ist,1:nstsv))
            End Do
            Call sorteval (vkl(1, ik), ra=occsv(1:nstsv, ik))
     ! write files
            Call putevecsv (ik, evecsv)
            Call putevalsv (ik, evalsv(1, ik))
            Call putoccsv (ik, occsv(1, ik))
        End Do
  ! read files
  ! write files
        Call writeeval
        Deallocate (evecsv)
      End If
End Subroutine rewritesorted
!EOC
