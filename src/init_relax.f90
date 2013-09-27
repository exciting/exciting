!BOP
! !ROUTINE: init_relax
! !INTERFACE:
!
!
Subroutine init_relax
! !USES:
      Use modinput
      Use modmain
!EOP
!BOC
      Implicit None
! local variables
      Integer :: ik, ispn

!____________________________
! lattice and symmetry set up

      Call findsymcrys                      ! find the crystal symmetries and shift atomic positions if required

      Call findsymsite                      ! find the site symmetries

      Call checkmt                          ! check for overlapping muffin-tins

      Call gencfun                          ! generate the characteristic function

      Call energynn                         ! determine the nuclear-nuclear energy

!_________________________________________________
! generate structure factors for G and G+k-vectors

      Call gensfacgp (ngvec, vgc, ngvec, sfacg)

      Do ik = 1, nkpt
          Do ispn = 1, nspnfv
              Call gensfacgp (ngk(ispn, ik), vgkc(:, :, ispn, ik), ngkmax, sfacgk(:, :, ispn, ik))
          End Do
      End Do

      Return
End Subroutine
!EOC
