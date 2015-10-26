!
!BOP
! !ROUTINE: putvnlmat
! !INTERFACE:
!
Subroutine putvnlmat 
! !USES:
      Use modmain
      use mod_hybrids
!
! !DESCRIPTION:
!   Writes the APW matrix elements of the non-local potential
!   into the file VNLMAT.OUT
!
! !REVISION HISTORY:
!   Created March 2015 (UW)
!
!EOP
!BOC

      Implicit None
  ! local variables
      Integer(8) :: recl,ik
  !$OMP CRITICAL
    Inquire (IoLength=Recl) nkpt, nmatmax ,vnlmat(:,:,1)
    Open (70, File='VNLMAT.OUT', Action='WRITE', Form='UNFORMATTED', &
   &   Access='DIRECT', status='REPLACE', Recl=Recl)
    do ik = 1, nkpt
        write(70, Rec=ik) nkpt, nmatmax ,vnlmat(:,:,ik)
    end do ! ik
    Close(70)
 !$OMP END CRITICAL
      Return
End Subroutine putvnlmat
