!
!BOP
! !ROUTINE: putvnlmat
! !INTERFACE:
!
subroutine putvnlmat()
! !USES:
  use modmain
  use mod_hybrids
  use modmpi
!
! !DESCRIPTION:
!   Writes the APW matrix elements of the non-local potential
!   into the file VNLMAT.OUT
!
! !REVISION HISTORY:
!   Created March 2015 (UW)
!   Changed May 2016 (DIN)
!
!EOP
!BOC

  implicit none
  integer(8) :: recl
  integer :: ik, ikfirst, iklast

!$OMP CRITICAL

  ikfirst = firstk(rank)
  iklast = lastk(rank)

  ! Save < APW' | \Sigma_x | APW >
  inquire(IoLength=Recl) nkpt, nmatmax ,vnlmat(:,:,ikfirst)
  open(70, File='VNLMAT.OUT', Action='WRITE', Form='UNFORMATTED', &
  &    Access='DIRECT', status='REPLACE', Recl=Recl)
  do ik = 1, nkpt
    ! check which rank should print
    if ((ik >= ikfirst).and.(ik <= iklast)) then
      write(70, Rec=ik) nkpt, nmatmax ,vnlmat(:,:,ik)
    end if
    call barrier
  end do ! ik
  close(70)

!$OMP END CRITICAL
  return
end subroutine
