!BOP
! !ROUTINE: writeqmtpts
! !INTERFACE:
!
subroutine writeqmtpts
! !USES:
  use modxs, only: nqmt, vqlmt, totalqcmt, vqcmt, ivgmt, vgcmt
  use m_getunit
  use m_genfilname
! !DESCRIPTION:
!   Writes the momentum transfer ${\bf Q}$-points
!
! !REVISION HISTORY:
!   Created 2017 (BA)
!EOP
!BOC

  implicit none

  ! Local variables
  integer :: iqmt, un
  character(256) :: filnam

  call getunit(un)
  Call genfilname(basename='QMTPOINTS', appfilext=.True., filnam=filnam)

  open(un, file=trim(filnam), action='WRITE', form='FORMATTED')
  write(un, '("# Momentum transfer vectors from qpointlist element.")')
  write(un, '("# Momentum transfer Qmt = Gmt + qmt, where qmt is in the unit cell.")')
  write(un, '("#")')
  write(un, '("#",1x,a6,3a8,10a18)')&
    & "iqmt", "ivgmt_1", "ivgmt_2", "ivgmt_3", "vqlmt_1", "vqlmt_2", "vqlmt_3",&
    & "vgcmt_1", "vgcmt_2", "vgcmt_3", "vqcmt_1", "vqcmt_2", "vqcmt_3", "|Q|"

  do iqmt = 1, nqmt
    write(un, '(4i8, 10E18.10)') iqmt, ivgmt(:,iqmt), vqlmt(:,iqmt),&
      & vgcmt(:,iqmt), vqcmt(:,iqmt), Norm2(totalqcmt(:,iqmt))
  end do

  close(un)

end subroutine writeqmtpts
!EOC
