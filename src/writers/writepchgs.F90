
! Copyright (C) 2010 S. Sagmeister and  C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: writepchgs
! !INTERFACE:
subroutine writepchgs(fnum,lmax)
! !USES:
      use modinput
      use modmain
! !DESCRIPTION:
!  Write partial charges to file.
!
! !REVISION HISTORY:
!   Created 2010 (Sagmeister)
!EOP
!BOC
  implicit none
  ! arguments
  integer, intent(in) :: fnum, lmax
  ! local variables
  integer :: ist,is,ia,ias,l,m,lm
  real(8) :: t1
  write(fnum,*)
  write(fnum,'("iteration number: ",i6)') iscl
  Write (fnum, '(I6, " : nstsv")') nstsv
  Write (fnum, '(f12.6, " : muffin-tin charge")') sum(chgpart)
  do ist=1,nstsv
    t1 = maxval(wkpt(:nkpt) * occsv(ist,:))
    if (abs(t1) .gt. input%groundstate%epsocc) Then
      Write (fnum, '(I6, " : state")') ist
      write(fnum,'("  sum over atoms and lm : ",f12.6)') sum(chgpart(:,:,ist))
      do is=1,nspecies
        do ia=1,natoms(is)
          ias=idxas(ia,is)
          Write (fnum, '("  Species : ", I4, " (", A, "), atom : ", I4)') &
            & is, trim &
            & (input%structure%speciesarray(is)%species%chemicalSymbol), &
            & ia
          write(fnum,'("  sum over lm : ",f12.6)') sum(chgpart(:,ias,ist))
          do l=0,lmax
            Write (fnum, '("   l-value : ", I4, " , sum over m : ",f12.6," ; m-components below")') &
              l, sum(chgpart(idxlm(l,-l):idxlm(l,l),ias,ist))
            write(fnum,'(100f12.6)') chgpart(idxlm(l,-l):idxlm(l,l),ias,ist)
          end do
        end do
      end do
    end if
  end do
  deallocate(chgpart)
end subroutine
!EOC
