! Copyright (C) 2014 exciting team 
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
! History:
! created July 2014 by Stefan Kontur
!
subroutine check_raman (imode, evec_phon, irep, active)
!
use mod_atoms
use mod_symmetry
use raman_symmetry
implicit none
! arguments
integer, intent(in) :: imode
real(8), intent(in) :: evec_phon (3*natmtot)
integer, intent(out) :: irep
logical, intent(out) :: active
! local variables
real(8) :: evec_phon_norm(3*natmtot), evec_phon_sop(3*natmtot)
real (8) :: rotvec (3)
integer :: ia, ja, i, j, k
!integer :: numsop
real (8) :: norm_factor
logical :: lassign
!
open( unit=13, file='RAMAN_SYM.OUT', status='old', action='write', position='append')
!
active = .true.
!
! check given mode
!
write(13, '(/," Phonon eigenvector of mode ",i3," :",/, &
         &    " --------------------------------",/)') imode
do ia = 1,natmtot
 write(13, '(" Atom ",i2," : ",3f10.4)') ia,evec_phon( (3*(ia-1)+1):(3*ia) )
enddo
write(13,*)
!
! normalize phonon eigenvector
norm_factor = 0.d0
do i = 1,3
   do ia = 1,natmtot
      norm_factor = norm_factor + evec_phon(3*(ia-1)+i)**2
   enddo
enddo
if (norm_factor .gt. eps) then
   evec_phon_norm = evec_phon / sqrt(norm_factor)
else
   evec_phon_norm = evec_phon
endif
!
lassign = .false.
!
! apply projection operators, using characters
do k = 1,cl
   evec_phon_sop = 0.d0
   do ia = 1,natmtot
      ! loop over SOPs
      do j = 1, numsop
        ja = atom_sop(ia, j)
!       call r3mv (symlatc(:, :, lsplsymc(j)), evec_phon((3*(ia-1)+1):(3*ia)), rotvec(:))
        call r3mv (sopmatc(:, :, j), evec_phon((3*(ia-1)+1):(3*ia)), rotvec(:))
        evec_phon_sop((3*(ja-1)+1):(3*ja)) = evec_phon_sop((3*(ja-1)+1):(3*ja)) &
 &       + dble(charact(class(j),k))*rotvec(:)
      enddo
   enddo
   norm_factor = 0.d0
   do i = 1,3
    do ia = 1,natmtot
      norm_factor = norm_factor + evec_phon_sop(3*(ia-1)+i)**2
    enddo
   enddo
   if (norm_factor .gt. eps) evec_phon_sop = evec_phon_sop / sqrt(norm_factor)
   if (all(abs(evec_phon_sop(:)-evec_phon_norm(:)) .lt. eps) .or. &
 &     all(abs(evec_phon_sop(:)+evec_phon_norm(:)) .lt. eps)) then
      ! the projected vector equals the original one
      lassign = .true.
      active = raman_active(k)
      irep = k
      write(13, '(" Phonon mode belongs to IREP   : ",i2," (",a4,")")') k, irep_ch(k)
      write(13, '(" Raman active                  : ",l2)') raman_active(k)
      write(13, '(" Triggering Raman computation... ",/)')
      exit
   endif
enddo
if (.not. lassign) write(13,*) ' No symmetry found for given phonon mode!'
!
close (13)
!
!
!
!
!220  format(' Atom ',i2,' : ',3f10.4)
!228  format(' No. degenerate modes combined : ',i2)
!230  format(' Phonon mode belongs to IREP   : ',i2,' (',a4,')')
!240  format(' Raman active                  : ',l2)
!
!
end subroutine check_raman
