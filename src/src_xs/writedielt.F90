! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!> Routine to write the dielectric tensor to file
!> Theory level can be independent particles (IPA) or 
!> random phase (RPA) approximation
subroutine writedielt(file_tag, n_freq, freq, diel_tens, theory_level)
  use mod_misc, only: filext
  use m_getunit
  use modinput, only: input
  use precision, only: sp,dp

  implicit none

  !> Filetag
  character(*), intent(in) :: file_tag
  !> Number of frequencies
  integer(sp), intent(in) :: n_freq
  !> Frequencies
  real(dp), intent(in) :: freq(n_freq)
  !> Dielectric tensor
  complex(dp), intent(in) :: diel_tens(3, 3, n_freq)
  !> Level of theory: IPA (theory_level=0) or RPA (theory_level=0)
  integer, intent(in) :: theory_level

  ! local variables
  !> File unit
  integer(sp) :: un
  !> Running index rows
  integer(sp) :: i_rows
  !> Running index frequencies
  integer(sp) :: i_freq

  call getunit(un)
  open(un, file=trim(file_tag)//trim(filext), form='formatted',&
    & action='write', status='replace')
  write(un,*)

  if (theory_level .eq. 0)  then    
    write(un, '(" (dielectric tensor, independent particle approximation)")')    
  else
    write(un, '(" (dielectric tensor, random phase approximation )")')
  end if
  write(un,*)

  ! Loop over frequencies
  do i_freq = 1, n_freq
    write(un, '(" frequency index and value: ",i6,f14.8)') i_freq, freq(i_freq)
    write(un, '(" real part, imaginary part below")')
    write(un, '(3f14.8,5x,3f14.8)')&
      & (dble(diel_tens(i_rows, :, i_freq)), aimag(diel_tens(i_rows, :, i_freq)), i_rows=1, 3)
    write(un,*)
  end do

  close(un)

end subroutine writedielt
