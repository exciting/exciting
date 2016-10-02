! Copyright (C) 2009 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: putscreen
! !INTERFACE:
subroutine putscreen(un, tq0, n, chi0, chi0h, chi0w)
! !USES:
  use mod_constants, only: krondelta
! !INPUT/OUTPUT PARAMETERS:
! In:
! integer :: un  ! Unit to wirte to 
! logical :: tp0 ! Flag if iq is the q=0 q-point
! integer :: n   ! Number of G+q vectros
! complex(8) :: ch0(n,n)    ! Body of RPA density-density response matrix 
! complex(8) :: ch0h(3,3)   ! Body of RPA density-density response matrix
! complex(8) :: ch0w(n,2,3) ! Wings of RPA density-density response matrix 
!
! !DESCRIPTION:
!   Writes the microscopic V-symmetrized Kohn-Sham dielectric function/tensor
!   $\tilde{\epsilon}^0_{\bf{GG'}}({\bf q},\omega) = \delta_{\bf{GG'}}
!   - \tilde{\chi}^0_{\bf{GG'}}({\bf q},\omega)$ for fixed frequency and q point
!   to a human readable text file.
!   Is used only for $\omega = 0$.
!
! !REVISION HISTORY:
!   Added to documentation scheme. (Aurich)
!EOP
!BOC

  implicit none

  ! Input parameters
  logical, intent(in) :: tq0
  integer, intent(in) :: un, n
  complex(8), intent(in) :: chi0(n, n), chi0h(3, 3), chi0w(n, 2, 3)

  ! Local variables
  integer :: ig1, ig2, i, j
  real(8) :: r1

  ! Loop over G+q and G'+q indices
  do ig1 = 1, n
    do ig2 = 1, n

      r1 = 0.d0
      if(ig1 .eq. ig2) r1 = 1.d0

      if(tq0) then

        ! Write head
        if((ig1 .eq. 1) .and. (ig2 .eq. 1)) then
          ! Write Cartesian directions with negative sign.
          ! Write real part, imaginary part and modulus of 1-chi
          write(un, '(i8,1x,i8,1x,E23.16,1x,E23.16,1x,E23.16)')&
            & ( (-i, -j, dble(krondelta(i, j))-chi0h(i, j),&
            & abs(dble(krondelta(i, j))-chi0h(i, j)), j=1, 3), i=1, 3)
        end if

        ! Write wings
        if((ig1 .eq. 1) .and. (ig2 .ne. 1)) then
          write(un, '(i8,1x,i8,1x,E23.16,1x,E23.16,1x,E23.16)')&
            & ( -i, ig2, -chi0w(ig2, 1, i), abs(-chi0w(ig2, 1, i)), i=1, 3)
        end if
        if((ig1 .ne. 1) .and. (ig2 .eq. 1)) then
          write(un, '(i8,1x,i8,1x,E23.16,1x,E23.16,1x,E23.16)')&
            & (ig1, -j, -chi0w(ig1, 2, j), abs(-chi0w(ig1, 2, j)), j=1, 3)
        end if

        ! Write body
        if((ig1 .ne. 1) .and. (ig2 .ne. 1)) then
          write(un, '(i8,1x,i8,1x,E23.16,1x,E23.16,1x,E23.16)')&
            & ig1, ig2, r1 - chi0(ig1, ig2), abs(r1-chi0(ig1, ig2))
        end if

      else

        ! Write full
        write(un, '(i8,1x,i8,1x,E23.16,1x,E23.16,1x,E23.16)')&
          & ig1, ig2, r1 - chi0(ig1, ig2), abs(r1-chi0(ig1, ig2))

      end if

    end do
  end do

  write(un,'("#",1x,a6,1x,a8,1x,a23,1x,a23,1x,a23)')&
    & "i1", "i2", "Re(esp_i1i2)", "Im(eps_i1i2)", "Abs(eps_i1i2)"
  write(un,'("#",1x,a)') "RPA static screening (microscopic dielectric tensor/function)"
  write(un,'("#",1x,a)') "frequency w = 0"
  if(tq0) then
    write(un,'("#",1x,a)') "q=0 : Head and wings of epsilon_GG'(q=0,w=0) present."
  end if
  write(un,'("#",1x,a)') "Positive integer index G+q vectors."
  write(un,'("#",1x,a)') "Negative integer index cartesian directrions for head and wings."
  write(un,'("#",1x,a)') "(Head has 2 catesian indices and no G+q index,"
  write(un,'("#",1x,a)') "while the wings have one cartesian and one G+q index)."


end subroutine
!EOC
