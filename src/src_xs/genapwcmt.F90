! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
subroutine genapwcmt(lmax, ngp, isti, istf, apwalm, evecfv, wfcmt)

  use modmain
  use modinput
  implicit none

  ! arguments
  integer, intent(in) :: lmax, isti, istf
  integer, intent(in) :: ngp
  complex(8), intent(in) :: apwalm(ngkmax, apwordmax, lmmaxapw, natmtot)
  complex(8), intent(in) :: evecfv(nmatmax, nstfv, nspnfv)
  complex(8), intent(out) :: wfcmt(istf-isti+1, apwordmax, (lmax+1)**2, natmtot)

  ! local variables
  integer :: ist, istc, is, ia, ias

  if(lmax .gt. input%groundstate%lmaxapw) then
    write(*,*)
    write(*, '("Error(genapwcmt): lmax > lmaxapw : ", i8)') lmax
    write(*,*)
    stop
  end if

  do istc = isti, istf
    ist = istc - isti + 1
    do is = 1, nspecies
      do ia = 1, natoms(is)
        ias = idxas(ia, is)
        call genapwcmt_part(lmax, ngp, ia, is, apwalm, evecfv(:, istc, 1),&
          & wfcmt(ist, :, :, ias))
      end do
    end do
  end do

  contains

    subroutine genapwcmt_part(lmax, ngp, ia, is, apwalm, evecfv, fcmt)

      implicit none

      ! arguments
      integer, intent(in) :: lmax, ia, is
      integer, intent(in) :: ngp
      complex(8), intent(in) :: apwalm(:, :, :, :)
      complex(8), intent(in) :: evecfv(:)
      complex(8), intent(out) :: fcmt(:, :)

      ! external functions
      complex(8) zdotu
      external zdotu

      ! local variables
      integer :: ias, l, m, lm, io

      ias = idxas(ia, is)
      do l = 0, lmax
        do m = - l, l
          lm = idxlm(l, m)
          do io = 1, apword(l, is)
             fcmt(io, lm) = zdotu(ngp, evecfv, 1, apwalm(1, io, lm, ias), 1)
          end do
        end do
      end do

    end subroutine genapwcmt_part

end subroutine genapwcmt
