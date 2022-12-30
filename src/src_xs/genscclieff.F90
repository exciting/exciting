! Copyright (C) 2008-2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genscclieff(iqr,  nmax, n, scieff)
  use modinput, only: input
  use constants, only: zzero
  use putgeteps0, only: geteps0_finite_q, geteps0_zero_q
  use modxs, only: eps0dirname, vqcr
  use m_genfilname
  use math_utils, only: all_zero
  use precision, only: dp

  implicit none

  ! Arguments
  integer, intent(in) :: iqr, n, nmax
  complex(8), intent(out) :: scieff(nmax, nmax)

  ! Local variables
  logical :: tq0
  complex(8), allocatable :: scrn(:, :), scrnw(:, :, :), scrnh(:, :)
  character(256) :: fneps0



  allocate(scrn(n, n), scrnw(n, 2, 3), scrnh(3, 3))
  scrn=zzero
  scrnw=zzero
  scrnh=zzero

  ! Check whether q=0, that is it checks whether mod_qpoint::vqc(:,iqr) has
  ! a norm smaller than epslat (mod_qpoint should contain the non reduced q grid)
  tq0 = all_zero(vqcr(:,iqr), tol=1e-12_dp)

  ! Read the Coulomb-symmetrized macroscopic dielectric function/tensor
  ! in RPA for requested q point and zero frequency. The Cartesian components
  ! of the head (G=G'=q=0) are also symmetrized w.r.t. the lattice symmetry.
  call genfilname(basename=trim(adjustl(eps0dirname))//'/'//'EPS0',&
   & appfilext=.true., iq=iqr, filnam=fneps0)

  ! Calculate effective screened interaction
  if(tq0) then
    ! Read form direct access file.
    call geteps0_zero_q(qvec=vqcr(:,iqr), iq=iqr, iw=1, w=0.0d0,&
          eps0=scrn, eps0wg=scrnw, eps0hd=scrnh, fname=fneps0,&
          debug=input%xs%dbglev>2)
    ! Averaging using Lebedev-Laikov spherical grids
    call angavsc0(n, nmax, scrnh, scrnw, scrn, scieff)
  else
    ! Read form direct access file.
    call geteps0_finite_q(qvec=vqcr(:,iqr), iq=iqr, iw=1, w=0.0d0,&
                eps0=scrn, fname=fneps0, debug=input%xs%dbglev>2)
    ! Averaging using numerical method and extrapolation
    call avscq(iqr, n, nmax, scrn, scieff)
  end if

  deallocate(scrn, scrnw, scrnh)

end subroutine genscclieff
