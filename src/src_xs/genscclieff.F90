! Copyright (C) 2008-2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genscclieff(iqr, iqrnr, nmax, n, scieff)
  use modinput, only: input
  use mod_constants, only: zzero
  use m_putgeteps0
  use modxs, only: eps0dirname
  use m_genfilname

  implicit none

  ! Arguments
  integer, intent(in) :: iqr, iqrnr, n, nmax
  complex(8), intent(out) :: scieff(nmax, nmax)

  ! Local variables
  logical :: tq0
  complex(8), allocatable :: scrn(:, :), scrnw(:, :, :), scrnh(:, :)
  character(256) :: fneps0

  ! External functions
  logical, external :: tqgamma

  allocate(scrn(n, n), scrnw(n, 2, 3), scrnh(3, 3))
  scrn=zzero
  scrnw=zzero
  scrnh=zzero

  ! Check whether q=0, that is it checks whether mod_qpoint::vqc(:,iqr) has
  ! a norm smaller than epslat (mod_qpoint should contain the non reduced q grid)
  tq0 = tqgamma(iqrnr)
  
  ! Read the Coulomb-symmetrized macroscopic dielectric function/tensor 
  ! in RPA for requested q point and zero frequency. The Cartesian components
  ! of the head (G=G'=q=0) are also symmetrized w.r.t. the lattice symmetry.
  if(input%xs%bse%beyond) then
    call genfilname(basename=trim(adjustl(eps0dirname))//'/'//'EPS0',&
     & appfilext=.true., iq=iqr, filnam=fneps0)
    ! Read form direct access file.
    call geteps0(reduced=.true., iq=iqr, iw=1, w=0.0d0,&
      & eps0=scrn, eps0wg=scrnw, eps0hd=scrnh, fname=fneps0)
  else
    ! Read from text file 
    call getscreen(iqr, n, scrnh, scrnw, scrn)
  end if

  ! Calculate effective screened interaction
  if(tq0) then
    ! Averaging using Lebedev-Laikov spherical grids
    call angavsc0(n, nmax, scrnh, scrnw, scrn, scieff)
  else
    ! Averaging using numerical method and extrapolation
    call avscq(iqr, n, nmax, scrn, scieff)
  end if

  deallocate(scrn, scrnw, scrnh)

end subroutine genscclieff
