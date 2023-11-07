! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: borncharges
! !INTERFACE:
!
!
subroutine borncharges( pol, disp, zstar, sumrule, mode)
  ! !USES:
  use modmpi
  use modinput
  use mod_lattice
  use mod_atoms

  ! !INPUT/OUTPUT PARAMETERS:
  !   pol     : polarization vectors for two different displacements of each atom in each direction (in,real(3,2,3,natmtot))
  !   disp    : displacement difference between both atom positions (in,real)
  !   zstar   : Born effective charge tensor for each atom (inout,real(3,3,0:natmtot))
  !   sumrule : if .true., the acoustic sum rule is imposed on the charge tensors (in,logical)
  !   mode    : either 'calc', 'write' or 'read' (in,character(*))
  ! !DESCRIPTION
  !   Computes the Born effective charge tensor from finite differences of polarization vectors
  !   $$ Z^{\kappa}_{\alpha\beta} = \frac{\Delta P^{\kappa}_{\alpha\beta}}{d} \;,$$
  !   where $\Delta P^{\kappa}_{:\beta}$ is the difference in the polarization with atom $\kappa$
  !   being at two different positions in direction $\beta$. The difference vector is mapped to the
  !   shortest possible one by subtracting an integer number of polarization quanta. The distance
  !   of both positions of atom $\kappa$ is given by $d$. \\
  !   If {\tt sumrule = .true.}, the acousting sum rule, i.e. $\sum_{\kappa} Z^{\kappa}_{\alpha\beta} = 0$, 
  !   is imposed on the result my deviding the actual deviation from 0 by the number of atoms and 
  !   subtracting it uniformly from all charge tensors. In this case, the portion that was subtracted
  !   from the original result is returned in {\tt zstar(:,:,0)}.\\
  !   If {\tt mode='write'}, {\tt zstar} is interpreted as input and written to the file {\tt ZSTAR.OUT}.\\
  !   If {\tt mode='read'}, the charge tensors are read from file and returned into {\tt zstar}.\\
  !   Otherwise, the charge tensors are calculated and returned into {\tt zstar}.
  !
  ! !REVISION HISTORY:
  !   Created July 2020 (SeTi)
!EOP
!BOC
  implicit none

  real(8), intent( in)      :: pol(3,2,3,natmtot), disp
  logical, intent( in)      :: sumrule
  character(*), intent( in) :: mode
  real(8), intent( inout)   :: zstar(3,3,0:natmtot)

  integer :: ia, is, ias, ip, vi(3), un
  real(8) :: vr(3)

  if( trim( mode) == 'write') then
    call put( zstar)
    return
  elseif( trim( mode) == 'read') then
    call get( zstar)
    return
  end if

  zstar = 0.d0
  do is = 1, nspecies
    do ia = 1, natoms( is)
      ias = idxas( ia, is)
      do ip = 1, 3
        call r3mv( ainv, pol(:,2,ip,ias)-pol(:,1,ip,ias), vr)
        vr = vr*omega
        call r3ws( 1.d-16, input%structure%crystal%basevect, vr, vi)
        call r3mv( input%structure%crystal%basevect, vr, zstar(ip,:,ias))
      end do
      zstar(:,:,ias) = zstar(:,:,ias)/disp
      zstar(:,:,0) = zstar(:,:,0) + zstar(:,:,ias)
    end do
  end do
  zstar(:,:,0) = zstar(:,:,0)/natmtot

  if( sumrule) then
    do ias = 1, natmtot
      zstar(:,:,ias) = zstar(:,:,ias) - zstar(:,:,0)
    end do
  end if
  
  return

  contains
    subroutine put( z)
      use mod_misc, only: filext
      use phonons_io_util, only: ph_io_write_borncharge
      real(8), intent( in) :: z(3,3,0:natmtot)

      logical :: success

      call ph_io_write_borncharge( z(:, :, 1:natmtot), 'ZSTAR'//trim( filext ), success, &
        sumrule_correction=z(:, :, 0) )
      call barrier
    end subroutine

    subroutine get( z)
      use mod_misc, only: filext
      use phonons_io_util, only: ph_io_read_borncharge
      real(8), intent( out) :: z(3,3,0:natmtot)
      
      logical :: success
      real(8) :: corr(3, 3)
      real(8), allocatable :: borncharge(:, :, :)

      call ph_io_read_borncharge( borncharge, 'ZSTAR'//trim( filext ), success, &
        sumrule_correction=corr )
      z(:, :, 0) = corr
      z(:, :, 1:natmtot) = borncharge

      if( .not. sumrule ) then
        do ias = 1, natmtot
          z(:,:,ias) = z(:,:,ias) + z(:,:,0)
        end do
      end if
    end subroutine
end subroutine borncharges
!EOC
