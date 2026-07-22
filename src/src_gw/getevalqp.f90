
!> Apply GW quasi-particle corrections to a set of Kohn-Sham eigenvalues.
!>
!> Reads the QP and KS eigenvalues stored by [[putevalqp]] (typically
!> `EVALQP.OUT`), forms the QP correction
!> $$ \Delta_{n\mathbf{k}} = E^\mathrm{QP}_{n\mathbf{k}} - \epsilon^\mathrm{KS}_{n\mathbf{k}} $$
!> on the k-grid of the file, Fourier-interpolates it onto the current k-grid
!> and adds it to `eqp2`, which on exit holds the QP energies measured relative
!> to the QP Fermi energy `eferqp`.
!>
!> @warning `eqp2` is `intent(inout)` and this routine is **additive**: the
!> interpolated correction is added to whatever is passed in. The caller must
!> therefore supply the *Kohn-Sham* eigenvalues of the target grid, and must do
!> so exactly once. Calling this routine on an array that already holds QP
!> energies applies the correction a second time and silently double-shifts the
!> spectrum. @endwarning
!>
!> Only bands in the range `[ibgw, nbgw]` read from the file are touched (capped
!> at `nstsv`); eigenvalues outside that window are left unchanged.
!>
!> Two code paths with *different* semantics exist:
!>
!> * `nkp1 > 1` — the general case described above, additive.
!> * `nkp1 == 1` (molecules, single k-point) — no interpolation is possible, so
!>   the QP energies are assigned *absolutely*, `eqp2 = eqp1`, and `eferqp` is
!>   not subtracted. This branch is idempotent, the general one is not.
!>   It requires `nkp2 == 1` and aborts otherwise.
!>
!> @note The interpolation target is `mod_kpoint::vkl`, not the `kvecs2` dummy
!> argument, and the output loop runs over `mod_kpoint::nkpt` rather than
!> `nkp2`. Callers must make sure that the current `vkl`/`nkpt` describe the grid
!> `eqp2` is dimensioned for; this is why BSE on top of GW is restricted to zero
!> momentum transfer, where the shifted and unshifted grids coincide. @endnote
!>
!> Side effects: sets the `mod_bands` grid variables and the `modgw` band limits
!> and Fermi energies (`ibgw`, `nbgw`, `eferqp`, `eferks`) from the file
!> contents. Via [[fourintp]] it may also reset `mod_symmetry::nsymcrys` to 1
!> (when `gw/@symmetryBandstructure` is false) without restoring it, so callers
!> that rely on the crystal symmetry afterwards have to save and restore it.
subroutine getevalqp(fname, nkp2, kvecs2, eqp2)
  use bandstructure,            only: fourintp
  use constants,                only: zzero
  use mod_bands,                only: nkp1, kvecs1, eks1, eqp1
  use mod_eigenvalue_occupancy, only: nstsv
  use mod_kpoint,               only: nkpt, vkl
  use mod_large_io,             only: inquire_large, open_direct_unformatted_large
  use modgw,                    only: ibgw, nbgw, eferqp, eferks
  use modmpi,                   only: terminate_if_false
  use precision,                only: dp, i32, long_int

  implicit none

  !> Name of the file holding the QP eigenvalues, e.g. `EVALQP.OUT`.
  character(*), intent(in)    :: fname
  !> Number of k-points `eqp2` is dimensioned for.
  integer(i32), intent(in)    :: nkp2
  !> k-points of the target grid, in lattice coordinates.
  !> Currently unused: the interpolation is carried out on `mod_kpoint::vkl`.
  real(dp),     intent(in)    :: kvecs2(3, nkp2)
  !> On entry the KS eigenvalues of the target grid, on exit the corresponding
  !> QP energies. See the warning above: the correction is *added*.
  real(dp),     intent(inout) :: eqp2(nstsv, nkp2)

  logical :: exist
  integer(i32) :: ik, ib, nb, nk, nqp, unit
  integer(long_int) :: recl
  real(dp), allocatable :: eqp(:)
  complex(dp), allocatable :: de1(:,:), de2(:,:)

  !-----------------------------------------------------------------------------
  ! Read the file
  !-----------------------------------------------------------------------------      
  inquire(File=fname, Exist=exist)
  call terminate_if_false( exist, 'ERROR(getevalqp): File ' // trim( fname ) // ' does not exist!')
      
  call inquire_large( recl, [nkp1, ibgw, nbgw] )
  call open_direct_unformatted_large( unit, trim(fname), "read", recl, "old" )
  read(unit, Rec=1) nkp1, ibgw, nbgw
  close(unit)
      
  allocate(kvecs1(1:3,nkp1))
  allocate(eqp1(ibgw:nbgw,nkp1))
  allocate(eks1(ibgw:nbgw,nkp1))
  
  call inquire_large( recl, [nkp1, ibgw, nbgw], kvecs1(1:3,1), &
                      eqp1(ibgw:nbgw,1), eks1(ibgw:nbgw,1), [eferqp, eferks] )
  
  call open_direct_unformatted_large( unit, trim(fname), "read", recl, "old" )
  
  nqp = nbgw-ibgw+1
  allocate(eqp(nqp))

  do ik = 1, nkp1
    read(unit, Rec=ik) nk, ib, nb, kvecs1(:,ik), &
         eqp1(ibgw:nbgw,ik), eks1(ibgw:nbgw,ik), &
         eferqp, eferks
  end do
  close(unit)

  !----------------------------------------------
  ! Special case of only one k-point (molecules)
  ! NOTE: unlike the interpolated branch below, this overwrites eqp2 instead of
  !       adding to it, and does not subtract eferqp.
  !----------------------------------------------
  if (nkp1==1) then
    if (nkp2==1) then
      do ib = ibgw, min(nbgw, nstsv)
        eqp2(ib,1) = eqp1(ib,1)
      end do
      deallocate(kvecs1, eqp1, eks1)
      return
    else
      write(*,*) 'ERROR(getevalqp):' 
      write(*,*) '  Interpolation is not possible!'
      write(*,*) '  EVALQP.OUT file contains data only for a single k-point.'
      stop
    end if
  end if 

  !-----------------------------------------------------------------------------
  ! Interpolate the energies
  ! The QP correction eqp1 - eks1, not the QP energy itself, is what gets
  ! interpolated and then added to the incoming KS eigenvalues.
  !-----------------------------------------------------------------------------
  allocate(de1(nkp1,ibgw:nbgw))
  do ik = 1, nkp1
    de1(ik,:) = cmplx(eqp1(ibgw:nbgw,ik)-eks1(ibgw:nbgw,ik), 0.d0, 8)
  enddo

  allocate(de2(nkpt,ibgw:nbgw))
  de2(:,:) = zzero
  
  call fourintp(de1, nkp1, kvecs1, de2, nkp2, vkl, nbgw-ibgw+1)

  do ib = ibgw, min(nbgw,nstsv)
     do ik = 1, nkpt
        eqp2(ib,ik) = eqp2(ib,ik) + dble(de2(ik,ib)) - eferqp
     enddo 
  enddo

  deallocate(kvecs1, eqp1, eks1)

end subroutine
