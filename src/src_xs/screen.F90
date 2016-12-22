! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: screen
! !INTERFACE:
subroutine screen
! !USES:
  use modmpi
!  use modxs, only: nwdf, unitout
  use m_genfilname
  use m_filedel

use mod_ematgrids
use m_ematqk
use mod_Gkvector, only: gkmax
use mod_eigensystem, only: nmatmax
use mod_eigenvalue_occupancy, only: nstfv
use mod_spin, only: nspnfv
use modxs
use m_b_ematqk
use mod_constants
use m_xsgauntgen
use m_findgntn0
use modinput
use mod_APW_LO
use m_writecmplxparts

use m_ematqk_orig
! !DESCRIPTION:
! This is a wrapper subroutine that backs up the number of 
! $ \omega $ points and calls the {\tt df} subroutine. Afterwards
! it restores them. It also sets the output file extensions for
! screening related quantities.
! 
! !REVISION HISTORY:
!   Added to documentation scheme. (Aurich)
!EOP
!BOC

  implicit none

  ! local variables
  integer :: nwdft


!  !!+++++++TEST HACK ++++++++++++!!
!  real(8) :: vqmtl(3)
!  complex(8), allocatable :: eveck(:,:,:), eveckq(:,:,:), ematsepp(:,:,:), emat(:,:,:), ematsag(:,:,:)
!  complex(8), allocatable :: evecks(:,:,:,:)
!  integer(4) :: ik, ikq, iq, myngq, igq
!  real(8) :: t0, t1, maxdiff
!  type(bcbs) :: bc
!
!  write(*,*) "Screen here"
!  call init0
!  call init1
!  call init2
!
!
!  !write(*,*) "Making grids"
!  vqmtl = 0.0d0
!  call ematgrids_init(vqmtl, gkmax)
!
!  call genfilname(dotext='_SCR.OUT', setfilext=.true.)
!  write(*,*) "Getting all evecs"
!  allocate(evecks(nmatmax, nstfv, nspnfv, kkpmt%kset%nkpt))
!  do ik = 1, kkqmt%kset%nkpt
!    write(*,*) "ik=", ik
!    call getevecfv(kkqmtset%kset%vkl(1:3, ik), gkset%vgkl(1:3, 1:gkset%ngkmax, 1, ik),&
!     & evecks(:,:,:,ik))
!  end do
!
!! CONTINUE HERE
!  !write(*,*) "Grids done"
!  !write(*,*) "Writing grids"
!  !call ematgrids_write_grids()
!  !write(*,*) "Writing done"
!
!
!  write(*,*) "Getting plane wave elements (Sepp modified)" 
!  !write(*,*) "Getting eigenvectors"
!  ik = 8
!  iq = 7
!  ikq = qset%ikiq2ikp_nr(ik,iq)
!  write(*,*) "ik=", ik
!  write(*,*) "vkl=", kkqmtset%kset%vkl(1:3, ik)
!  write(*,*) "iq=", iq
!  write(*,*) "vql=", qset%qset%vkl(1:3, iq)
!  write(*,*) "ikq=", ikq
!  write(*,*) "vkql=", kkqmtset%kqmtset%vkl(1:3, ikq)
!  write(*,*) "ngk=", gkset%ngk(1,ik)
!  write(*,*) "ngkq=", gkqmtset%ngk(1,ikq)
!  write(*,*) "ngq=", gqset%ngk(1,iq)
!
!  !write(*,*) "nmatmax,nstfv,nspnfv", nmatmax,nstfv,nspnfv
!  allocate(eveck(nmatmax, nstfv, nspnfv))
!  allocate(eveckq(nmatmax, nstfv, nspnfv))
!  call getevecfv(kkqmtset%kset%vkl(1:3, ik), gkset%vgkl(1:3, 1:gkset%ngkmax, 1, ik),&
!   & eveck)
!  call getevecfv(kkqmtset%kqmtset%vkl(1:3, ikq), gkqmtset%vgkl(1:3, 1:gkqmtset%ngkmax, 1, ikq),&
!   & eveckq)
!  !write(*,*) "Writing eigenvectors"
!  !call writecmplxparts('eveck',revec=dble(eveck),&
!  !  & imvec=aimag(eveck), veclen=size(eveck))
!  !call writecmplxparts('eveckq',revec=dble(eveckq),&
!  !  & imvec=aimag(eveckq), veclen=size(eveckq))
!  myngq = gqset%ngk(1,iq)
!  allocate(emat(nstfv, nstfv, myngq))
!  write(*,*) "Calculate emat"
!  call timesec(t0)
!  !write(*,*) "emat init"
!  call emat_initialize(lmaxapw_ = input%groundstate%lmaxmat) 
!  !write(*,*) "emat init done"
!  !emat = zzero
!  call emat_gen(ik, iq, eveck(:,:,1),eveckq(:,:,1),emat)
! ! write(*,*) "Writing emat to text file"
! ! call emat_write_plain('EMAT', emat)
! ! write(*,*) "Done"
!
!  call emat_finalize()
!  call timesec(t1)
!  write(*,*) "Finalize emat"
!  write(*,*) "s=", t1-t0
!  write(*,*) "done"
!
!!  write(*,*) "Getting plane wave elements (Sepp)"
!!  allocate(ematsepp(nstfv, nstfv, myngq))
!!  ematsepp=zzero
!!  do igq = 1, myngq
!!    write(*,*) "init"
!!    call emat_init(qset%qset%vkl(1:3,iq), gset%ivg(1:3, gqset%igkig(igq, 1, iq)),&
!!      & input%groundstate%lmaxmat, input%xs%lmaxemat)
!!    write(*,*) "gen"
!!    call emat_genemat(ik, 1, nstfv, 1, nstfv, eveck, eveckq, ematsepp(:,:,igq))
!!    write(*,*) "des"
!!    call emat_destroy
!!  end do
!!  write(*,*) "Maximal absolute difference:"
!!  write(*,*) maxval(abs(emat-ematsepp))
!!  write(*,*) "Writing emat to text file"
!!  call emat_write_plain('EMAT_SEPP', ematsepp)
!!  write(*,*) "Done"
!
!  write(*,*) "finalize ematgrids"
!  call ematgrids_finalize()
!  write(*,*) "Deallocating done"
!
!!  write(*,*) "Getting plane wave elements (Sag)" 
!!  call init0
!!  call init1
!!  call init2
!!  call xssave0
!!
!!  allocate(ematsag(nstfv, nstfv, myngq))
!!  call timesec(t0)
!!
!!  call xsgauntgen(max(input%groundstate%lmaxapw, lolmax),&
!!    & input%xs%lmaxemat, max(input%groundstate%lmaxapw, lolmax))
!!  call findgntn0(max(input%xs%lmaxapwwf, lolmax),&
!!    & max(input%xs%lmaxapwwf, lolmax), input%xs%lmaxemat, xsgnt)
!!
!!  call ematrad(iq)
!!  call ematqalloc
!!
!!  ematsag = zzero
!!  bc%n1 = nstfv
!!  bc%n2 = nstfv 
!!  bc%il1 = 1
!!  bc%il2 = 1 
!!  bc%iu1 = nstfv
!!  bc%iu2 = nstfv 
!!  call b_ematqk(iq, ik, ematsag, bc)
!!  call timesec(t1)
!!  write(*,*) "s=", t1-t0
!!  write(*,*) "done"
!!
!!  write(*,*) "Writing sagmeisters emat to text file"
!!  call emat_write_plain('EMAT_SAG', ematsag)
!!  write(*,*) "Done"
!!
!!  write(*,*) "Maximal absolute difference:"
!!  write(*,*) maxval(abs(emat(:,:,:)-ematsag(:,:,:)))
!
!  write(*,*) "Screen bye!"
!
!  call terminate
!  !!+++++++TEST HACK END ++++++++++++!!

  ! Back up number of energy points
  nwdft = nwdf

  ! Change file extension variable in mod_misc to '_SCR.OUT' for the
  ! following calculations
  call genfilname(dotext='_SCR.OUT', setfilext=.true.)

  ! Call dielectric function with only one frequency point
  ! df is a wrapper for dfq(iq) which for 'screen' sets nwdf=1, i.e. static screening
  call df

  ! Alternative for checking only:
  nwdf = nwdft

  if(rank == 0) then
    write(unitout, '(a)') "Info(screen): Screening finished"
  end if
end subroutine screen
!EOC
