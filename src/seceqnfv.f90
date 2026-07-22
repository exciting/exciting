!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: seceqnfv
! !INTERFACE:
!
!
Subroutine seceqnfv(ik, ispn, nmatp, ngp, igpig, vgpc, apwalm, sfacgp, tpgpc, cdft_maximum_overlap, evalfv, evecfv)
  
      use modinput,                  only: input
      use precision,                 only: dp, i32
      use cdft,                      only: set_overlap_times_psi_gs
      use constants,                 only: zzero, zone
      use modmpi,                    only: mpiglobal
      use mod_Gkvector,              only: ngkmax
      use mod_APW_LO,                only: apwordmax, apword
      use mod_atoms,                 only: natmtot, natoms, nspecies, idxas
      use mod_muffin_tin,            only: lmmaxapw, idxlm
      use mod_eigensystem,           only: nmatmax, h1on, h1aa, h1loa, h1lolo, mt_hscf, MTRedirect, MaxAPWs
      use mod_eigenvalue_occupancy,  only: nstfv
      use mod_potential_and_density, only: ex_coef
      use modfvsystem,               only: evsystem, newsystem, deletesystem, solvewithlapack
      use mod_hybrids,               only: vnlmat
      use mod_misc,                  only: task
      use mGGA_eigensystem,          only: gen_mGGA_H_and_S, mGGA_H, mGGA_S
      use modmpi,                    only: terminate
      use mod_selfconsistent_gw,     only: is_gw_selfconsistent_flavour, qsgw, gw_first_iteration
      use mod_qsgw,                  only: read_optimized_vxc_to_a_file
      use mod_secular_equation_inversion_symmetry, only: transform_eigenvectors_inversion_symmetry, &
                                                         get_lo_transformation_matrix_inv_sym, &
                                                         solve_secular_equation_inversion_symmetry

  ! !INPUT/OUTPUT PARAMETERS:
  !   nmatp  : order of overlap and Hamiltonian matrices (in,integer(i32))
  !   ngp    : number of G+k-vectors for augmented plane waves (in,integer(i32))
  !   igpig  : index from G+k-vectors to G-vectors (in,integer(i32)(ngkmax))
  !   vgpc   : G+k-vectors in Cartesian coordinates (in,real(3,ngkmax))
  !   apwalm : APW matching coefficients
  !            (in,complex(ngkmax,apwordmax,lmmaxapw,natmtot))
  !   sfacgp : structure factors of G+p-vectors (out,complex(ld,natmtot))
  !   tpgpc  : (theta, phi) coordinates of G+p-vectors (out,real(2,ngkmax))
  !   evalfv : first-variational eigenvalues (out,real(nstfv))
  !   evecfv : first-variational eigenvectors (out,complex(nmatmax,nstfv))
  ! !DESCRIPTION:
  !   Solves the secular equation,
  !   $$ (H-\epsilon O)b=0, $$
  !   for the all the first-variational states of the input $k$-point.
  !
  ! !REVISION HISTORY:
  !   Created March 2004 (JKD)
  !   Revised Oct 2024 (Ronaldo)
  !EOP
  !BOC
      Implicit None
  ! arguments
      integer(i32), intent(in) :: ik
      integer(i32), intent(in) :: ispn
      integer(i32), intent(in) :: nmatp
      integer(i32), intent(in) :: ngp
      integer(i32), intent(in) :: igpig (ngkmax)
      real(dp), intent(in) :: vgpc (3, ngkmax)
      complex(dp), intent(in) :: apwalm (ngkmax, apwordmax, lmmaxapw, natmtot)
      complex(dp), intent(in) :: sfacgp (ngkmax, natmtot)
      real(dp), intent(in) :: tpgpc (2, ngkmax)
      !> If .true., then a constrained DFT calculation with the maximum overlap method
      !> is performed. In this case, the product of the overlap matrix with the GS wavefunctions
      !> must be evaluated and saved
      logical, intent(in) :: cdft_maximum_overlap
      real(dp), intent(out) :: evalfv (nstfv)
      complex(dp), intent(out) :: evecfv (nmatmax, nstfv)
      ! local variables
      type(evsystem) :: system
      logical :: packed
      integer(i32) :: ist
      complex(dp), allocatable :: zm(:,:),zm2(:,:)
      complex(dp), allocatable :: lo_transformation_matrix_inv_sym(:,:)
      real(dp), allocatable :: evec_real(:,:)
      complex(dp), allocatable :: vxcopt(:,:)

      !character( len=64) :: fname
      ! apwi related variables for storing matching coefficients in a more convenient way
      integer(i32) :: is,ia,ias,l,io,m,ifun,lm

  !----------------------------------------!
  !     Hamiltonian and overlap set up     !
  !----------------------------------------!

      packed = input%groundstate%solver%packedmatrixstorage

      if ( associated(input%groundstate%mgga)) then
          Call newsystem (system, packed, nmatp)
          call gen_mGGA_H_and_S( ik, system%hamilton%za, system%overlap%za )

      else if ((input%groundstate%solver%type.ne.'Davidson').or.(input%groundstate%solver%constructHS)) then
        Call newsystem (system, packed, nmatp)
        h1on=(input%groundstate%ValenceRelativity.eq.'iora*')
        call MTRedirect(mt_hscf%main,mt_hscf%spinless)
        Call hamiltonsetup (system, ngp, apwalm, igpig, vgpc)
        Call overlapsetup (system, ngp, apwalm, igpig, vgpc)

        ! If the maximum overlap method is used in a CDFT calculation
        if( cdft_maximum_overlap ) call set_overlap_times_psi_gs( ik, system%overlap%za )

        !------------------------------------------------------------------------!
        !   If Hybrid potential is used apply the non-local exchange potential   !
        !------------------------------------------------------------------------!
        if (task == 7) then
          system%hamilton%za(:,:) = system%hamilton%za(:,:) + &
                                    ex_coef*vnlmat(1:nmatp,1:nmatp,ik)
        end if

        !------------------------------------------------------------------------!
        !   If QSGW here we add the optimized exchange-correlation potential     !
        !------------------------------------------------------------------------!
        if (.not. gw_first_iteration() .and. is_gw_selfconsistent_flavour(qsgw)) then
          call read_optimized_vxc_to_a_file(vxcopt, ik, input%gw%taskGroup%outputFormat, .true.)
          !$omp parallel workshare
          system%hamilton%za(:,:) = system%hamilton%za(:,:) + vxcopt(:,:)
          !$omp end parallel workshare
        end if

      else !expecting matrix-free Davidson here
        nullify(system%hamilton%za)
        nullify(system%overlap%za)
        nullify(system%hamilton%ca)
        nullify(system%overlap%ca)
        nullify(system%hamilton%zap)
        nullify(system%overlap%zap)
        nullify(system%hamilton%cap)
        nullify(system%overlap%cap)
        nullify(system%hamilton%ipiv)
        nullify(system%overlap%ipiv)        
        if (task == 7) then 
          call terminate("seqeqn: cannot run hybrid calculations with iterative"// & 
                         "eigensolver without constructing the Hamiltonian matrix explicitly")
        endif
        
        if (is_gw_selfconsistent_flavour(qsgw)) then
          call terminate("seqeqn: cannot run QSGW calculations with iterative"// & 
                         "eigensolver without constructing the Hamiltonian matrix explicitly")
        end if
        
        call MTRedirect(mt_hscf%main,mt_hscf%spinless)        
      endif

      ! Rearranging matching coefficients should be generalised beyond the Davidson solver and moved to seceqn
      if (input%groundstate%solver%type.eq.'Davidson') then 
        allocate(system%apwi(MaxAPWs(),ngp,natmtot))
        system%apwi=zzero
        Do is = 1, nspecies
          Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            ifun=0
            Do l = 0, input%groundstate%lmaxmat
              Do m = - l, l
                lm = idxlm (l, m)
                Do io = 1, apword (l, is)
                  ifun=ifun+1
                  system%apwi(ifun,1:ngp,ias)=apwalm(1:ngp, io, lm, ias)
                End Do
              End Do
            End Do
          End Do
        End Do

      endif

  !------------------------------------!
  !     solve the secular equation     !
  !------------------------------------!
      if (input%groundstate%solver%type == 'Lapack') then
        Call solvewithlapack(system,nstfv,evecfv,evalfv)
      elseif (input%groundstate%solver%type == 'Davidson') then
        call singularcomponents(mpiglobal, system, ik)
        call davidson(system,nstfv,evecfv,evalfv,ik)
        Call deletesystem (system)
        deallocate(system%apwi)
      elseif (input%groundstate%solver%type == 'inversionsymmetry') then
        lo_transformation_matrix_inv_sym = get_lo_transformation_matrix_inv_sym(ispn,ik)
        Call solve_secular_equation_inversion_symmetry(system, nmatp, ngp, nstfv, nmatmax, &
                lo_transformation_matrix_inv_sym, evalfv, evec_real)
        Call transform_eigenvectors_inversion_symmetry(ngp, lo_transformation_matrix_inv_sym, evec_real, evecfv)
      endif

      if (task == 7) call kinetic_energy(ik,evecfv,apwalm,ngp,vgpc,igpig)

if (input%groundstate%ValenceRelativity == 'iora*') then
! normalise large components
      Call newsystem (system, packed, nmatp)
      h1aa=0d0
      h1loa=0d0
      h1lolo=0d0
      h1on=.false.
!      Call hamiltonandoverlapsetup (system, ngp, apwalm, igpig, vgpc)
      Call overlapsetup (system, ngp, apwalm, igpig, vgpc)
      call olprad
      allocate(zm(nmatp,nstfv))
      allocate(zm2(nstfv,nstfv))


      call zgemm('N', &           ! TRANSA = 'C'  op( A ) = A**H.
                 'N', &           ! TRANSB = 'N'  op( B ) = B.
                  nmatp, &          ! M ... rows of op( A ) = rows of C
                  nstfv, &           ! N ... cols of op( B ) = cols of C
                  nmatp, &          ! K ... cols of op( A ) = rows of op( B )
                  zone, &          ! alpha
                  system%overlap%za, &           ! A
                  nmatp,&           ! LDA ... leading dimension of A
                  evecfv, &           ! B
                  nmatmax, &          ! LDB ... leading dimension of B
                  zzero, &          ! beta
                  zm, &  ! C
                  nmatp &      ! LDC ... leading dimension of C
                  )
      call zgemm('C', &           ! TRANSA = 'C'  op( A ) = A**H.
                 'N', &           ! TRANSB = 'N'  op( B ) = B.
                  nstfv, &          ! M ... rows of op( A ) = rows of C
                  nstfv, &           ! N ... cols of op( B ) = cols of C
                  nmatp, &          ! K ... cols of op( A ) = rows of op( B )
                  zone, &          ! alpha
                  evecfv, &           ! A
                  nmatmax,&           ! LDA ... leading dimension of A
                  zm, &           ! B
                  nmatp, &          ! LDB ... leading dimension of B
                  zzero, &          ! beta
                  zm2, &  ! C
                  nstfv &      ! LDC ... leading dimension of C
                  )

!     write(*,*) zm2(1:2,1:2)
!      do ist=1,nstfv
!        write(*,*) zm2(ist,ist)
!      enddo
!      write(*,*)
      do ist=1,nstfv
        evecfv(:,ist)=evecfv(:,ist)/sqrt(abs(zm2(ist,ist)))
      enddo
      deallocate(zm,zm2)
      Call deletesystem (system)
endif

End Subroutine seceqnfv
!EOC
