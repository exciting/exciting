! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
! Copyright (C) Exciting Code, SOL group. 2020

! HISTORY:
!
! Created April 2019 (Ronaldo Rodrigues Pela)
! Improved documentation: July 2021 (Ronaldo)
! Reference: https://doi.org/10.1088/2516-1075/ac0c26

!> This module deals with the momentum matrix considering as basis (L)APW+lo
module rttddft_pmat
  use modmpi
  use modxs, only: ripaa, ripalo, riploa, riplolo
  use mod_gkvector, only: ngk, ngkmax, gkc, vgkc, igkig
  use mod_kpoint, only: nkpt
  use mod_eigensystem, only: nmat, nmatmax, idxlo
  use mod_atoms, only: nspecies, natoms, idxas, natmtot
  use mod_APW_LO, only: apword, apwordmax, nlorb, lorbl, nlotot, nlomax, lolmax
  use mod_muffin_tin, only: idxlm, lmmaxapw
  use modinput, only: input
  use mod_gvector, only: ivg, ivgig, cfunig
  use rttddft_GlobalVariables, only: pmat, apwalm, pmatmt
  use constants, only: zzero, zone, zi
  use precision, only: dp

  implicit none
  private 

  public :: Obtain_Pmat_LAPWLOBasis
contains

  !> Here, we calculate the momentum matrix elements considering as basis 
  !> (L)APW+lo. We copied most of the code from `/src/src_xs/genpmatxs.F90`, 
  !> but there the basis are the KS-wavefunctions
  subroutine Obtain_Pmat_LAPWLOBasis( make_hermitian, evaluate_pmat_mt )
    !> If .True., for each `ik`, force the x, y, and z components of `pmat` to be hermitian
    logical,intent(in)        :: make_hermitian
    logical,intent(in)        :: evaluate_pmat_mt

    integer                   :: ik, first_kpt, last_kpt

    pmat(:,:,:,:) = zzero
    if ( evaluate_pmat_mt ) pmatmt(:,:,:,:,:) = zzero


    if(allocated(ripaa)) deallocate(ripaa)
    allocate(ripaa(apwordmax, lmmaxapw, apwordmax, lmmaxapw,natmtot, 3))
    if(nlotot .gt. 0) then
      if(allocated(ripalo)) deallocate(ripalo)
      allocate(ripalo(apwordmax, lmmaxapw, nlomax,-lolmax:lolmax, natmtot, 3))
      if(allocated(riploa)) deallocate(riploa)
      allocate(riploa(nlomax,-lolmax:lolmax, apwordmax, lmmaxapw, natmtot, 3))
      if(allocated(riplolo)) deallocate(riplolo)
      allocate(riplolo(nlomax,-lolmax:lolmax, nlomax,-lolmax:lolmax, natmtot, 3))
    end if

    ! Calculate gradient of radial functions times spherical harmonics
    call pmatrad

    call distribute_loop(mpi_env_k, nkpt, first_kpt, last_kpt)
    
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(NONE), PRIVATE(ik) SHARED(make_hermitian,pmat,pmatmt) &
!$OMP& SHARED(rank,first_kpt,last_kpt,evaluate_pmat_mt,apwalm,ripaa,ripalo,riploa,riplolo)
!$OMP DO
#endif
    do ik = first_kpt, last_kpt
      if( evaluate_pmat_mt ) then
        call generate_pmat_ik( ik, apwalm(:,:,:,:,ik), make_hermitian, pmat(:,:,:,ik), &
          pmatmt(:,:,:,:,ik) )
      else
        call generate_pmat_ik( ik, apwalm(:,:,:,:,ik), make_hermitian, pmat(:,:,:,ik) )
      end if
    end do
#ifdef USEOMP
!$OMP END DO NOWAIT
!$OMP END PARALLEL
#endif

    deallocate(ripaa)
    if(nlotot .gt. 0) then
      deallocate(ripalo)
      deallocate(riploa)
      deallocate(riplolo)
    end if

  end subroutine Obtain_Pmat_LAPWLOBasis

  !> Obtain the momentum matrix elements for a given a `k-point`
  subroutine generate_pmat_ik(ik, apwalmk, make_hermitian, pmatk, pmat_mt)
    !> Index of the `k-point` considered
    integer, intent(in)       :: ik
    !> apwalmk:  The matching coefficients for the (L)APW's (coefficient that
    !>           smoothly connect LAPW's and their MT counterparts)
    !>           Dimensions: ngkmax,apwordmax,lmmaxapw,natmtot
    complex(dp), intent(in)   :: apwalmk(:, :, :, :)
    !> If .True., force the x, y, and z components of `pmatk` to be hermitian
    logical,intent(in)        :: make_hermitian
    !> Momentum matrix elements for this given `k-point`. 
    !> Dimensions: `nmatmax`, `nmatmax`, `3`.
    complex(dp), intent(out)  :: pmatk(:, :, :)
    !> MT part of the momentum matrix
    complex(dp), intent(out), optional   :: pmat_mt(nmatmax,nmatmax,3,natmtot)


    integer                 :: is, ia, ias, io, io1, io2, ilo, ilo1, ilo2
    integer                 :: ig1, igp1, ig2, igp2, ig, iv1 (3), iv (3)
    integer                 :: j, l1, m1, lm1, l2, m2, lm2
    integer                 :: nmatp, ng
    logical                 :: evaluate_pmat_mt
    complex(dp), Allocatable :: zv (:), aux(:,:)

    nmatp = nmat(1,ik)
    ng = ngk(1,ik)

    allocate (zv(ng))
    allocate (aux(ng,ng))

    evaluate_pmat_mt = present( pmat_mt )

    ! loop over species and atoms
    Do is = 1, nspecies
      Do ia = 1, natoms (is)
        ias = idxas (ia, is)
        !---------------------------!
        !     APW-APW contribution  !
        !---------------------------!
        Do j = 1, 3 ! (loop over the x,y,z components)
          aux(:,:) = zzero
          Do l1 = 0, input%groundstate%lmaxapw
            Do m1 = - l1, l1
              lm1 = idxlm (l1, m1)
              Do io1 = 1, apword (l1, is)
                zv (:) = zzero
                Do l2 = 0, input%groundstate%lmaxapw
                  Do m2 = - l2, l2
                    lm2 = idxlm (l2, m2)
                    Do io2 = 1, apword (l2, is)
                      zv(1:ng) = zv(1:ng) + ripaa(io1,lm1,io2,lm2,ias,j)*apwalmk(1:ng,io2,lm2,ias)
                    End Do
                  End Do
                End Do
                Call zoutpr(ng,ng,zone,apwalmk(1:ng,io1,lm1,ias),zv(1:ng),aux(1:ng,1:ng))
              End Do
            End Do
          End Do
          pmatk(1:ng,1:ng,j) = pmatk(1:ng,1:ng,j) + aux(1:ng,1:ng)
          if ( evaluate_pmat_mt ) pmat_mt(1:ng,1:ng,j,ias) = aux(1:ng,1:ng)
        End Do !Do j = 1, 3
        if (nlotot .gt. 0) then
        !--------------------------------------!
        !     APW-local-orbital contribution   !
        !--------------------------------------!
          Do j = 1, 3
            Do l1 = 0,input%groundstate%lmaxapw ! (loop over lapw)
              Do m1 = - l1, l1
                lm1 = idxlm(l1, m1)
                Do io = 1, apword(l1, is)
                  Do ilo = 1, nlorb(is) ! (loop over local orbitals)
                    l2 = lorbl(ilo, is)
                    Do m2 = - l2, l2
                      lm2 = idxlm(l2, m2)
                      zv(1:ng) = ripalo(io,lm1,ilo,m2,ias,j)*conjg(apwalmk(1:ng,io,lm1,ias))
                      pmatk(1:ng,ng+idxlo(lm2,ilo,ias),j) = pmatk(1:ng,ng+idxlo(lm2,ilo,ias),j) + &
                        & zv(1:ng)
                      if ( evaluate_pmat_mt ) pmat_mt(1:ng,ng+idxlo(lm2,ilo,ias),j,ias) = &
                          & pmat_mt(1:ng,ng+idxlo(lm2,ilo,ias),j,ias) + zv(1:ng)
                    End Do
                  End Do
                End Do
              End Do
            End Do
          End Do
          !--------------------------------------!
          !     local-orbital-APW contribution   !
          !--------------------------------------!
          Do j = 1, 3
            Do ilo = 1, nlorb (is)
              l1 = lorbl (ilo, is)
              Do m1 = - l1, l1
                lm1 = idxlm (l1, m1)
                Do l2 = 0, input%groundstate%lmaxapw
                  Do m2 = - l2, l2
                    lm2 = idxlm(l2, m2)
                    Do io = 1, apword(l2, is)
                      zv(1:ng) = riploa(ilo,m1,io,lm2,ias,j)*apwalmk(1:ng,io,lm2,ias)
                      pmatk(ng+idxlo(lm1,ilo,ias),1:ng,j) = &
                        & pmatk(ng+idxlo(lm1,ilo,ias),1:ng,j) + zv(1:ng)
                      if ( evaluate_pmat_mt ) pmat_mt(ng+idxlo(lm1,ilo,ias),1:ng,j,ias) = &
                          & pmat_mt(ng+idxlo(lm1,ilo,ias),1:ng,j,ias) + zv(1:ng)
                    End Do
                  End Do
                End Do
              End Do
            End Do
          End Do
          !------------------------------------------------!
          !     local-orbital-local-orbital contribution   !
          !------------------------------------------------!
          Do j = 1, 3
            Do ilo1 = 1, nlorb (is)
              l1 = lorbl (ilo1, is)
              Do m1 = - l1, l1
                  lm1 = idxlm (l1, m1)
                  Do ilo2 = 1, nlorb (is)
                    l2 = lorbl (ilo2, is)
                    Do m2 = - l2, l2
                      lm2 = idxlm (l2, m2)
                      pmatk(ng+idxlo(lm1,ilo1,ias),ng+idxlo(lm2,ilo2,ias),j) = &
                        & riplolo(ilo1,m1,ilo2,m2,ias,j)
                      if ( evaluate_pmat_mt ) pmat_mt(ng+idxlo(lm1,ilo1,ias),ng+idxlo(lm2,ilo2,ias),j,ias) = &
                          & riplolo(ilo1,m1,ilo2,m2,ias,j)
                    End Do
                  End Do
              End Do
            End Do
          End Do
        ! end case of local orbitals
        End If
      ! end loop over atoms and species
      End Do
    End Do

    pmatk(:,:,1) = -zi*pmatk(:,:,1)
    pmatk(:,:,3) = -zi*pmatk(:,:,3)
    if ( evaluate_pmat_mt ) then
      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas(ia,is)
          pmat_mt(:,:,1,ias) = -zi*pmat_mt(:,:,1,ias)
          pmat_mt(:,:,3,ias) = -zi*pmat_mt(:,:,3,ias)
        end do
      end do
    end if

    !  calculate momentum matrix elements in the interstitial region
    do j = 1, 3
      do igp1 = 1, ng
        ig1 = igkig(igp1,1,ik)
        iv1 (:) = ivg (:, ig1)
        do igp2 = 1, ng
          ig2 = igkig(igp2,1,ik)
          iv (:) = iv1 (:) - ivg (:, ig2)
          ig = ivgig (iv(1), iv(2), iv(3))
          pmatk(igp1,igp2,j) = pmatk(igp1,igp2,j) + vgkc(j,igp2,1,ik)*cfunig(ig)
        end do
      end do
    end do

    ! This forces the matrix to be hermitian
    if( make_hermitian ) then
      do j = 1, 3 
        do ig1 = 1, nmatp
          do ig2 = ig1+1, nmatp
            pmatk(ig2,ig1,j) = conjg(pmatk(ig1,ig2,j))
          end do
        end do
      end do
    end if

  end subroutine

end module rttddft_pmat
