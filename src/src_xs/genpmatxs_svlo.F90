!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!>  Calculates the momentum matrix elements in second-variational mode.
!> This routine computes:
!> \f$ p_{ij} = \langle \Psi_{i,\mathbf{k}} \mid -i \nabla \mid \Psi_{j,\mathbf{k}} \rangle \f$
!> for muffin-tin plus interstitial contributions, using local orbitals
!> in a second-variational scheme.

subroutine genpmatxs_svlo(ngp, igpig, vgpc, evecfv, evecsv, pmat)

   use precision,  only: dp
   use modinput
   use modmain
   use modxs,       only: apwcmt, ripaa, ripalo, riploa, riplolo
   use svlo,         only: get_num_of_basis_functions_sv
   use matrix_contraction,  only: contract_A_and_C_with_B_complex_dp 
   use m_ematqk,       only: emat_ccket  
   use mod_variation, only: variation_multiplication
   implicit none

  ! Explanation for each argument
  !> ngp:Number of (G + p) vectors
   Integer, Intent (In) :: ngp
   !> igpig (dimension ngkmax):Index mapping from (G + p)-vectors to G-vectors
   Integer, Intent (In) :: igpig (ngkmax)
   !> vgpc (dimension (3, ngkmax)):(G + p)-vectors in Cartesian coordinates
   real(dp), Intent (In) :: vgpc (3, ngkmax)
   !> evecfv (dimension (nmatmax, nstfv)):First-variational eigenvectors
   complex(dp), Intent (In) :: evecfv (nmatmax, nstfv)
   !> evecsv (dimension (nstsv, nstsv)):Second-variational eigenvectors
   complex(dp), Intent (In) :: evecsv (nstsv, nstsv)
   !> pmat (dimension (3, nstsv, nstsv)):Final momentum matrix elements
   complex(dp), Intent (Out) :: pmat (3, nstsv, nstsv)
  ! local variables
      Integer :: is, ia, ias, ist, jst
      Integer :: ist1, l1, m1, lm1, l3, m3, lm3, io, io1, io2, ilo, &
     & ilo1, ilo2
      Integer :: j, lo_index
      Integer :: igp1, igp2, ig1, ig2, ig, iv1 (3), iv (3)
      Integer :: num_of_basis_functions_sv 
  ! allocatable arrays
      complex(dp), Allocatable :: pm (:, :, :)
      complex(dp), Allocatable :: cfunt (:, :), h (:, :), pmt (:, :)
      complex(dp), Allocatable :: evecfv1 (:, :), evecfv2 (:, :)
      complex(dp), Allocatable :: zv2 (:), zv3(:,:)
      complex(dp), Allocatable :: zv2_lo (:)
      complex(dp), Allocatable :: locmt_zones (:, :, :, :)
  ! external functions
      complex(dp) zfmtinp
      External zfmtinp

      num_of_basis_functions_sv = get_num_of_basis_functions_sv()  

      Allocate (zv2(nstfv))
      Allocate (zv2_lo(nlotot))
      Allocate (cfunt(ngp, ngp))
      Allocate (h(ngp, nstfv))
      Allocate (pmt(nstfv, nstfv))
      Allocate (evecfv1(nstfv, ngp), evecfv2(ngp, nstfv))
      Allocate (pm(num_of_basis_functions_sv, num_of_basis_functions_sv, 3))
      Allocate (locmt_zones(nlotot, nlomax,-lolmax:lolmax, natmtot))
      locmt_zones (:, :, :, :) = zzero


      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            Do ilo1 = 1, nlorb (is)
               l1 = lorbl (ilo1, is)
               Do m1 = - l1, l1
                  lm1 = idxlm (l1, m1)
                  lo_index = idxlo (lm1, ilo1, ias)
                  locmt_zones(lo_index, ilo1, m1, ias) = zone
               End Do
            End Do
         End Do
      End Do
  ! set the momentum matrix elements to zero
      pm (:, :, :) = 0.0_dp
  ! loop over species and atoms
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
        !---------------------------!
        !     APW-APW contribution  !
        !---------------------------!
            Do j = 1, 3
               Do l1 = 0, input%groundstate%lmaxapw
                  Do m1 = - l1, l1
                     lm1 = idxlm (l1, m1)
                     Do io1 = 1, apword (l1, is)
                        zv2 (:) = zzero
                        Do l3 = 0, input%groundstate%lmaxapw
                           Do m3 = - l3, l3
                              lm3 = idxlm (l3, m3)
                              Do io2 = 1, apword (l3, is)
                                 Call zaxpy (nstfv, zone*ripaa(io1, &
                                & lm1, io2, lm3, ias, j), apwcmt(1, &
                                & io2, lm3, ias), 1, zv2, 1)
                              End Do
                           End Do
                        End Do
                        Call zoutpr (nstfv, nstfv, zone, apwcmt(1, io1, &
                       & lm1, ias), zv2, pm(1:nstfv, 1:nstfv, j))
                     End Do
                  End Do
               End Do
            End Do
            If (nlotot .Gt. 0) Then
           !--------------------------------------!
           !     APW-local-orbital contribution   !
           !--------------------------------------!
               Do j = 1, 3
                  Do l3 = 0, input%groundstate%lmaxapw
                     Do m3 = - l3, l3
                        lm3 = idxlm (l3, m3)
                        Do io = 1, apword (l3, is)
                           zv2_lo (:) = zzero
                           Do ilo = 1, nlorb (is)
                              l1 = lorbl (ilo, is)
                              Do m1 = - l1, l1
                                 lm1 = idxlm (l1, m1)
                                 Call zaxpy (nlotot, zone*ripalo(io, &
                                & lm3, ilo, m1, ias, j), locmt_zones(1, ilo, &
                                & m1, ias), 1, zv2_lo, 1)
                              End Do
                           End Do
                           Call zoutpr (nstfv, nlotot, zone, apwcmt(1, &
                          & io, lm3, ias), zv2_lo, pm(1:nstfv, nstfv+1:num_of_basis_functions_sv, j))
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
                        zv2 (:) = zzero
                        Do l3 = 0, input%groundstate%lmaxapw
                           Do m3 = - l3, l3
                              lm3 = idxlm (l3, m3)
                              Do io = 1, apword (l3, is)
                                 Call zaxpy (nstfv, zone*riploa(ilo, &
                                & m1, io, lm3, ias, j), apwcmt(1, io, &
                                & lm3, ias), 1, zv2, 1)
                              End Do
                           End Do
                        End Do
                        Call zoutpr (nlotot, nstfv, zone, locmt_zones(1, ilo, &
                       & m1, ias), zv2, pm(nstfv+1:num_of_basis_functions_sv, 1:nstfv, j))
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
                        zv2_lo (:) = zzero
                        Do ilo2 = 1, nlorb (is)
                           l3 = lorbl (ilo2, is)
                           Do m3 = - l3, l3
                              lm3 = idxlm (l3, m3)
                              Call zaxpy (nlotot, zone*riplolo(ilo1, m1, &
                             & ilo2, m3, ias, j), locmt_zones(1, ilo2, m3, &
                             & ias), 1, zv2_lo, 1)
                           End Do
                        End Do
                        Call zoutpr (nlotot, nlotot, zone, locmt_zones(1, ilo1, &
                       & m1, ias), zv2_lo, pm(nstfv+1:num_of_basis_functions_sv, nstfv+1:num_of_basis_functions_sv, j))
                     End Do
                  End Do
               End Do
           ! end case of local orbitals
            End If
        ! end loop over atoms and species
         End Do
      End Do
  ! multiply y-component with imaginary unit
      pm (:, :, 2) = zi * pm (:, :, 2)
  !  calculate momentum matrix elements in the interstitial region
      Forall (ist1=1:nstfv)
         evecfv1 (ist1, :) = conjg (evecfv(1:ngp, ist1))
      End Forall
      evecfv2 (:, :) = evecfv (1:ngp, :)
      Do j = 1, 3
         Do igp1 = 1, ngp
            ig1 = igpig (igp1)
            iv1 (:) = ivg (:, ig1)
            Do igp2 = 1, ngp
               ig2 = igpig (igp2)
               iv (:) = iv1 (:) - ivg (:, ig2)
               ig = ivgig (iv(1), iv(2), iv(3))
               cfunt (igp1, igp2) = zi * vgpc (j, igp2) * cfunig (ig)
            End Do
         End Do
         Call zgemm ('n', 'n', ngp, nstfv, ngp, zone, cfunt, ngp, &
        & evecfv2, ngp, zzero, h, ngp)
         Call zgemm ('n', 'n', nstfv, nstfv, ngp, zone, evecfv1, nstfv, &
        & h, ngp, zzero, pmt, nstfv)
         pm (1:nstfv, 1:nstfv, j) = pm (1:nstfv, 1:nstfv, j) + pmt (:, :)
      End Do
  ! multiply by -i and set lower triangular part
      Do ist = 1, num_of_basis_functions_sv
         Do jst = ist, num_of_basis_functions_sv
            pm (ist, jst, :) = - zi * pm (ist, jst, :)
            pm (jst, ist, :) = conjg (pm(ist, jst, :))
         End Do
      End Do
  ! compute the second-variational momentum matrix elements
   If (input%groundstate%tevecsv) Then
     allocate(zv3(nstsv,nstsv))
     do j=1,3
        call variation_multiplication( &
             evecsv, pm(:,:,j), evecsv, zv3, &
              dimA=nstsv, dimB=nstsv, startA=1, startB=1 )
      pmat(j,:,:) = zv3(:,:)
      end do
   deallocate(zv3)
   else
    do j=1,3
      pmat(j,:,:) = pm(:,:,j)
     end do
   end if
   deallocate(pm, cfunt, h, pmt, evecfv1, evecfv2, zv2, zv2_lo, locmt_zones)
end Subroutine genpmatxs_svlo

