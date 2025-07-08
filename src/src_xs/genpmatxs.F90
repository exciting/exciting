!
!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!> Calculates the momentum matrix elements
!> \[
!>   p_{ij}=\langle\Psi_{i,{\bf k}}|-i\nabla|\Psi_{j,{\bf k}}\rangle.
!> \]
!> The gradient is applied explicitly only to the radial functions and
!> corresponding spherical harmonics for the muffin-tin part. In the
!> interstitial region the gradient is evaluated analytically.
!> Parts taken from the routine `[[genpmat]]`.
subroutine genpmatxs (ngp, igpig, vgpc, evecfv, evecsv, pmat)
      use precision, only: dp
      use modinput
      use modmain
      use modxs, Only: apwcmt, locmt, ripaa, ripalo, riploa, riplolo
      use svlo, only: get_num_of_basis_functions_sv   
      use xlapack, only: matrix_multiply
      use matrix_contraction, only: contract_A_and_C_with_B_complex_dp
      use m_ematqk, only: emat_ccket
      use mod_variation, only: variation_multiplication
! !REVISION HISTORY:
!   Created April 2008 (Sagmeister)
  implicit none
  ! Explanation for each argument
  !> ngp:Number of (G + p) vectors
  integer,   intent(in)           :: ngp
   !> igpig (dimension ngkmax):Index from (G + p) vectors to G vectors
  integer,   intent(in)           :: igpig(ngkmax)
  !> vgpc (dimension (3, ngkmax)):(G + p) vectors in Cartesian coordinates
  real(dp),  intent(in)           :: vgpc(3, ngkmax)    
  !> evecfv (dimension (nmatmax, nstfv)):First-variational eigenvectors
  complex(dp), intent(in)         :: evecfv(nmatmax, nstfv)
  !> evecsv (dimension (nstsv, nstsv)):Second-variational eigenvectors
  complex(dp), intent(in)         :: evecsv(nstsv, nstsv)
  !> pmat (dimension (3, nstsv, nstsv)):Momentum matrix elements
  complex(dp), intent(out)        :: pmat(3, nstsv, nstsv)

  ! local variables
  integer :: ispn, is, ia, ias, ist, jst
  integer :: ist1, l1, m1, lm1, l3, m3, lm3, io, io1, io2, ilo, &
             ilo1, ilo2
  integer :: i, j, k, l
  integer :: igp1, igp2, ig1, ig2, ig, iv1 (3), iv (3)
  complex(dp) :: zt1, zv (3)
  integer :: num_of_basis_functions_sv  
  ! allocatable arrays
  complex(dp), allocatable :: wfmt (:, :, :)
  complex(dp), allocatable :: gwfmt (:, :, :, :)
  complex(dp), allocatable :: pm (:, :, :)
  complex(dp), allocatable :: cfunt (:, :), h (:, :), pmt (:, :)
  complex(dp), allocatable :: evecfv1 (:, :), evecfv2 (:, :)
  complex(dp), allocatable :: zv2 (:), zv3(:,:)

  ! external functions
  complex(dp) zfmtinp
  External zfmtinp

  num_of_basis_functions_sv = get_num_of_basis_functions_sv()

  allocate (zv2(nstfv))
  allocate (wfmt(lmmaxapw, nrcmtmax, nstfv))
  allocate (gwfmt(lmmaxapw, nrcmtmax, 3, nstfv))
  allocate (cfunt(ngp, ngp))
  allocate (h(ngp, nstfv))
  allocate (pmt(nstfv, nstfv))
  allocate (evecfv1(nstfv, ngp), evecfv2(ngp, nstfv))
  allocate (pm(nstfv, nstfv, 3))
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
                       & lm1, ias), zv2, pm(1, 1, j))
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
                           zv2 (:) = zzero
                           Do ilo = 1, nlorb (is)
                              l1 = lorbl (ilo, is)
                              Do m1 = - l1, l1
                                 lm1 = idxlm (l1, m1)
                                 Call zaxpy (nstfv, zone*ripalo(io, &
                                & lm3, ilo, m1, ias, j), locmt(1, ilo, &
                                & m1, ias), 1, zv2, 1)
                              End Do
                           End Do
                           Call zoutpr (nstfv, nstfv, zone, apwcmt(1, &
                          & io, lm3, ias), zv2, pm(1, 1, j))
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
                        Call zoutpr (nstfv, nstfv, zone, locmt(1, ilo, &
                       & m1, ias), zv2, pm(1, 1, j))
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
                        zv2 (:) = zzero
                        Do ilo2 = 1, nlorb (is)
                           l3 = lorbl (ilo2, is)
                           Do m3 = - l3, l3
                              lm3 = idxlm (l3, m3)
                              Call zaxpy (nstfv, zone*riplolo(ilo1, m1, &
                             & ilo2, m3, ias, j), locmt(1, ilo2, m3, &
                             & ias), 1, zv2, 1)
                           End Do
                        End Do
                        Call zoutpr (nstfv, nstfv, zone, locmt(1, ilo1, &
                       & m1, ias), zv2, pm(1, 1, j))
                     End Do
                  End Do
               End Do
           End If
        End Do
     End Do
  ! multiply y-component with imaginary unit
      pm (:, :, 2) = zi * pm (:, :, 2)
  ! calculate momentum matrix elements in the interstitial region
  evecfv1 = conjg(transpose(evecfv(1:ngp, :)))
  evecfv2 (:, :) = evecfv (1:ngp, :)
  Do j = 1, 3
    cfunt(:,:) = 0.0_dp   
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
     h = zzero
     call matrix_multiply(cfunt, evecfv2, h)
     pmt = zzero
     call matrix_multiply(evecfv1, h, pmt)
     pm(:,:,j) = pm(:,:,j) + pmt
  End Do
  
  ! multiply by -i and set lower triangular part
      Do ist = 1, nstfv
         Do jst = ist, nstfv
            pm (ist, jst, :) = - zi * pm (ist, jst, :)
            pm (jst, ist, :) = conjg (pm(ist, jst, :))
         End Do
      End Do
  ! compute the second-variational momentum matrix elements
      If (input%groundstate%tevecsv) Then
        allocate(zv3(nstsv,nstsv))
        do j=1,3
          call variation_multiplication( &
          evecsv, pm(:,:,j), evecsv, zv3,         &
          dimA=nstsv, dimB=nstsv, startA=1, startB=1 )
          pmat(j,:,:) = zv3(:,:)
         end do
   deallocate(zv3)
      else
        do j=1,3
           pmat(j,:,:) = pm(:,:,j)
        end do
      End if
   deallocate(wfmt, gwfmt, pm, cfunt, h, pmt, evecfv1, evecfv2)
end subroutine genpmatxs


