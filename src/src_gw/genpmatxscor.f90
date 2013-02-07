!
!
!
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: genpmatxs
! !INTERFACE:
!
!
Subroutine genpmatxscor(ik,ngp,apwalm,evecfv,evecsv,pmatc)
! !USES:

      Use modinput
      Use modmain
      Use modxs, Only: apwcmt
      Use modgw

! !INPUT/OUTPUT PARAMETERS:
!   ngp    : number of G+p-vectors (in,integer)
!   igpig  : index from G+p-vectors to G-vectors (in,integer(ngkmax))
!   vgpc   : G+p-vectors in Cartesian coordinates (in,real(3,ngkmax))
!   evecfv : first-variational eigenvector (in,complex(nmatmax,nstfv))
!   evecsv : second-variational eigenvectors (in,complex(nstsv,nstsv))
!   pmat   : momentum matrix elements (out,complex(3,nstsv,nstsv))
! !DESCRIPTION:
!   Calculates the momentum matrix elements
!   $$ p_{ij}=\langle\Psi_{i,{\bf k}}|-i\nabla|\Psi_{j,{\bf k}}\rangle. $$
!   The gradient is applied explicitly only to the radial functions and
!   corresponding spherical harmonics for the muffin-tin part. In the
!   interstitial region the gradient is evaluated analytically.
!   Parts taken from the routine {\tt genpmat}.
!
! !REVISION HISTORY:
!
!   Created July 2011 based on the routine src_xs/genpmatxs.f90 (by S.Sagmeister)
!
!EOP
!BOC
      Implicit None
!     arguments
      integer, intent(in)     :: ik
      integer, intent(in)     :: ngp
      complex(8), intent(in)  :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
      complex(8), intent(in)  :: evecfv(nmatmax,nstfv)
      complex(8), intent(in)  :: evecsv(nstsv,nstsv)
      complex(8), intent(out) :: pmatc(3,ncg,nstsv)
!     local variables
      Integer :: is, ia, ias, ist, j
      Integer :: l1, m1, lm1, l3, m3, lm3, io, icor, icg
      
      real(8)    :: arg
      complex(8) :: phs
      
      logical :: lhermit=.true.

!     allocatable arrays
      Complex(8), Allocatable :: pmcv(:,:,:)
      Complex(8), Allocatable :: pmvc(:,:,:)
      allocate(pmcv(nstsv,nclm,3))
      allocate(pmvc(nstsv,nclm,3))

!     loop over species and atoms
      icg = 0
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas(ia, is)
           !--------------------------------------!
           !     APW-core-orbital contribution   !
           !--------------------------------------!
            pmvc(:,:,:) = zzero
            Do j = 1, 3
               Do l3 = 0, input%groundstate%lmaxapw
                  Do m3 = -l3, l3
                     lm3 = idxlm(l3, m3)
                     Do io = 1, apword(l3,is)
                        Do icor = 1, ncore(is)
                           l1 = spl(icor,is)
                           Do m1 = -l1, l1
                              lm1 = idxlm(l1,m1)
                              
                              do ist = 1, nstfv
                                 pmvc(ist,lm1,j) = pmvc(ist,lm1,j) + &
                                &  ripacor(io,lm3,icor,lm1,ias,j)  * &
                                &  conjg(apwcmt(ist,io,lm3,ias))
                              end do
                              
                           End Do
                        End Do
                     End Do
                  End Do
               End Do
            End Do
           !--------------------------------------!
           !     core-orbital-APW contribution   !
           !--------------------------------------!
            pmcv(:,:,:) = zzero
            Do j = 1, 3
               Do icor = 1, ncore(is)
                  l1 = spl(icor,is)
                  Do m1 = -l1, l1
                     lm1 = idxlm(l1, m1)
                     Do l3 = 0, input%groundstate%lmaxapw
                        Do m3 = - l3, l3
                           lm3 = idxlm (l3, m3)
                           Do io = 1, apword (l3, is)
                              
                              do ist = 1, nstfv
                                 pmcv(ist,lm1,j) = pmcv(ist,lm1,j)       + &
                                &  zone*ripcora(icor,lm1,io,lm3,ias,j) * &
                                &  apwcmt(ist,io,lm3,ias)
                              end do

                           End Do
                        End Do
                     End Do
                  End Do
               End Do
            End Do

           ! multiply y-component with imaginary unit
           pmvc(:,:,2) = zi*pmvc(:,:,2)
           pmcv(:,:,2) = zi*pmcv(:,:,2)

           ! multiply by -i
           pmvc(:,:,:) = -zi*pmvc(:,:,:)
           pmcv(:,:,:) = -zi*pmcv(:,:,:)
           
           ! add the core wavefunction phases
           arg = atposl(1,ia,is)*vkl(1,ik)+ &
                 atposl(2,ia,is)*vkl(2,ik)+ &
                 atposl(3,ia,is)*vkl(3,ik)
           
           phs = cmplx(cos(2.0d0*pi*arg),sin(2.0d0*pi*arg),8)
           pmvc(:,:,:) = pmvc(:,:,:)*phs
           
           phs = cmplx(cos(2.0d0*pi*arg),-sin(2.0d0*pi*arg),8)
           pmcv(:,:,:) = pmcv(:,:,:)*phs
           
           ! the output array
           do icor = 1, ncore(is)
             icg = icg+1 ! counter over total amount of core states

             if(lhermit)then
               pmatc(1,icg,:) = 0.5d0*(pmcv(:,icor,1)+conjg(pmvc(:,icor,1)))
               pmatc(2,icg,:) = 0.5d0*(pmcv(:,icor,2)+conjg(pmvc(:,icor,2)))
               pmatc(3,icg,:) = 0.5d0*(pmcv(:,icor,3)+conjg(pmvc(:,icor,3)))
             else           
               pmatc(1,icg,:) = pmcv(:,icor,1)
               pmatc(2,icg,:) = pmcv(:,icor,2)
               pmatc(3,icg,:) = pmcv(:,icor,3)
             end if
           
           end do ! icor

        ! end loop over atoms and species
         End Do
      End Do

      deallocate(pmcv,pmvc)
      return
End Subroutine genpmatxscor
!EOC
