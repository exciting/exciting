! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: bandstr
! !INTERFACE:
!
!
subroutine wann_ph (dynmat, ndynmat, nqco, nqcored, qvcco, xqirrdeg, dynmatD, nqfi, qvcfi, nqred)
  ! !USES:
  use mod_wannier
  use m_wsweight
  use mod_constants, only: twopi, zzero, zi

  ! !DESCRIPTION:
  !
  ! !REVISION HISTORY:
  !EOP
  !BOC
  implicit none
  INTEGER, intent( in) :: nqco                                 ! Nr of Q-points on the COarse mesh (irr)
  INTEGER, intent( in) :: nqred (3)                            ! Nr of Q-points on the COarse mesh along x,y,z
  INTEGER, intent( in) :: nqcored                              ! Nr of Q-points on the COarse mesh (red) 
  INTEGER, intent( in) :: nqfi                                 ! Nr of Q-points on the FIne meth (red)
  INTEGER, intent( in) :: ndynmat                              ! 3*N_atoms (size of the dynmat)
  INTEGER, intent( in) :: xqirrdeg (nqco)                      ! degeneracy of the irreducible q-point on the coarse mesh  
  COMPLEX(8), intent( in)  :: dynmat  ( ndynmat,ndynmat,nqco)  ! DYNamical MATrix on the coarse gri
  COMPLEX(8), intent( out) :: dynmatD ( ndynmat,ndynmat,nqfi)  ! DYNamical MATrix on the Dense grid
  REAL(8), intent( in) :: qvcco( 3, nqco)                      ! Q-Vectors (in Cartesian coord) on the COarse mesh (irr)
  REAL(8), intent( in) :: qvcfi( 3, nqfi)                      ! Q-Vectors (in Cartesian coord) on the FIne mesh (red)
  COMPLEX(8), ALLOCATABLE  :: dynmatW (:,:,:)  ! DYNamical MATrix in Wannier representation
  
  INTEGER :: nrpt, ia, ix, iy, iz, iqnr, ik, ir
  REAL(8) :: vc(3), energy, rdotk 
  COMPLEX(8) :: ftweight

  REAL(8),    ALLOCATABLE :: rptc(:,:)
  REAL(8),    ALLOCATABLE :: rptl(:,:)
  COMPLEX(8), ALLOCATABLE :: auxmat(:,:), auxmat2(:,:,:), wanme(:,:,:), evectmp(:,:)

  print *, ' Wannier interpolation of the phonons'
  !--------------------------------
  ! generate set of lattice vectors 
  !--------------------------------
  nrpt = nqcored
  IF (ALLOCATED(rptc)) DEALLOCATE(rptc)
  ALLOCATE(rptc(3,nrpt)) 
  rptc(:,:) = 0.d0 
  IF (ALLOCATED(rptl)) DEALLOCATE(rptl)
  ALLOCATE(rptl(3,nrpt))
  rptl(:,:) = 0.d0 
  !
  ia = 0
  !print *, nqred(3) , nqred(2), nqred(1)
  DO iz = -nqred(3)/2, -nqred(3)/2+nqred(3)-1
    DO iy = -nqred(2)/2, -nqred(2)/2+nqred(2)-1
      DO ix = -nqred(1)/2, -nqred(1)/2+nqred(1)-1
        ia = ia + 1
        rptl( :, ia) = (/ dble( ix), dble( iy), dble( iz)/)
        call r3mv( input%structure%crystal%basevect, rptl( :, ia), rptc( :, ia))
      END DO
    END DO
  END DO
  
  !-------------------------------------
  ! Bring dynmat in Wannier representation 
  !-------------------------------------
  if(allocated(dynmatW)) deallocate(dynmatW)
  allocate(dynmatW(ndynmat,ndynmat,nrpt))
  dynmatW(:,:,:) = zzero
  DO iqnr = 1, nqco
    DO ir = 1, nrpt 
      rdotk = dot_product(rptc(:,ir),qvcco(:,iqnr))
      ftweight = exp(-zi*rdotk)*dble(xqirrdeg(iqnr))/ dble(nqcored)
      dynmatW(:,:,ir) = dynmatW(:,:,ir) + dynmat (:,:,iqnr)*ftweight   ! need the weight of irreducible points
    END DO
  END DO
  !
  !TODO one should check the decay of dynmatW
  !
  !---------------------------------------------
  ! Bring dynmat back to Bloch on the dense grid 
  !---------------------------------------------
  dynmatD (:,:,:) = zzero
  DO ir = 1, nrpt
    DO iqnr = 1, nqfi
      rdotk = twopi * dot_product( qvcfi(:,iqnr), rptc(:,ir))
      ftweight = exp(zi*rdotk) !/dble( ndegen(ir))
      dynmatD (:,:,iqnr) = dynmatD (:,:,iqnr) + ftweight * dynmatW(:,:,ir) 
      !if(iqnr.eq.1 .or. iqnr.eq.nqfi) print *, ' CP4: ir ', ir, 'iqnr',  iqnr  
    ENDDO
  ENDDO

  DEALLOCATE (dynmatW)
  DEALLOCATE (rptl)
  DEALLOCATE (rptc) 

  return
END subroutine wann_ph
