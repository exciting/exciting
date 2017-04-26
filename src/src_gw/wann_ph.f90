!
subroutine wann_ph (nqcored,qvlcored,nqfi,qvlfi,nqred,wphdense)
  use mod_wannier
  use m_wsweight
  use mod_constants, only: twopi, zzero, zi
  use mod_dynmat,    only: AMU_RY, amass, nat, ityp, ndynmat, dynmat, dynmatD, ntyp

  implicit none
  REAL(8), intent( out) :: wphdense(ndynmat,nqfi)
  INTEGER, intent( in) :: nqred (3)                            ! Nr of Q-points on the COarse mesh along x,y,z
  INTEGER, intent( in) :: nqcored                              ! Nr of Q-points on the COarse mesh (red) 
  INTEGER, intent( in) :: nqfi                                 ! Nr of Q-points on the FIne meth (red)
  COMPLEX(8), ALLOCATABLE  :: dynmatW (:,:,:)  ! DYNamical MATrix in Wannier representation
  REAL(8), intent( in) :: qvlcored ( 3, nqcored)               ! Q-Vectors (in Lattice coord) on the COarse mesh (red)
  REAL(8), intent( in) :: qvlfi    ( 3, nqfi)                  ! Q-Vectors (in Lattice coord) on the FIne mesh (red)
  
  INTEGER :: nrpt, ia, ix, iy, iz, iqnr, ik, ir, fid
  REAL(8) :: vc(3), energy, rdotk 
  COMPLEX(8) :: ftweight

  INTEGER,    ALLOCATABLE :: rptl(:,:)
  COMPLEX(8), ALLOCATABLE :: auxmat(:,:), auxmat2(:,:,:), wanme(:,:,:), evectmp(:,:)
  integer nat3, na, nta, ntb, nb, i, j 

  print *, ' Wannier interpolation of the phonons'

  !--------------------------------
  ! generate set of lattice vectors 
  !--------------------------------
  nrpt = nqcored

  IF (ALLOCATED(rptl)) DEALLOCATE(rptl)
  ALLOCATE(rptl(3,nrpt))
  rptl(:,:) = 0.d0 
  !
  ia = 0
  DO iz = -nqred(3)/2, -nqred(3)/2+nqred(3)-1
    DO iy = -nqred(2)/2, -nqred(2)/2+nqred(2)-1
      DO ix = -nqred(1)/2, -nqred(1)/2+nqred(1)-1
        ia = ia + 1
        rptl( :, ia) = (/ ix, iy, iz/)
      END DO
    END DO
  END DO
   
  !-------------------------------------
  ! Bring dynmat in Wannier representation 
  !-------------------------------------
  if(allocated(dynmatW)) deallocate(dynmatW)
  allocate(dynmatW(ndynmat,ndynmat,nrpt))
  dynmatW(:,:,:) = zzero
  DO iqnr = 1, nqcored
    DO ir = 1, nrpt
      rdotk = twopi * dot_product(qvlcored(:,iqnr),dble(rptl(:,ir)))
      ftweight = exp(-zi*rdotk) !/ dble(nqcored)
      dynmatW(:,:,ir) = dynmatW(:,:,ir) + dynmat (:,:,iqnr)*ftweight   ! need the weight of irreducible points
    END DO
  END DO
  !TODO one should check the decay of dynmatW
  !
  !---------------------------------------------
  ! Bring dynmat back to Bloch on the dense grid
  !---------------------------------------------
  dynmatD (:,:,:) = zzero
  DO ir = 1, nrpt
    DO iqnr = 1, nqfi
      rdotk = twopi * dot_product(  qvlfi(:,iqnr), dble(rptl(:,ir)))
      ftweight = exp(zi*rdotk) / nrpt 
      dynmatD (:,:,iqnr) = dynmatD (:,:,iqnr) + ftweight * dynmatW(:,:,ir)
    ENDDO
  ENDDO

  wphdense(:,:) = 0.d0 
  DO iqnr = 1, nqfi
    CALL dyndiag2 (nat,ntyp,amass,ityp,dynmatD(:,:,iqnr),wphdense(:,iqnr))
  END DO

!  CALL getunit(fid)
!  OPEN (fid,file='rdotk2.OUT',action='Write',status='Unknown')
!  DO ir = 1, nrpt
!    rdotk = dot_product(rptl(:,ir),qvlfi(:,2))
!    write(fid,103) ir, rptl(:,ir), qvlfi(:,2),  rdotk
!  END DO
!  CLOSE(fid)

  IF(ALLOCATED(dynmatW))DEALLOCATE (dynmatW)
  IF(ALLOCATED(rptl))DEALLOCATE (rptl)
 
  RETURN
END subroutine wann_ph
