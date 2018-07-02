!--------------------------------------
MODULE mod_dynmat 
  !
  IMPLICIT NONE
  !
  INTEGER :: nat 
  INTEGER :: ntyp
  INTEGER :: ndynmat                   ! size of the dynamical matrix 
  INTEGER :: nqredtot                  ! number of reducible points in the BZ (coarse grid)
  INTEGER :: nqirr                     ! number of irreducible points in the BZ  (coarse grid)
  INTEGER :: nqred (3)                 ! nqx, nqy, nqz
  !
  COMPLEX(8),  ALLOCATABLE :: dynmat   (:,:,:) ! the dynamical matrix 
  COMPLEX(8),  ALLOCATABLE :: dynmatD  (:,:,:) ! the dynamical matrix on the dense grid
  REAL   (8),  ALLOCATABLE :: ev       (:,:)   ! the eigenvalues of the dynamical matrix 
  REAL   (8),  ALLOCATABLE :: wphdense (:,:)   ! the eigenvalues of the dynamical matrix on the dense grid 
  REAL   (8),  ALLOCATABLE :: xqcirr   (:,:)   ! reducible q-point mesh in Cartesian coordinates in units of 2\pi/alat 
  REAL   (8),  ALLOCATABLE :: xqcred   (:,:)   ! reducible q-point mesh in Cartesian coordinates in units of 2\pi/alat 
  REAL   (8),  ALLOCATABLE :: xqlred   (:,:)   ! reducible q-point mesh in lattice coordinates in units of 2\pi/alat 
  REAL   (8),  ALLOCATABLE :: amass    (:)
  !
  INTEGER,     ALLOCATABLE :: ityp  (:)
  INTEGER,     ALLOCATABLE :: xqirrdeg(:) ! degeneracy of the irreducible q-point 
  INTEGER,     ALLOCATABLE :: irr2red (:)  ! map index between irreducible and reducible set of q points 
  !            
  REAL(8),     PARAMETER   :: RY_TO_THZ  = 3289.8419608358563
  REAL(8),     PARAMETER   :: RY_TO_CMM1 = 109737.31570111268
  REAL(8),     PARAMETER   :: AMU_RY     = 911.44424213227251
  REAL(8)              :: tpiba 
  !
END MODULE mod_dynmat 
!--------------------------------------

SUBROUTINE read_dyn 

  !----------------------------------------------------------------------------
  !
  USE mod_dynmat, ONLY : ndynmat, nqred, nqredtot, xqcred, xqcirr, nqirr, irr2red, xqirrdeg, dynmat
  ! 
  IMPLICIT NONE
  !
  CHARACTER(LEN=:), ALLOCATABLE :: id
  CHARACTER(LEN=100) :: tmpstr, fname

  INTEGER :: iqirr, iqred, iqred1
 
  tmpstr='diam.dyn'
  

  WRITE(6,*) ' -- Importing the dynamical matrix -- '  

  ! read the total number of q points 
  write(fname,*)'diam.dyn0'
  OPEN(unit=2,file=fname,action='read')
  READ(2,*) nqred
  READ(2,*) nqirr 
  write(6,*)' I found nqirr = ', nqirr
  write(6,*)' I found nqred = ', nqred
  nqredtot = nqred(1)*nqred(2)*nqred(3) 

  ALLOCATE(xqcred(3,nqredtot))
  ALLOCATE(irr2red(nqredtot))
  ALLOCATE(xqcirr(3,nqirr),xqirrdeg(nqirr))

  iqred = 0 
  iqred1 = 0 
  DO iqirr = 1, nqirr  ! loop over all q points 
    ! 
    ! define file names for reading the dynamical matrix 
    if(iqirr > 9 ) write(fname,101)'diam.dyn',iqirr
    if(iqirr < 10) write(fname,102)'diam.dyn',iqirr
    101 format(A8,i2)
    102 format(A8,i1)
    !
    ! import the q-mesh
    CALL load_qmesh (fname,iqirr,iqred) 

    ! read the dynamical matrix from file for each q 
    CALL load_dynmat_red(fname,iqirr,iqred1) 
    !
  ENDDO ! loop over q points found in diam.dyn0 

  !PRINT*, ' I found ', iqred , ' reducible q points '
END SUBROUTINE read_dyn
!----------------------------------------------------------

!----------------------------------------------------------
SUBROUTINE load_qmesh (fname, iqirr, iqred)
  ! 
  USE mod_dynmat, ONLY : xqcred, nqredtot, irr2red, nqirr, xqcirr, xqirrdeg, tpiba
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=100), INTENT(IN)  :: fname
  INTEGER,            INTENT(IN)  :: iqirr 
  INTEGER,            INTENT(OUT) :: iqred 
  !
  ! local variables
  INTEGER :: i , ios, counter
  LOGICAL :: first 
  CHARACTER(LEN=75) :: line
  !

  !IF(iqred==0) PRINT*, 'I am trying to load the q-mesh '
  OPEN(unit=1,file=fname,action='read')
  
  READ(1,*)
  ios=0
  counter = 0
  first = .true.
  DO WHILE (ios==0) 
    !
    READ(1,'(a)',iostat=ios) line
    !
    IF(ios==0 .AND. line(6:14).EQ.'Dynamical') THEN
      !
      counter = counter + 1 
      READ(1,*) 
      READ(1,'(a)',iostat=ios) line
      iqred = iqred + 1 
      READ(line(11:75),*) (xqcred(i,iqred),i=1,3)
      xqcred(:,iqred) = xqcred(:,iqred) * tpiba
      irr2red(iqred) = iqirr
      IF(first) xqcirr(:,iqirr)=xqcred(:,iqred)
      !
      first = .false.
    ENDIF
    !
  ENDDO 
  xqirrdeg(iqirr) = counter 
  !
  IF(iqred .GT. nqredtot) THEN 
    PRINT*, 'There are too many k-points in the dynamical matrix files, check it! Dunno how to deal with this shit. STOP!  '
    STOP 
  ENDIF
  !
  CLOSE(1)
  RETURN
  !
END SUBROUTINE load_qmesh
!----------------------------------------------------------

!----------------------------------------------------------
SUBROUTINE load_dynmat_red (fname, iqirr, iqred)
  ! 
  USE mod_dynmat, ONLY : ndynmat, dynmat, ev, nqirr, nqredtot, tpiba, nat, ityp, amass, ntyp
  USE mod_constants, ONLY : twopi 
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=100), INTENT(IN) :: fname
  INTEGER,            INTENT(IN) :: iqirr 
  INTEGER,            INTENT(OUT) :: iqred

  LOGICAL, SAVE                 :: first =.TRUE.

  REAL (8),         ALLOCATABLE :: tau   (:,:)
  COMPLEX(8),      ALLOCATABLE :: phiq  (:,:,:,:)
  CHARACTER(LEN=3), ALLOCATABLE :: atm(:)

  REAL(8) :: epsil(3,3)
  REAL(8) :: phir(3), phii(3)
  REAL(8) :: xq(3,48), celldm(6), at(3,3) 
  REAL(8) :: tau1(3), amass1, at1(3,3), celldm1(6), q2
  REAL(8), PARAMETER :: eps8=1.D-8
  REAL(8), PARAMETER :: amu_ry=911.4442421322725

  ! I/O variables
  LOGICAL :: lrigid
  INTEGER :: nqs, ibrav, i1, j1
  ! 
  ! local variables
  INTEGER :: ntyp1,nat1,ibrav1,ityp1
  INTEGER :: i, j, na, nb, nt, ios, nta, ntb
  CHARACTER(LEN=3)  :: atm1

  CHARACTER(LEN=75) :: line

  OPEN(unit=1,file=fname,action='read')
  first = .true. 
  
  READ(1,*)
  READ(1,*)

  IF (first) THEN
     !
     ! read cell information from file
     !
     READ(1,*) ntyp,nat,ibrav,(celldm(i),i=1,6)
     !
     tpiba = twopi / celldm(1)
     !
     ALLOCATE (atm(ntyp),tau(3,nat) )
     ALLOCATE ( phiq (3,3,nat,nat) )
     IF(iqirr .EQ. 1)THEN
       ALLOCATE ( ityp(nat),amass(ntyp)) 
     ENDIF

     if (ibrav==0) then
        read (1,'(a)') atm1  ! for compatibility
        read (1,*) ((at(i,j),i=1,3),j=1,3)
     end if
     DO nt = 1,ntyp
        READ(1,*) i,atm(nt),amass(nt)
     END DO
     DO na=1,nat
        READ(1,*) i,ityp(na),(tau(j,na),j=1,3)
     END DO
     !
     !
     IF (iqirr .EQ. 1)THEN
       ndynmat=3*nat 
       IF(ALLOCATED(dynmat))DEALLOCATE(dynmat) 
       ALLOCATE(dynmat(ndynmat,ndynmat,nqredtot)) 
       IF(ALLOCATED(ev))DEALLOCATE(ev)
       ALLOCATE(ev(ndynmat,nqredtot))
       dynmat(:,:,:)=(0.d0,0.d0)
       ev(:,:)=0.d0
     ENDIF
     !
     first=.FALSE.
     !
  ELSE
     !
     ! check cell information with previous one
     !
     READ(1,*) ntyp1,nat1,ibrav1,(celldm1(i),i=1,6)
     if (ibrav==0) then
         read (1,'(a)') atm1 ! for compatibility
         read (1,*) ((at1(i,j),i=1,3),j=1,3)
     end if
     DO nt = 1,ntyp
        READ(1,*) i,atm1,amass1
     END DO
     DO na=1,nat
        READ(1,*) i,ityp1,(tau1(j),j=1,3)
     END DO
  END IF
  amass = amass / AMU_RY 
  !
  ! REPORT READ DATA 
  IF(iqirr.eq.1)THEN
  WRITE(6,*) 'READING THE REDUCIBLE DYNAMICAL MATRIX I FOUND THE FOLLOWING STUFF: ' 
  WRITE(6,*) ' Number of atomic species   : ',  ntyp
  WRITE(6,*) ' Atomic species             : ',  atm
  WRITE(6,*) ' Number of atoms            : ',  nat
  WRITE(6,*) ' Atomic masses found        : ',  amass
  WRITE(6,*) ' Bravais Lattice index (QE) : ',  ibrav 
  WRITE(6,*) ' twopi / alat               : ',  tpiba 
  ENDIF
  !
  nqs = 0
  !
  100 CONTINUE
  READ(1,*,iostat=ios)
  !
  IF(ios==0) READ(1,'(a)',iostat=ios) line
  !
  IF (ios/=0 .or. line(6:14).NE.'Dynamical') GOTO 201 
  ! 
  nqs = nqs + 1
  READ(1,*)
  READ(1,'(a)') line
  READ(line(11:75),*) (xq(i,nqs),i=1,3)
  READ(1,*)
  !
  DO na=1,nat
    DO nb=1,nat
      READ(1,*) i,j
      DO i=1,3
        READ (1,*) (phir(j),phii(j),j=1,3)
        DO j = 1,3
          phiq (i,j,na,nb) = CMPLX(phir(j),phii(j))
        END DO
      END DO
    END DO
  END DO

  iqred = iqred + 1 
  DO na = 1,nat
    DO nb = 1,nat
      DO i= 1,3
        DO j= 1,3
          dynmat((na-1)*3+i,(nb-1)*3+j,iqred) = phiq(i,j,na,nb)
        END DO
      END DO
    END DO
  END DO


  CALL dyndiag2 (nat,ntyp,amass,ityp,dynmat(:,:,iqred),ev(:,iqred))
  !CALL print_ev2 (ev(:,iqred)) 

  ! NOTA BENE: 
  ! If the following GOTO is activated, the subrouting does the diagonalization for all
  ! the reducible points (they are all listed in the dynamical matrix file), if it's commented
  ! only the first point is considered (the eigenvalues are the same for symmetry, but eigenvectors 
  ! should be rotated using the symmetry operations). 
  !
  GOTO 100 
  !
201  CLOSE(1)
  !
END SUBROUTINE load_dynmat_red 

!----------------------------------------------------------

!----------------------------------------------------------
SUBROUTINE print_ev2 (ev)
  !
  USE mod_dynmat, ONLY : ndynmat, RY_TO_THZ, RY_TO_CMM1
  !
  IMPLICIT NONE
  !
  REAL(8),    INTENT(IN) :: ev (ndynmat)
  INTEGER             :: iw 

  PRINT*,' '
  DO iw=1, ndynmat
    WRITE(6,*) iw , SQRT(ABS(ev(iw))) * RY_TO_THZ, ' [THz] ', SQRT(ABS(ev(iw))) * RY_TO_CMM1, ' [cm-1] '
  ENDDO
  
END SUBROUTINE print_ev2

!----------------------------------------------------------
SUBROUTINE print_ev3 (ev,nq,vkl,str)
  !
  USE mod_dynmat, ONLY : ndynmat, RY_TO_THZ, RY_TO_CMM1
  USE m_getunit
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)    :: nq
  REAL(8),    INTENT(IN) :: ev (ndynmat,nq)
  REAL(8),    INTENT(IN) :: vkl(3,nq)
  INTEGER                :: iw , iq, fid
  CHARACTER(1)           :: str


  CALL getunit (fid)
  IF (str .EQ. 'd') OPEN(fid,file='ph-dense.OUT',action='Write',status='Unknown')
  IF (str .EQ. 'c') OPEN(fid,file='ph-coarse.OUT',action='Write',status='Unknown')

  DO iq=1, nq
    WRITE(fid,*)'  ' 
    WRITE(fid,103) iq , nq , vkl(:,iq)  
    DO iw=1, ndynmat 
      WRITE(fid,*) iw , SQRT(ABS(ev(iw,iq))) * RY_TO_THZ, ' [THz] ', SQRT(ABS(ev(iw,iq))) * RY_TO_CMM1, ' [cm-1] ' 
    ENDDO
  ENDDO
  103 format(' iq = ',3x,i4,'/',i4,' | q = ( ',3f8.4,' ) ')

  CLOSE(fid)
  
  RETURN
END SUBROUTINE print_ev3

!-----------------------------------------------------------------------
subroutine dyndiag2 (nat,ntyp,amass,ityp,dyn,ev)
  !-----------------------------------------------------------------------
  !
  !   diagonalise the dynamical matrix
  !   On input:  amass = masses, in amu
  !   On output: ev = energies, z = displacements
  !
  USE mod_dynmat, ONLY: AMU_RY !, dynmat
  
  implicit none
  ! input
  integer nat, ntyp, ityp(nat)
  REAL(8), INTENT (OUT) :: ev (3*nat)   ! the eigenvalues of the dynamical matrix 
  !complex(8) dyn(3,3,nat,nat)
  complex(8) dyn(3*nat,3*nat)
  real(8) amass(ntyp)
  ! output
  !complex(8) z(3*nat,3*nat)
  ! local
  real(8) diff, dif1, difrel
  integer nat3, na, nta, ntb, nb, ipol, jpol, i, j
  complex(8), allocatable :: dyn2(:,:)
  !
  !  fill the two-indices dynamical matrix
  !
  nat3 = 3*nat
  allocate(dyn2 (nat3, nat3)) 
  dyn2 = dyn
  !
  ! impose hermiticity
  diff = 0.d0
  difrel=0.d0
  do i = 1,nat3
     dyn2(i,i) = CMPLX( DBLE(dyn2(i,i)),0.d0)
     do j = 1,i - 1
        dif1 = abs(dyn2(i,j)-CONJG(dyn2(j,i)))
        if ( dif1 > diff .and. &
             max ( abs(dyn2(i,j)), abs(dyn2(j,i))) > 1.0d-6) then
           diff = dif1
           difrel=diff / min ( abs(dyn2(i,j)), abs(dyn2(j,i)))
        end if
        dyn2(i,j) = 0.5d0* (dyn2(i,j)+CONJG(dyn2(j,i)))
        dyn2(j,i) = CONJG(dyn2(i,j))
     end do
  end do
  !if ( diff > 1.d-6 ) write (6,'(5x,"Max |d(i,j)-d*(j,i)| = ",f9.6,/,5x, &
  !     & "Max |d(i,j)-d*(j,i)|/|d(i,j)|: ",f8.4,"%")') diff, difrel*100
  !
  !  divide by the square root of masses
  !
  do na = 1,nat
     nta = ityp(na)
     do nb = 1,nat
        ntb = ityp(nb)
        do ipol = 1,3
           do jpol = 1,3
             dyn2((na-1)*3+ipol, (nb-1)*3+jpol) = &
                  dyn2((na-1)*3+ipol, (nb-1)*3+jpol) / &
                  (amu_ry*sqrt(amass(nta)*amass(ntb)))
          end do
       end do
    end do
 end do
 !
 !
 !print*, 'check 4;'
 !  diagonalisation
 !
 call cdiagh2(nat3,dyn2,nat3,ev)
 !
 deallocate(dyn2)
 !
 !  displacements are eigenvectors divided by sqrt(amass)
 !
! do i = 1,nat3
!    do na = 1,nat
!       nta = ityp(na)
!       do ipol = 1,3
!          z((na-1)*3+ipol,i) = z((na-1)*3+ipol,i)/ sqrt(amu_ry*amass(nta))
!       end do
!    end do
! end do
 !
 return
end subroutine dyndiag2
!
!-----------------------------------------------------------------------
subroutine cdiagh2 (n,h,ldh,e)
  !-----------------------------------------------------------------------
  !
  !   calculates all the eigenvalues and eigenvectors of a complex
  !   hermitean matrix H . On output, the matrix is unchanged
  !
  implicit none
  !
  ! on INPUT
  integer          n,       &! dimension of the matrix to be diagonalized
       &           ldh       ! leading dimension of h, as declared
  ! in the calling pgm unit
  complex(8), intent (in)  :: h(ldh,n)  ! matrix to be diagonalized
  !
  ! on OUTPUT
  real(8)     e(n)      ! eigenvalues
  !complex(8)  v(ldh,n)  ! eigenvectors (column-wise)
  !
  ! LOCAL variables (LAPACK version)
  !
  integer          lwork,   &! aux. var.
       &           ILAENV,  &! function which gives block size
       &           nb,      &! block size
       &           info      ! flag saying if the exec. of libr. routines was ok
  !
  real(8), allocatable::   rwork(:)
  complex(8), allocatable:: work(:)
  !
  !     check for the block size
  !
  nb = ILAENV( 1, 'ZHETRD', 'U', n, -1, -1, -1 )
  if (nb.lt.1) nb=max(1,n)
  if (nb.eq.1.or.nb.ge.n) then
     lwork=2*n-1
  else
     lwork = (nb+1)*n
  endif
  !print*, 'check 5.0;', nb
  !
  !if ( lwork .LT. 10 ) lwork = 10
  if(allocated(work))deallocate(work)
  allocate(work (lwork))
  if(allocated(rwork))deallocate(rwork)
  allocate(rwork (3*n-2))

  !print*, 'check 5;'
  call ZHEEV('N','U',n,h,ldh,e,work,lwork,rwork,info)
  deallocate(rwork)
  deallocate(work)
  !
  return
end subroutine cdiagh2
