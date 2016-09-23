!---------------------------------------------------------------
!BOP
!
!!MODULE: fouri
!
module fouri
!
!!DESCRIPTION:
!
! This module declares the R-vectors and stars for the Fourier transform
!

!!PUBLIC TYPES:
    integer(4) :: nrr                       ! Number of R vectors
    integer(4) :: nst                       ! Number of stars
    integer(4) :: nirk                      ! Number of irreducible k-points
    integer(4) :: irdivk                    ! Division of irred. k-points 
    integer(4), allocatable :: irk(:)       ! Integer index of irred. k-points 
    integer(4), allocatable :: rindex(:,:) 
    integer(4), allocatable :: rst(:,:)
    integer(4), allocatable :: kir(:,:)

    real(8)    :: rmax
    real(8)    :: rbas(3,3)                 ! basis vectors
      
    complex(8), allocatable :: ek(:,:)
    complex(8), allocatable :: er(:,:)
    complex(8), allocatable :: ekpl(:,:)

    logical :: setrindex_done = .false.
      
end module fouri
!EOP

!---------------------------------------------------------------

!BOP
!
!!ROUTINE: setrindex
!
!!INTERFACE:
!
subroutine setrindex
!
!!DESCRIPTION:
!
! Sets the indexes of the real space lattice vectors used for  
! Fourier interpolation ordered by increasing length
!
!!USES:
    use modinput
    use modmain
    use fouri

!!LOCAL VARIABLES:
    implicit none     
    integer(4) :: i         ! (Counter): runs over coordinates
    integer(4) :: ippw      ! (Counter): runs over plane waves
    integer(4) :: ir1       ! (Counter): run over x-coord of G
    integer(4) :: ir2       ! (Counter): run over y-coord of G 
    integer(4) :: ir3       ! (Counter): run over z-coord of G 
    integer(4) :: jppw      ! (Counter): runs over plane waves
    integer(4) :: nr1       ! Max. ir1
    integer(4) :: nr2       ! Max. ir2
    integer(4) :: nr3       ! Max. ir3
    integer(4) :: nr       ! Maximum numb,er of plane waves in intipw
    integer(4) :: nrmin,nrmax,ist,j,isym,lspl
    integer(4), dimension(3) :: irvec ! Integer coordinates of the G-vector
    integer(4), allocatable  :: invrindex(:,:,:)
      
    real(8) :: rr
    real(8), dimension(3) :: rvec     ! Cartesian coordinates of the R-vector
    integer(4), allocatable :: rind(:,:) ! Temporary storage for rindex
    real(8), allocatable :: rlen(:)   ! Temporary storage for all the  qpg's
      
    integer :: ierr
    character(len=10)::sname="setrindex"


! !REVISION HISTORY:
!
! Created May 2004 by RGA
! Last Modified 9. Nov. 2005 by RGA
! Revisited July 2011 by DIN
!
!EOP
!BOC

    ! real-space cutoff parameter (maximal length)
    if (.not.associated(input%gw)) then
      rmax = 120.d0
    else
      rmax = input%gw%rmax
    end if

    do i = 1, 3
      rvec(i) = sqrt(rbas(1,i)**2+rbas(2,i)**2+rbas(3,i)**2)
    enddo
      
    nr1 = 2*nint(rmax/rvec(1))
    nr2 = 2*nint(rmax/rvec(2))
    nr3 = 2*nint(rmax/rvec(3))
    nr = (2*nr1+1)*(2*nr2+1)*(2*nr3+1)

    !write(6,*) '--------- R vectors generation ----------'
    !write(6,*) '  Parameters'
    !write(6,12) rmax,nr1,nr2,nr3,nr

    allocate(rlen(1:nr),rind(3,1:nr))
      
    ippw = 0
    do ir1 = -nr1, nr1
      irvec(1) = ir1
      do ir2 = -nr2, nr2
        irvec(2) = ir2
        do ir3 = -nr3, nr3
          irvec(3) = ir3
          ! Transform to cartesian coordinates
          do i = 1, 3
            rvec(i) = dble(irvec(1))*rbas(i,1)+ &
            &         dble(irvec(2))*rbas(i,2)+ &
            &         dble(irvec(3))*rbas(i,3)
          end do
          rr = sqrt(rvec(1)*rvec(1)+rvec(2)*rvec(2)+rvec(3)*rvec(3))
          if (rr<=rmax) then
            ippw = ippw + 1
            rlen(ippw) = rr
            rind(1:3,ippw) = irvec(1:3)
          endif
        enddo
      enddo
    enddo
    nrr = ippw

    ! sort by increasing length using shell algorithm
    call shelsort(nrr,rind,rlen)

    allocate(rindex(3,1:nrr))
    rindex(:,1:nrr) = rind(:,1:nrr)

    allocate(rst(2,1:nrr))
    rst(:,:) = 0
      
    ! generate the inverse of rindex      
    nrmax = maxval(rindex)+1
    nrmin = minval(rindex)-1
    allocate(invrindex(nrmin:nrmax,nrmin:nrmax,nrmin:nrmax))
    invrindex(nrmin:nrmax,nrmin:nrmax,nrmin:nrmax) = 0
    do ippw = 1, nrr
      ir1 = rindex(1,ippw)
      ir2 = rindex(2,ippw)
      ir3 = rindex(3,ippw)
      invrindex(ir1,ir2,ir3) = ippw
    enddo  

    rst(1,1) = 1
    rst(2,1) = nsymcrys
    ist = 1

    do ippw = 2, nrr
      if (rst(1,ippw)==0) then
        ist = ist+1
        rst(1,ippw) = ist
        do isym = 1, nsymcrys
          lspl = lsplsymc(isym)
          irvec(:) = matmul(symlat(:,:,lspl),rindex(:,ippw))
          jppw = invrindex(irvec(1),irvec(2),irvec(3))
          if (jppw<=0) then
            write(6,11) ippw, rindex(1:3,ippw), isym, irvec(1:3)
          else  
            rst(1,jppw) = ist
            rst(2,jppw) = rst(2,jppw)+1
          endif  
        enddo ! isym
      endif
    enddo ! ippw
    nst = ist

    !write(6,*) " nrr=", nrr
    !write(6,*) " nst=", nst

    deallocate(invrindex,rlen,rind)

    10 format(i6,3i4,4f12.6,2i4)
    11 format('(setrindex) ',2(i6,3i4)) 
    12 format('rmax = ',f15.10,' nr1 =',i4,' nr3 =',i4,' nr3 =',i4,' nr =',i8)
 
end subroutine
!EOC
