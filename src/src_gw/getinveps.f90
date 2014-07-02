!BOP
!
! !ROUTINE: getinveps
!
! !INTERFACE:
    subroutine getinveps(iq)

! !DESCRIPTION:
!
! This subroutine calculates the inverted dielectric function
!
! !USES:

    use modmain
    use modgw
    use modmpi

! !INPUT PARAMETERS:
      
    implicit none

    integer(4), intent(in) :: iq  ! Index of the q vector for which the
!                                   polarization matrix is calculated
    integer(4) :: iqp
    integer(4) :: isym
    integer(4) :: Recl
    integer(4) :: im
    integer(4) :: iom
    
    real(8)    :: tstart,tend

    complex(8), allocatable :: eps(:,:), tmat(:,:)
    complex(8), allocatable :: rmat(:,:), temp(:,:)
    character(128)::sbuffer
! !EXTERNAL ROUTINES: 

    external zgemm

! !INTRINSIC ROUTINES: 
      
    intrinsic conjg
    intrinsic cpu_time
        
! !REVISION HISTORY:
!
! Created Dec 2011 by DIN
!
!EOP
!BOC

    call cpu_time(tstart)
    
    if(allocated(inveps))deallocate(inveps)
    allocate(inveps(matsizmax,matsizmax,nomeg))

    iqp=indkpq(iq,1)
    inquire(IoLength=Recl) inveps
    write(sbuffer,*)procofindex (iqp, nqpt)
    open(44,file='INVEPS'//trim(adjustl(sbuffer))//'.OUT',action='READ',form='UNFORMATTED', &
   &  access='DIRECT',status='OLD',recl=Recl)

    

    read(44,rec=iqp-firstofset(procofindex (iqp, nqpt),nqpt)+1) inveps

    close(44)
    
    ! if the q-point is non-reduced, the dielectric matrix is calculated by symmetry
    !
    ! Note, iq --isym--> iqp
    !
    isym=iksymq(iq,1)

    if (isym.ne.1) then

      ! Matrix representation of the group symmetry operations in MB
      call genmbrotmat(iq,isym)

!-----------------------------------------------------------------
!     Transform R_{ij} to the eigenvectors of the coulomb matrix
!-----------------------------------------------------------------
      allocate(rmat(mbsiz,mbsiz))      
      allocate(temp(matsiz,mbsiz))
      
      call zgemm('n','n',matsiz,mbsiz,matsiz, &
     &           zone,rotmat,matsiz,vbas,matsiz,zzero,temp,matsiz)

!     here barcvm should be for Rq point
      iqp=idikpq(iqp,1) ! position of iqp in the full q-grid
      call diagsgi(iqp)
      call calcmpwipw(iqp)
      call calcbarcmb(iqp)
      call setbarcev(input%gw%BareCoul%barcevtol)
      
      call zgemm('c','n',mbsiz,mbsiz,matsiz, &
     &           zone,vbas,matsiz,temp,matsiz,zzero,rmat,mbsiz)
      
      deallocate(temp,rotmat)
      deallocate(vbas,barcvm)

!-----------------------------------------------------------------
!     Rotate the dielectric function matrix
!-----------------------------------------------------------------
      allocate(eps(mbsiz,mbsiz))
      allocate(tmat(mbsiz,mbsiz))
     
      do iom = 1, nomeg
     
        eps(1:mbsiz,1:mbsiz)=inveps(1:mbsiz,1:mbsiz,iom)
     
        call zgemm('n','n',mbsiz,mbsiz,mbsiz, &
       &     zone,eps,mbsiz,rmat,mbsiz,zzero,tmat,mbsiz)
      
        call zgemm('c','n',mbsiz,mbsiz,mbsiz, &
       &     zone,rmat,mbsiz,tmat,mbsiz,zzero,eps,mbsiz)
    
        inveps(1:mbsiz,1:mbsiz,iom)=eps(1:mbsiz,1:mbsiz)
     
      end do ! iom
      
      deallocate(tmat)
      deallocate(eps)
      deallocate(rmat)
      
    end if ! isym

    ! data used for treating q->0 singularities
    if (Gamma) then

      if(allocated(head))deallocate(head)
      allocate(head(1:nomeg))

      open(42,file='INVHEAD.OUT',action='READ',form='UNFORMATTED',status='OLD')
      read(42) head
      close(42)

      if(allocated(epsw1))deallocate(epsw1)
      allocate(epsw1(mbsiz,nomeg))
      open(43,file='INVWING1.OUT',action='READ',form='UNFORMATTED',status='OLD')
      read(43) epsw1
      close(43)
      
      if(allocated(epsw2))deallocate(epsw2)
      allocate(epsw2(mbsiz,nomeg))
      open(43,file='INVWING2.OUT',action='READ',form='UNFORMATTED',status='OLD')
      read(43) epsw2
      close(43)

    end if

!--------------------------------------------
!   \epsilon^{-1}_{ij}-\delta_{ij}
!--------------------------------------------

    do iom = 1, nomeg
      if (Gamma) head(iom)=head(iom)-zone
      do im=1,mbsiz
        inveps(im,im,iom)=inveps(im,im,iom)-zone
      enddo  
    enddo ! iom

    call cpu_time(tend)
    if(tend.lt.0.0d0)write(fgw,*)'warning, tend < 0'
    call write_cputime(fgw,tend-tstart,'GETINVEPS')

    end subroutine
!EOC
