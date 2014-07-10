!BOP
! !ROUTINE: task_dos

! !INTERFACE:
subroutine task_dos

!!DESCRIPTION:
!
! Calculate G0W0 DOS data.
! 
!!USES:
    use modinput
    use modmain
    use modgw
    implicit none

!!LOCAL VARIABLES:
    integer(4) :: ik, ist, nstsv0
    integer(4) :: iw, n, nsk(3)
    integer(4) :: nsmdos, nwdos, ngrdos
    real(8) :: winddos(2)
    real(8) :: dw, egap
    real(8), allocatable :: e(:,:)
    real(8), allocatable :: f(:,:)
    real(8), allocatable :: w(:)
    real(8), allocatable :: g(:)

!!REVISION HISTORY:
! Created Mar 2014 by DIN

!EOP
!BOC
    !-----------------
    ! Initialization
    !-----------------
    call init0
    call init1
    
    ! read KS energies 
    do ik = 1, nkpt
      call getevalsv(vkl(:,ik),evalsv(:,ik))
    end do
   
    ! read QP energies from file and perform Fourier interpolation (if required)
    call getevalqp(nkpt,vkl,evalsv)
    
    ! generate the k-grid for the tetrahedron integration library
    !input%groundstate%stypenumber = -1
    !call init1
    ! calculate Fermi energy
    !n = min(int(chgval/2.d0)+10,nstsv)
    !call fermi(nkpt,n,evalsv(ibgw:n,:),ntet,tnodes,wtet,tvol, &
    !&          chgval,.false.,efermi,egap)
    call occupy

    ! GW number of states
    nstsv = nbgw-ibgw+1
    allocate(e(nstsv,nkpt))
    ! shift Efermi to zero
    e(1:nstsv,:) = evalsv(ibgw:nbgw,:)-efermi
      
    ! dos element of properties
    if (.not.associated(input%properties%dos)) &
    &  input%properties%dos => getstructdos(emptynode)
    nsmdos = input%properties%dos%nsmdos
    nwdos = input%properties%dos%nwdos
    ngrdos = input%properties%dos%ngrdos
    winddos(:) = input%properties%dos%winddos(:)    
    
    ! generate energy grid
    allocate(w(nwdos))
    dw = (winddos(2)-winddos(1))/dble(nwdos)
    do iw = 1, nwdos
      w(iw) = dw*dble(iw-1)+winddos(1)
    end do
      
    ! number of subdivisions used for interpolation
    nsk(:) = max(ngrdos/input%groundstate%ngridk(:),1)
    
    ! diagonal of spin density matrix for weight
    allocate(f(nstsv,nkpt))
    f(:,:) = 1.d0
    
    ! BZ integration
    allocate(g(nwdos)) ! DOS
    call brzint(nsmdos, input%groundstate%ngridk, nsk, ikmap, &
    &  nwdos, winddos, nstsv, nstsv, e, f, g)

    ! output file
    open(50,file='DOS-QP.OUT',action='WRITE',form='FORMATTED')
    do iw = 1, nwdos
      write(50,'(2G18.10)') w(iw), g(iw)
    end do
    close(50)  
    
    deallocate(e,w,f,g)
 
    return
end subroutine
!!EOC
