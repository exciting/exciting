
module mod_xsf_format
  use modinput
  use modmain, only : natmtot, nspecies, natoms, atposc, spzn
  implicit none

  real(8), parameter :: bohr2ang = 0.529177d0
  private bohr2ang

contains

!-------------------------------------------------------------------------------
  subroutine write_structure_xsf(fname)
    implicit none
    ! input/output
    character*(*), intent(in) :: fname
    ! local
    integer :: is, ia
    open(80,file=trim(fname),status='Unknown',action='Write')
    write(80,*) 'CRYSTAL'
    write(80,*) 'PRIMVEC'
    write(80,*) input%structure%crystal%basevect(:,1)*bohr2ang
    write(80,*) input%structure%crystal%basevect(:,2)*bohr2ang
    write(80,*) input%structure%crystal%basevect(:,3)*bohr2ang
    write(80,*) 'PRIMCOORD'
    write(80,*) natmtot, ' 1'
    do is = 1, nspecies
    do ia = 1, natoms(is)
      ! write(80,'(A2,3F14.8)') &
      ! &  trim(input%structure%speciesarray(is)%species%chemicalSymbol), &
      ! &  atposc(:,ia,is)*bohr2ang
      write(80,'(I5,3F14.8)') &
      &  abs(int(spzn(is))),  &
      &  atposc(:,ia,is)*bohr2ang
    end do
    end do
    write(80,*)
    close(80)
    return
  end subroutine

!-------------------------------------------------------------------------------
  subroutine write_supercell_xsf(fname,nx,ny,nz)
    implicit none
    ! input/output
    character*(*), intent(in) :: fname
    integer, intent(in) :: nx(2)
    integer, intent(in) :: ny(2)
    integer, intent(in) :: nz(2)
    ! local
    integer :: is, ia, i, j, k
    real(8) :: pos(3)
    open(80,file=trim(fname),status='Unknown',action='Write')
    write(80,*) 'CRYSTAL'
    write(80,*) 'PRIMVEC'
    write(80,*) (nx(2)-nx(1)+1)*input%structure%crystal%basevect(:,1)*bohr2ang
    write(80,*) (ny(2)-ny(1)+1)*input%structure%crystal%basevect(:,2)*bohr2ang
    write(80,*) (nz(2)-nz(1)+1)*input%structure%crystal%basevect(:,3)*bohr2ang
    write(80,*) 'PRIMCOORD'
    write(80,*) natmtot*(nx(2)-nx(1)+1)*(ny(2)-ny(1)+1)*(nz(2)-nz(1)+1), ' 1'
    do is = 1, nspecies
    do ia = 1, natoms(is)
      do k = nz(1), nz(2)
      do j = ny(1), ny(2)
      do i = nx(1), nx(2)
        pos = atposc(:,ia,is)+                               &
        &     dble(i)*input%structure%crystal%basevect(:,1)+ &
        &     dble(j)*input%structure%crystal%basevect(:,2)+ &
        &     dble(k)*input%structure%crystal%basevect(:,3)
        write(80,'(A2,3F14.8)') &
        &  trim(input%structure%speciesarray(is)%species%chemicalSymbol), &
        &  pos*bohr2ang
      end do
      end do
      end do
    end do
    end do
    write(80,*)
    close(80)
    return
  end subroutine
  
  !-------------------------------------------------------------------------------
  subroutine write_2d_xsf(fname,label,boxl,igrid,npt,rdata)
    implicit none
    ! input/output
    character*(*) :: fname
    character*(*) :: label
    real(8), intent(in) :: boxl(3,3)
    integer, intent(in) :: igrid(3)
    integer, intent(in) :: npt
    real(8), intent(in) :: rdata(npt)
    ! local
    integer :: ip, i, is, ia
    real(8) :: boxc(3,3)
    character(30) :: frmt
  
    do i = 1, 3
      call r3mv(input%structure%crystal%basevect, boxl(i,:), boxc(i,:))
    end do

    open(80,file=trim(fname),status='Unknown',action='Write',access='Append')
    write(80,*) 'BEGIN_BLOCK_DATAGRID_2D'
    write(80,*) trim(label)
    write(80,*) 'BEGIN_DATAGRID_2D'
    write(80,*) igrid(1)+1, &
    &           igrid(2)+1
    write(80,*) boxc(1,:)*bohr2ang
    write(80,*) boxc(2,:)*bohr2ang
    write(80,*) boxc(3,:)*bohr2ang
    frmt = "(8e20.10)"
    write(80,trim(frmt)) rdata
    write(80,*) 'END_DATAGRID_2D'
    write(80,*) 'END_BLOCK_DATAGRID_2D'
    close(80)
  
    return
  end subroutine

  !-------------------------------------------------------------------------------
  subroutine write_3d_xsf(fname,label,boxl,igrid,npt,rdata)
    implicit none
    ! input/output
    character*(*) :: fname
    character*(*) :: label
    real(8), intent(in) :: boxl(4,3)
    integer, intent(in) :: igrid(3)
    integer, intent(in) :: npt
    real(8), intent(in) :: rdata(npt)
    ! local
    integer :: ip, i, is, ia
    real(8) :: boxc(4,3)
    character(30) :: frmt

    do i = 1, 4
      call r3mv(input%structure%crystal%basevect, boxl(i,:), boxc(i,:))
    end do

    open(80,file=trim(fname),status='Unknown',action='Write',access='Append')
    write(80,*) 'BEGIN_BLOCK_DATAGRID_3D'
    write(80,*) trim(label)
    write(80,*) 'BEGIN_DATAGRID_3D'
    write(80,*) igrid(1)+1, &
    &           igrid(2)+1, &
    &           igrid(3)+1
    write(80,*) boxc(1,:)*bohr2ang
    write(80,*) boxc(2,:)*bohr2ang
    write(80,*) boxc(3,:)*bohr2ang
    write(80,*) boxc(4,:)*bohr2ang
    frmt = "(8e20.10)"
    write(80,trim(frmt)) rdata
    write(80,*) 'END_DATAGRID_3D'
    write(80,*) 'END_BLOCK_DATAGRID_3D'
    close(80)
  
    return
  end subroutine
  
end module
