!BOP
!
!!ROUTINE: plotevecprod
!
!!INTERFACE:
    subroutine plotevecprod(ik,jk,ib1,ib2,atom1,atom2)

!!DESCRIPTION:
!
!This subroutine calculates the real space representation of the product
!of two eigenvectors in the line joining the two given atoms for ploting. 

!!USES:
    use modmain
    use modgw
    use mod_rpath

!!INPUT PARAMETERS:
    implicit none
    integer(4), intent(in) :: ik
    integer(4), intent(in) :: jk
    integer(4), intent(in) :: ib1 ! THe band index
    integer(4), intent(in) :: ib2 ! THe band index
    integer(4), intent(in) :: atom1 ! the atom used as origin
    integer(4), intent(in) :: atom2 ! the atom used as final position
      
!!LOCAL VARIABLES:
    integer(4) :: ia, ia1, ia2
    integer(4) :: ir, jr, kr
    integer(4) :: is, is1, is2
    integer(4) :: i, j
    complex(8), allocatable :: evecprod1(:)
    complex(8), allocatable :: evecprod2(:)
    character(len=64) :: filename
    
    call boxmsg(6,'-','PLOTEVECPROD')
      
    write(*,*) 'Parameters:'
    write(*,*) '1 k-point number (iik): ', ik
    write(*,*) '2 k-point number (jjk): ', jk
    write(*,*) '1 band index (ib1): ', ib1
    write(*,*) '2 band index (ib2): ', ib2
    write(*,*) 'atom 1 (at1): ', atom1
    write(*,*) 'atom 2 (at2): ', atom2
    write(*,*)

    if ((atom1<1).or.(atom1>natmtot)) stop 'atom1 is wrong'
    if ((atom2<1).or.(atom2>natmtot)) stop 'atom2 is wrong'

    ! Set the name of the output file
5   format('evecs-',i4,'-',i4,'-',i4,'-',i4,'-',i4,'-',i4,'.out')      
    write(filename,5) ik, jk, ib1, ib2, atom1, atom2
    call str_strip(filename)
    open(unit=71,file=filename,status='unknown')

6   format('eprod-',i4,'-',i4,'-',i4,'-',i4,'-',i4,'-',i4,'.out')      
    write(filename,6) ik, jk, ib1, ib2, atom1, atom2
    call str_strip(filename)
    open(unit=72,file=filename,status='unknown')
    
    allocate(evecprod1(rpath%nptot))
    allocate(evecprod2(rpath%nptot))
    call calcevecplot(ik,ib1,evecprod1)
    call calcevecplot(jk,ib2,evecprod2)
      
    do i = 1, rpath%nptot
      write(71,'(5g18.10)') rpath%r(i,1), real(evecprod1(i)), aimag(evecprod1(i)), &
      &                     real(evecprod2(i)), aimag(evecprod2(i))
      evecprod1(i) = evecprod1(i)*conjg(evecprod2(i))
      write(72,'(4g18.10)') rpath%r(i,1), real(evecprod1(i)), aimag(evecprod1(i)), &
      &                     abs(evecprod1(i))   
    end do
    
    deallocate(evecprod1)
    deallocate(evecprod2)
    
    close(71)
    close(72)  
      
    return
end subroutine  
!EOC
