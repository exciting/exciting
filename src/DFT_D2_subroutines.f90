module DFT_D2_subroutines
  implicit none

contains
  
  subroutine loadoldpar(c6ab,r0ab)
    use DFT_D2_parameters, only : au_to_ang, J_to_au, N_A, max_elem
    implicit none
    real(8) :: r0(max_elem),c6(max_elem)
    real(8) :: r0ab(max_elem,max_elem)
    real(8) :: c6ab(max_elem,max_elem)
    
    integer :: i,j
    
    ! the published radii in S.Grimme, J.Comput.Chem. 27, (2006), 1787-1799 (tab 1)
    ! refer to the following values multiplied by 1.1 (rs6 in this code)
    
    r0(1:max_elem) = (/ 0.91d0,0.92d0,&
         0.75d0,1.28d0,1.35d0,1.32d0,1.27d0,1.22d0,1.17d0,1.13d0, &
         1.04d0,1.24d0,1.49d0,1.56d0,1.55d0,1.53d0,1.49d0,1.45d0,&
         1.35d0,1.34d0,&
         1.42d0,1.42d0,1.42d0,1.42d0,1.42d0,&
         1.42d0,1.42d0,1.42d0,1.42d0,1.42d0,&
         1.50d0,1.57d0,1.60d0,1.61d0,1.59d0,1.57d0,&
         1.48d0,1.46d0,&
         1.49d0,1.49d0,1.49d0,1.49d0,1.49d0,&
         1.49d0,1.49d0,1.49d0,1.49d0,1.49d0,&
         1.52d0,1.64d0,1.71d0,1.72d0,1.72d0,1.71d0,&
         1.638d0,1.602d0,1.564d0,1.594d0,1.594d0,1.594d0,1.594d0,&
         1.594d0,1.594d0,1.594d0,1.594d0,1.594d0,1.594d0,1.594d0,&
         1.594d0,1.594d0,1.594d0,&
         1.625d0,1.611d0,1.611d0,1.611d0,1.611d0,1.611d0,1.611d0,&
         1.611d0,&
         1.598d0,1.805d0,1.767d0,1.725d0,1.823d0,1.810d0,1.749d0/)
    
    c6(1:max_elem) = (/0.14d0,0.08d0,&
         1.61d0,1.61d0,3.13d0,1.75d0,1.23d0,0.70d0,0.75d0,0.63d0,&
         5.71d0,5.71d0,10.79d0,9.23d0,7.84d0,5.57d0,5.07d0,4.61d0,&
         10.8d0,10.8d0,10.8d0,10.8d0,10.8d0,&
         10.8d0,10.8d0,10.8d0,10.8d0,10.8d0,10.8d0,10.8d0,16.99d0,&
         17.10d0,16.37d0,12.64d0,12.47d0,12.01d0,24.67d0,24.67d0,&
         24.67d0,24.67d0,24.67d0,24.67d0,24.67d0,24.67d0,24.67d0,&
         24.67d0,24.67d0,24.67d0,37.32d0,38.71d0,38.44d0,31.74d0,&
         31.50d0,29.99d0,315.275d0,226.994d0,176.252d0,&
         140.68d0,140.68d0,140.68d0,140.68d0,140.68d0,140.68d0,140.68d0,&
         140.68d0,140.68d0,140.68d0,140.68d0,140.68d0,140.68d0,140.68d0,&
         105.112d0,&
         81.24d0,81.24d0,81.24d0,81.24d0,81.24d0,81.24d0,81.24d0,&
         57.364d0,57.254d0,63.162d0,63.540d0,55.283d0,57.171d0,56.64d0 /)
    
    c6ab = -1
    do i=1,max_elem
       do j=1,i
          r0ab(i,j)=(r0(i)+r0(j))/au_to_ang
          r0ab(j,i)=(r0(i)+r0(j))/au_to_ang
          c6ab(i,j)=dsqrt(c6(i)*c6(j))
          c6ab(j,i)=dsqrt(c6(i)*c6(j))
       enddo
    enddo
    !convert to au
    c6ab = c6ab * 1d6/J_to_au/(au_to_ang**6)/N_A
    
  end subroutine loadoldpar
  
  subroutine getnumofatoms(numofatoms)
    use modinput
    implicit none
    integer :: icount, numofatoms
    numofatoms = 0
    do icount = 1,size(input%structure%speciesarray)
       numofatoms = numofatoms + size(input%structure%speciesarray(icount)%species%atomarray)
    enddo
  end subroutine getnumofatoms
  
  subroutine getatomdata(numofatoms,xyz,iz)
    use modinput
    implicit none
    integer :: numofatoms
    real(8) :: xyz(3,numofatoms)
    integer :: iz(numofatoms)
    integer :: atomcount
    integer :: icount, jcount, kcount, numberofspecies, numofatomsperspecies, atomicnumber_
    character (len=2) :: elementname
    atomcount = 1
    numberofspecies = size(input%structure%speciesarray)
    do icount = 1,numberofspecies
       numofatomsperspecies = size(input%structure%speciesarray(icount)%species%atomarray)
       elementname = input%structure%speciesarray(icount)%species%speciesfile
       if ( elementname(2:2) .eq. '.')  elementname(2:2) = ' '
       call atomicnumber(elementname, atomicnumber_,'ea')
       iz(atomcount:(atomcount+numofatomsperspecies-1)) = atomicnumber_
       do jcount = 1,numofatomsperspecies
          !        xyz(:,atomcount) = 0
          !        do kcount = 1,3
          !           xyz(:,atomcount) =   xyz(:,atomcount) + &
          !                input%structure%speciesarray(icount)%species%atomarray(jcount)%atom%coord(kcount) * &
          !                input%structure%crystal%basevect(:,kcount)
          !        enddo !kcount
          xyz(:,atomcount) = matmul(input%structure%crystal%basevect, &
               input%structure%speciesarray(icount)%species%atomarray(jcount)%atom%coord(:))
          atomcount = atomcount + 1
       enddo !jcount
    enddo !icount
  end subroutine getatomdata
  
  subroutine atomicnumber(element, atomicnumber_, direction)
    implicit none
    character (len=2) :: element
    integer :: atomicnumber_,i
    character (len=2) :: direction! element name to atomic number: 'ea'
    ! atomic number to element name: 'ae'
    character (len=2), dimension(86) :: all_elements
    all_elements = (/"H ","He",&
         "Li","Be","B ","C ","N ","O ","F ","Ne",&
         "Na","Mg","Al","Si","P ","S ","Cl","Ar",&
         "K ","Ca","Sc","Ti","V ","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",&
         "Rb","Sr","Y ","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I ","Xe",&
         "Cs","Ba",&
         "La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb",&
         "Lu","Hf","Ta","W ","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn"/)
    if (direction == 'ae') then
       if (atomicnumber_ < 1 .or. atomicnumber_ > 86) then
          write(*,*)'SUBROUTINE ATOMICNUMBER: input "atomicnumber_" is not within the allowed range.'
          return
       endif
       element=all_elements(atomicnumber_)
    elseif (direction == 'ea') then
       atomicnumber_ = -1
       do i=1,86
          if (all_elements(i) == element) then
             atomicnumber_=i
             exit
          endif
       enddo
       if (atomicnumber_ == -1) then
          write(*,*)'SUBROUTINE ATOMICNUMBER: input "element" is not included in name list'
          write(*,*)'("element" is a string of length 2, i.e. the string for hydrogen is "H ".'
          write(*,*)'The first letter is always a capital one)'
          return
       endif
    else
       write(*,*)'SUBROUTINE ATOMICNUMBER: input "direction" must either be string "ae" or string "ea"'
       return
    endif
  end subroutine atomicnumber
  
  function cross(a,b)
    real(8) :: cross(3)
    real(8), intent(in) :: a(3), b(3)
    
    cross(1) = a(2)*b(3) - a(3)*b(2)
    cross(2) = a(3)*b(1) - a(1)*b(3)
    cross(3) = a(1)*b(2) - a(2)*b(1)
  end function cross
  
  subroutine getlatticerepetition(latrep)
    use modinput
    use DFT_D2_parameters, only : cutoff
    implicit none
    integer :: latrep(3)
    real(8) :: vec(3,3)
    integer :: icount
      
    !orthogonal system
    vec(:,1) = cross(input%structure%crystal%basevect(:,2),input%structure%crystal%basevect(:,3))
    vec(:,2) = cross(input%structure%crystal%basevect(:,1),input%structure%crystal%basevect(:,3))
    vec(:,3) = cross(input%structure%crystal%basevect(:,1),input%structure%crystal%basevect(:,2))
    do icount = 1,3
       vec(:,icount) = vec(:,icount)/sqrt(dot_product(vec(:,icount),vec(:,icount)))!normalize
       latrep(icount) = int(abs(cutoff/(dot_product(input%structure%crystal%basevect(:,icount),vec(:,icount))))) + 1
    enddo
  end subroutine getlatticerepetition
  
end module DFT_D2_subroutines
