subroutine getedata(etype,nparam,ename,pname)
! get eos name and number of parameters
implicit none
! arguments
integer, intent(in) :: etype
integer, intent(out) :: nparam
character(256), intent(out) :: ename(2)
character(256), intent(out) :: pname(*)
select case(etype)
case(1)
  ename(1)="Universal EOS"
  ename(2)="Vinet P et al., J. Phys.: Condens. Matter 1, p1941 (1989)"
  nparam=4
  pname(1)="V0"
  pname(2)="E0"
  pname(3)="B0"
  pname(4)="B0'"
case(2)
  ename(1)="Murnaghan EOS"
  ename(2)="Murnaghan F D, Am. J. Math. 49, p235 (1937)"
  nparam=4
  pname(1)="V0"
  pname(2)="E0"
  pname(3)="B0"
  pname(4)="B0'"
case(3)
  ename(1)="Birch-Murnaghan 3rd-order EOS"
  ename(2)="Birch F, Phys. Rev. 71, p809 (1947)"
  nparam=4
  pname(1)="V0"
  pname(2)="E0"
  pname(3)="B0"
  pname(4)="B0'"
case(4)
  ename(1)="Birch-Murnaghan 4th-order EOS"
  ename(2)="Birch F, Phys. Rev. 71, p809 (1947)"
  nparam=5
  pname(1)="V0"
  pname(2)="E0"
  pname(3)="B0"
  pname(4)="B0'"
  pname(5)="B0''"
case(5)
  ename(1)="Natural strain 3rd-order EOS"
  ename(2)="Poirier J-P and Tarantola A, Phys. Earth Planet Int. 109, p1 (1998)"
  nparam=4
  pname(1)="V0"
  pname(2)="E0"
  pname(3)="B0"
  pname(4)="B0'"
case(6)
  ename(1)="Natural strain 4th-order EOS"
  ename(2)="Poirier J-P and Tarantola A, Phys. Earth Planet Int. 109, p1 (1998)"
  nparam=5
  pname(1)="V0"
  pname(2)="E0"
  pname(3)="B0"
  pname(4)="B0'"
  pname(5)="B0''"
case(7)
  ename(1)="Cubic polynomial in (V-V0)"
  ename(2)=""
  nparam=4
  pname(1)="V0"
  pname(2)="E0"
  pname(3)="B0"
  pname(4)="B0'"
case default
  write(*,*)
  write(*,'("Error(getedata): etype not defined : ",I4)') etype
  write(*,*)
  stop
end select
return
end subroutine
