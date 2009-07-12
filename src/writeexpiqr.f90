
subroutine writeexpiqr
use modmain
implicit none
! local variables
integer nk,ik,jk,i,j
real(8) vecqc(3),a,b
! allocatable arrays
complex(8), allocatable :: emat(:,:)
! initialise universal variables
call init0
call init1
! allocate the matrix elements array for < i,k+G+q | exp(iq.r) | j,k >
allocate(emat(nstsv,nstsv))
! read in the density and potentials from file
call readstate
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
! number of k-points to write out
if (kstlist(1,1).le.0) then
  nk=nkpt
else
  nk=nkstlist
end if
open(50,file='EXPIQR.OUT',action='WRITE',form='FORMATTED')
write(50,*)
write(50,'("q-vector (lattice coordinates) :")')
write(50,'(3G18.10)') vecql
write(50,'("q-vector (Cartesian coordinates) :")')
call r3mv(bvec,vecql,vecqc)
write(50,'(3G18.10)') vecqc
write(50,*)
write(50,'(I8," : number of k-points")') nk
write(50,'(I6," : number of states per k-point")') nstsv
do jk=1,nk
  if (kstlist(1,1).le.0) then
    ik=jk
  else
    ik=kstlist(1,jk)
  end if
  if ((ik.le.0).or.(ik.gt.nkpt)) then
    write(*,*)
    write(*,'("Error(writeexpiqr): k-point out of range : ",I8)') ik
    write(*,*)
    stop
  end if
  write(50,*)
  write(50,'(" k-point (lattice coordinates) :")')
  write(50,'(3G18.10)') vkl(:,ik)
  write(50,*)
  write(50,'(" k-point (Cartesian coordinates) :")')
  write(50,'(3G18.10)') vkc(:,ik)
  call genexpiqr(ik,emat)
  do i=1,nstsv
    write(50,*)
    write(50,'(I6," : state i; state j, <...>, |<...>|^2 below")') i
    do j=1,nstsv
      a=dble(emat(i,j))
      b=aimag(emat(i,j))
      write(50,'(I6,3G18.10)') j,a,b,a**2+b**2
    end do
  end do
! end loop over k-points
end do
close(50)
write(*,*)
write(*,'("Info(writeexpiqr)")')
write(*,'(" < i,k+q | exp(iq.r) | j,k > matrix elements written to&
 & EXPIQR.OUT")')
write(*,'(" for the q-vector in vecql and all k-points in kstlist")')
write(*,*)
deallocate(emat)
return
end subroutine

