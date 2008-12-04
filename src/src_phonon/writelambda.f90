
subroutine writelambda(wq,gq)
use modmain
implicit none
! arguments
real(8), intent(in) :: wq(3*natmtot,nqpt)
real(8), intent(in) :: gq(3*natmtot,nqpt)
! local variables
integer ik,iq,i
real(8) t1,t2
! get the eigenvalues and occupancies from file
do ik=1,nkpt
  call getevalsv(vkl(:,ik),evalsv(:,ik))
  call getoccsv(vkl(:,ik),occsv(:,ik))
end do
! compute the density of states at the Fermi energy
call occupy
open(50,file='LAMBDAQ.OUT',action='WRITE',form='FORMATTED')
write(50,*)
write(50,'(I4," : total number of atoms")') natmtot
write(50,'(I6," : number of q-points")') nqpt
write(50,*)
do iq=1,nqpt
  write(50,'(I6," : q-point")') iq
  write(50,'(3G18.10," : q-vector (lattice coordinates)")') vql(:,iq)
  write(50,'(3G18.10," : q-vector (Cartesian coordinates)")') vqc(:,iq)
  do i=1,3*natmtot
    t1=pi*fermidos*wq(i,iq)**2
    if (t1.gt.1.d-8) then
      t2=gq(i,iq)/t1
    else
      t2=0.d0
    end if
    write(50,'(I4,G18.10)') i,t2
  end do
  write(50,*)
end do
close(50)
write(*,*)
write(*,'("Info(writelambda):")')
write(*,'(" wrote electron-phonon coupling constants for all q-points to&
 & LAMBDAQ.OUT")')
write(*,*)
return
end subroutine

