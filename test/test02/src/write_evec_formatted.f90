subroutine write_evec_formatted()
use modmain
implicit none
  complex(8), allocatable :: evecfv(:,:,:)
integer ik
allocate(evecfv(nmatmax,nstfv,nspnfv))
 open(50,file="EIGVECTORFV.OUT")
do ik=1,nkpt
call getevecfv(vkl(:,ik),vgkl(:,:,ik,1),evecfv)
write(50,*)evecfv
end do
close(50)
deallocate (evecfv)
end subroutine
