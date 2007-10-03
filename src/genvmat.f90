subroutine genvmat(vmt,vir,vmat)
use modmain
implicit none
! arguments
real(8), intent(in) :: vmt(lmmaxvr,nrmtmax,natmtot)
real(8), intent(in) :: vir(ngrtot)
complex(8), intent(out) :: vmat(nstsv,nstsv,nkpt)
! local variables
integer ik
! local arrays
real(8), allocatable :: evalsvl(:)
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: evecfv(:,:)
complex(8), allocatable :: evecsv(:,:)
complex(8), allocatable :: wfmt(:,:,:,:,:)
complex(8), allocatable :: wfir(:,:,:)
! allocate local arrays
allocate(evalsvl(nstsv))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(evecfv(nmatmax,nstfv))
allocate(evecsv(nstsv,nstsv))
allocate(wfmt(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfir(ngrtot,nspinor,nstsv))
do ik=1,nkpt
! get the eigenvectors and values from file
  call getevalsv(vkl(1,ik),evalsvl)
  call getevecfv(vkl(1,ik),vgkl(1,1,ik,1),evecfv)
  call getevecsv(vkl(1,ik),evecsv)
! find the matching coefficients
  call match(ngk(ik,1),gkc(1,ik,1),tpgkc(1,1,ik,1),sfacgk(1,1,ik,1),apwalm)
! calculate the  wavefunctions for all states
  call genwfsv(.false.,ngk(ik,1),igkig(1,ik,1),evalsvl,apwalm,evecfv,evecsv, &
   wfmt,wfir)
  call genvmatk(vmt,vir,wfmt,wfir,vmat(1,1,ik))
end do
deallocate(evalsvl,apwalm,evecfv,evecsv,wfmt,wfir)
return
end subroutine

