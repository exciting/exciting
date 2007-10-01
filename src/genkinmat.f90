subroutine genkinmat(ik,evecfv,evecsv)
use modmain
implicit none
! arguments
integer, intent(in) :: ik
complex(8), intent(in) :: evecfv(nmatmax,nstfv)
complex(8), intent(in) :: evecsv(nstsv,nstsv)
! local variables
integer ist
! allocatable arrays
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: wfmt(:,:,:,:,:)
complex(8), allocatable :: wfir(:,:,:)
complex(8), allocatable :: c(:,:)
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(wfmt(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfir(ngrtot,nspinor,nstsv))
allocate(c(nstsv,nstsv))
! find the matching coefficients
call match(ngk(ik,1),gkc(1,ik,1),tpgkc(1,1,ik,1),sfacgk(1,1,ik,1),apwalm)
! calculate the wavefunctions for all states of the input k-point
call genwfsv(.false.,ngk(ik,1),igkig(1,ik,1),evalsv(1,ik),apwalm,evecfv, &
 evecsv,wfmt,wfir)
! compute effective potential matrix elements
call genvmatk(veffmt,veffir,wfmt,wfir,kinmatc(1,1,ik))
kinmatc(:,:,ik)=-kinmatc(:,:,ik)
! add second-variational states along the diagonal
do ist=1,nstsv
  kinmatc(ist,ist,ik)=kinmatc(ist,ist,ik)+evalsv(ist,ik)
end do
! rotate kinetic matrix elements to Cartesian basis
call zgemm('N','C',nstsv,nstsv,nstsv,zone,kinmatc(1,1,ik),nstsv,evecsv, &
 nstsv,zzero,c,nstsv)
call zgemm('N','N',nstsv,nstsv,nstsv,zone,evecsv,nstsv,c,nstsv,zzero, &
 kinmatc(1,1,ik),nstsv)
deallocate(apwalm,wfmt,wfir,c)
return
end subroutine

