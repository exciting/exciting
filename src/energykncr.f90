subroutine energykncr
use modmain
implicit none
integer is,ia,ias,ist,ir
! allocatable local arrays
real(8), allocatable :: rfmt(:,:)
! external functions
real(8) rfmtinp
external rfmtinp
! allocate local arrays
allocate(rfmt(lmmaxvr,nrmtmax))
! calculate the kinetic energy for core states
engykncr=0.d0
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
! sum of core eigenvalues
    do ist=1,spnst(is)
      if (spcore(ist,is)) engykncr=engykncr+spocc(ist,is)*evalcr(ist,ias)
    end do
! core density
    do ir=1,nrmt(is)
      rfmt(1,ir)=rhocr(ir,ias)/y00
    end do
    engykncr=engykncr-rfmtinp(1,0,nrmt(is),spr(1,is),lmmaxvr,rfmt, &
     veffmt(1,1,ias))
  end do
end do
deallocate(rfmt)
return
end subroutine

