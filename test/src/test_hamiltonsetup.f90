/*
 * test_hamiltonsetup.f90
 *
 *  Created on: Aug 11, 2008
 *      Author: chm
 */

subroutine test_hamiltonsetup(apwalm)
 use modfvsystem
 use modmain
 use sclcontroll

implicit none
complex(8),intent(in)::apwalm(:,:,:,:,:)

integer:: ik=1,ispn=1
type (evsystem)::system



call match(ngk(ik,ispn),gkc(1,ik,ispn),tpgkc(1,1,ik,ispn), &
          sfacgk(1,1,ik,ispn),apwalm(1,1,1,1,ispn))
call newsystem(system,packedmatrixstorage,nmat(ik,ispn))
call hamiltonandoverlapsetup(system,ngk(ik,ispn),apwalm,igkig(1,ik,ispn),vgkc(1,1,ik,ispn))
call  deleteystem(system)

end subroutine
