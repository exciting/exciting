/*
 * test_arpacksolver.f90
 *
 *  Created on: Aug 11, 2008
 *      Author: chm
 */

subroutine test_arpacksolver (apwalm)
 use modfvsystem
 use modmain
 use sclcontroll

implicit none
complex(8),intent(in)::apwalm(:,:,:,:,:)
integer:: ik=1,ispn=1
type (evsystem)::system




call newsystem(system,packedmatrixstorage,nmat(ik,ispn))
call hamiltonandoverlapsetup(system,ngk(ik,ispn),apwalm,igkig(1,ik,ispn),vgkc(1,1,ik,ispn))
call deleteystem(system)


end subroutine
