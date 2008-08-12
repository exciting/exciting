subroutine iterativearpacksecequn(ik,ispn,apwalm,vgpc,evalfv,evecfv)

  !USES:
  use modmain
  use modmpi
  use sclcontroll
  use modfvsystem
  ! !INPUT/OUTPUT PARAMETERS:
  !   ik     : k-point number (in,integer)
  !   ispn   : first-variational spin index (in,integer)
  !   apwalm : APW matching coefficients
  !            (in,complex(ngkmax,apwordmax,lmmaxapw,natmtot))
  !   vgpc   : G+k-vectors in Cartesian coordinates
  !   evalfv : first-variational eigenvalues (out,real(nstfv))
  !   evecfv : first-variational eigenvectors (out,complex(nmatmax,nstfv))
  ! !DESCRIPTION:
  ! This routine will perform several ARPACK iterations

  !BOC
  implicit none
#ifdef DEBUG
  !include declarations for timing output of ARPACK
#include "./debugf90.h"
#endif
  ! arguments
  integer,intent(in)       :: ik
  integer,intent(in)       :: ispn
  real(8),intent(in)       :: vgpc(3,ngkmax)
  complex(8),intent(in)    :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
  real(8),intent(inout)    :: evalfv(nstfv,nspnfv)
  complex(8),intent(inout) :: evecfv(nmatmax,nstfv,nspnfv)

  ! local variables
type (evsystem)::system
  complex(8)::sigma
  integer ::n
  real:: cpu0,cpu1,cpu2
  Complex(8)::                 zero, one
  parameter         (zero = (0.0D+0, 0.0D+0) ,one = (1.0D+0, 0.0D+0) )
  !IO vars
  integer::koffset,recl
  character(256):: outfilenamestring,filetag
  external outfilenamestring




 if(lowesteval.eq.-1.d0) then
  call minenergy(sigma)
 else
 sigma=dcmplx(lowesteval,0)
 endif


  !##################
  !setup hamiltonian#
  !##################

  call newsystem(system,packedmatrixstorage,n)
  call hamiltonandoverlapsetup(system,ngk(ik,ispn),apwalm,igkig(1,ik,ispn),vgpc)
  call cpu_time(cpu0)

  call arpacksolve(system,sigma,evecfv(:,:,ispn), evalfv(:,ispn))
  call deleteystem(system)

  return
end subroutine iterativearpacksecequn
!EOC
