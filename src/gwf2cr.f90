
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gwf2cr(gwf2mt)
use modmain
implicit none
! arguments
real(8), intent(inout) :: gwf2mt(lmmaxvr,nrmtmax,natmtot)
! local variables
integer is,ia,ias,ist,i
integer l,m,lm,ir,itp
real(8) t1
! automatic arrays
complex(8) zftp(lmmaxvr)
! allocatable arrays
complex(8), allocatable :: zfmt(:,:)
complex(8), allocatable :: gzfmt(:,:,:)
allocate(zfmt(lmmaxvr,nrmtmax))
allocate(gzfmt(lmmaxvr,nrmtmax,3))
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do ist=1,spnst(is)
      if (spcore(ist,is).and.(spk(ist,is).eq.spl(ist,is)+1)) then
        l=spl(ist,is)
        do m=-l,l
          lm=idxlm(l,m)
          zfmt(:,1:nrmt(is))=0.d0
          do ir=1,nrmt(is)
            zfmt(lm,ir)=rwfcr(ir,1,ist,ias)/spr(ir,is)
          end do
          call gradzfmt(lmaxvr,nrmt(is),spr(:,is),lmmaxvr,nrmtmax,zfmt,gzfmt)
          do i=1,3
            do ir=1,nrmt(is)
              call zgemv('N',lmmaxvr,lmmaxvr,zone,zbshtvr,lmmaxvr, &
               gzfmt(:,ir,i),1,zzero,zftp,1)
              do itp=1,lmmaxvr
                t1=dble(zftp(itp))**2+aimag(zftp(itp))**2
! factor of 2 from spin
                gwf2mt(itp,ir,ias)=gwf2mt(itp,ir,ias)+2.d0*t1
              end do
            end do
          end do
        end do
      end if
    end do
! end loops over atoms and species
  end do
end do
deallocate(zfmt,gzfmt)
return
end subroutine

