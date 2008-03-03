
! Copyright (C) 2007 F. Bultmark, F. Cricchio, L. Nordstrom and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine symdmatlu(dmat)
use modmain
implicit none
! arguments
complex(8), intent(inout) :: dmat(lmmaxlu,lmmaxlu,nspinor,nspinor,natmtot)
! local variables
integer isym,lspl,lspn
integer ispn,jspn
integer is,ia,ja,ias,jas
integer lm1,lm2
real(8) det,t1
! automatic arrays
logical done(natmmax)
! allocatable arrays
complex(8), allocatable :: zflm(:,:)
complex(8), allocatable :: ulm(:,:,:)
complex(8), allocatable :: su2(:,:,:)
complex(8), allocatable :: dm1(:,:,:,:,:)
complex(8), allocatable :: dm2(:,:)
complex(8), allocatable :: dm3(:,:,:,:)
complex(8), allocatable :: dm4(:,:)
complex(8), allocatable :: dm5(:,:)
if (ldapu.eq.0) return
! allocate local arrays
allocate(zflm(lmmaxlu,lmmaxlu))
allocate(ulm(lmmaxlu,lmmaxlu,nsymlat))
allocate(su2(nspinor,nspinor,nsymlat))
allocate(dm1(lmmaxlu,lmmaxlu,nspinor,nspinor,natmmax))
allocate(dm2(lmmaxlu,lmmaxlu))
allocate(dm3(lmmaxlu,lmmaxlu,nspinor,nspinor))
allocate(dm4(nspinor,nspinor),dm5(nspinor,nspinor))
! setup a complex unit matrix for (l,m) components
zflm(:,:)=0.d0
do lm1=1,lmmaxlu
  zflm(lm1,lm1)=1.d0
end do
do isym=1,nsymlat
! construct (l,m) rotation matrix for each lattice symmetry
  call rotzflm(symlatc(1,1,isym),lmaxlu,lmmaxlu,lmmaxlu,zflm,ulm(1,1,isym))
! construct SU(2) matrix for proper rotation of spinor components
! (note that rotsu2 uses only the proper part of the rotation matrix)
  if (spinpol) then
    call rotsu2(symlatc(1,1,isym),det,su2(1,1,isym))
  end if
end do
t1=1.d0/dble(nsymcrys)
do is=1,nspecies
! make copy of the density matrices
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    dm1(:,:,:,:,ia)=dmat(:,:,:,:,ias)
  end do
  done(:)=.false.
  do ia=1,natoms(is)
    if (.not.done(ia)) then
      ias=idxas(ia,is)
      dmat(:,:,:,:,ias)=0.d0
      do isym=1,nsymcrys
        lspl=lsplsymc(isym)
        lspn=lspnsymc(isym)
! equivalent atom index (symmetry rotates atom ja into atom ia)
        ja=ieqatom(ia,is,isym)
        jas=idxas(ja,is)
! apply (l,m) symmetry matrix as U*D*conjg(U')
        do ispn=1,nspinor
          do jspn=1,nspinor
            call zgemm('N','N',lmmaxlu,lmmaxlu,lmmaxlu,zone,ulm(1,1,lspl), &
             lmmaxlu,dm1(1,1,ispn,jspn,ja),lmmaxlu,zzero,dm2,lmmaxlu)
            call zgemm('N','C',lmmaxlu,lmmaxlu,lmmaxlu,zone,dm2,lmmaxlu, &
             ulm(1,1,lspl),lmmaxlu,zzero,dm3(1,1,ispn,jspn),lmmaxlu)
          end do
        end do
! apply SU(2) symmetry matrix as U*D*conjg(U') and add
        if (spinpol) then
          do lm1=1,lmmaxlu
            do lm2=1,lmmaxlu
              dm4(:,:)=dm3(lm1,lm2,:,:)
              call z2mm(su2(1,1,lspn),dm4,dm5)
              call z2mmct(dm5,su2(1,1,lspn),dm4)
              dmat(lm1,lm2,:,:,ias)=dmat(lm1,lm2,:,:,ias)+dm4(:,:)
            end do
          end do
        else
          dmat(:,:,1,1,ias)=dmat(:,:,1,1,ias)+dm3(:,:,1,1)
        end if
! end loop over crystal symmetries
      end do
! normalise
      dmat(:,:,:,:,ias)=t1*dmat(:,:,:,:,ias)
      done(ia)=.true.
! rotate into equivalent atoms
      do isym=1,nsymcrys
        ja=ieqatom(ia,is,isym)
        if (.not.done(ja)) then
          jas=idxas(ja,is)
          lspl=lsplsymc(isym)
          lspn=lspnsymc(isym)
! apply (l,m) symmetry matrix as conjg(U')*D*U (rotates atom ia into atom ja)
          do ispn=1,nspinor
            do jspn=1,nspinor
              call zgemm('C','N',lmmaxlu,lmmaxlu,lmmaxlu,zone,ulm(1,1,lspl), &
               lmmaxlu,dmat(1,1,ispn,jspn,ias),lmmaxlu,zzero,dm2,lmmaxlu)
              call zgemm('N','N',lmmaxlu,lmmaxlu,lmmaxlu,zone,dm2,lmmaxlu, &
               ulm(1,1,lspl),lmmaxlu,zzero,dm3(1,1,ispn,jspn),lmmaxlu)
            end do
          end do
! apply SU(2) symmetry matrix as conjg(U')*D*U
          if (spinpol) then
            do lm1=1,lmmaxlu
              do lm2=1,lmmaxlu
                dm4(:,:)=dm3(lm1,lm2,:,:)
                call z2mctm(su2(1,1,lspn),dm4,dm5)
                call z2mm(dm5,su2(1,1,lspn),dm4)
                dmat(lm1,lm2,:,:,jas)=dm4(:,:)
              end do
            end do
          else
            dmat(:,:,1,1,jas)=dm3(:,:,1,1)
          end if
          done(ja)=.true.
        end if
      end do
    end if
  end do
end do
deallocate(zflm,ulm,su2)
deallocate(dm1,dm2,dm3,dm4,dm5)
return
end subroutine

