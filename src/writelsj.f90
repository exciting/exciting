
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writelsj
use modmain
implicit none
! local variables
integer ispn,is,ia,nr,ir
integer ik,ist,kst,i,j,n
real(8) xl(3),xs(3),xj(3)
complex(8) zt1
! allocatable arrays
logical, allocatable :: done(:)
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: evecfv(:,:)
complex(8), allocatable :: evecsv(:,:)
complex(8), allocatable :: wfmt1(:,:,:)
complex(8), allocatable :: wfmt2(:,:,:)
complex(8), allocatable :: wfmt3(:,:,:,:)
complex(8), allocatable :: wfmt4(:,:,:,:)
complex(8), allocatable :: zlflm(:,:)
! external functions
complex(8) zfmtinp
external zfmtinp
! initialise universal variables
call init0
call init1
allocate(done(nstfv))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(evecfv(nmatmax,nstfv))
allocate(evecsv(nstsv,nstsv))
allocate(wfmt1(lmmaxapw,nrmtmax,nstfv))
allocate(wfmt2(lmmaxapw,nrmtmax,2))
allocate(wfmt3(lmmaxapw,nrmtmax,2,3))
allocate(wfmt4(lmmaxapw,nrmtmax,2,3))
allocate(zlflm(lmmaxapw,3))
! read density and potentials from file
call readstate
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
open(50,file='LSJ.OUT',action='WRITE',form='FORMATTED')
write(50,*)
write(50,'("(expectation values computed only over the muffin-tin)")')
! begin loop over k-point and state list
do kst=1,nkstlist
  ik=kstlist(1,kst)
  j=kstlist(2,kst)
  if ((ik.le.0).or.(ik.gt.nkpt)) then
    write(*,*)
    write(*,'("Error(writelsj): k-point out of range : ",I8)') ik
    write(*,*)
    stop
  end if
  if ((j.le.0).or.(j.gt.nstsv)) then
    write(*,*)
    write(*,'("Error(writelsj): state out of range : ",I8)') j
    write(*,*)
    stop
  end if
! get the eigenvectors from file
  call getevecfv(vkl(1,ik),vgkl(1,1,ik,1),evecfv)
  call getevecsv(vkl(1,ik),evecsv)
! find the matching coefficients
  call match(ngk(ik,1),gkc(1,ik,1),tpgkc(1,1,ik,1),sfacgk(1,1,ik,1),apwalm)
  write(50,*)
  write(50,*)
  write(50,'("k-point : ",I6,3G18.10)') ik,vkl(:,ik)
  write(50,'("State : ",I6)') j
  write(50,*)
  do is=1,nspecies
    nr=nrmt(is)
    n=lmmaxapw*nr
    do ia=1,natoms(is)
      write(50,'("Species : ",I4," (",A,"), atom : ",I4)') is, &
       trim(spsymb(is)),ia
      done(:)=.false.
! generate spinor wavefunction from second-variational eigenvectors
      wfmt2(:,:,:)=0.d0
      i=0
      do ispn=1,nspinor
        do ist=1,nstfv
          i=i+1
          zt1=evecsv(i,j)
          if (abs(dble(zt1))+abs(aimag(zt1)).gt.epsocc) then
            if (.not.done(ist)) then
              call wavefmt(1,lmaxapw,is,ia,ngk(ik,1),apwalm,evecfv(1,ist), &
               lmmaxapw,wfmt1(1,1,ist))
              done(ist)=.true.
            end if
! add to spinor wavefunction
            call zaxpy(n,zt1,wfmt1(1,1,ist),1,wfmt2(1,1,ispn),1)
          end if
        end do
      end do
! apply the L operator
      do ispn=1,nspinor
        do ir=1,nr
          call lopzflm(lmaxapw,wfmt2(1,ir,ispn),lmmaxapw,zlflm)
          do i=1,3
            wfmt3(:,ir,ispn,i)=zlflm(:,i)
          end do
        end do
      end do
! apply the S operator
      do ir=1,nr
! S_x
        wfmt4(:,ir,1,1)=0.5d0*wfmt2(:,ir,2)
        wfmt4(:,ir,2,1)=0.5d0*wfmt2(:,ir,1)
! S_y
        wfmt4(:,ir,1,2)=-0.5d0*zi*wfmt2(:,ir,2)
        wfmt4(:,ir,2,2)=0.5d0*zi*wfmt2(:,ir,1)
! S_z
        wfmt4(:,ir,1,3)=0.5d0*wfmt2(:,ir,1)
        wfmt4(:,ir,2,3)=-0.5d0*wfmt2(:,ir,2)
      end do
! compute the expectation values of L, S and J
      xl(:)=0.d0
      xs(:)=0.d0
      do i=1,3
        do ispn=1,nspinor
          zt1=zfmtinp(lmaxapw,nr,spr(1,is),lmmaxapw,wfmt2(1,1,ispn), &
           wfmt3(1,1,ispn,i))
          xl(i)=xl(i)+dble(zt1)
          zt1=zfmtinp(lmaxapw,nr,spr(1,is),lmmaxapw,wfmt2(1,1,ispn), &
           wfmt4(1,1,ispn,i))
          xs(i)=xs(i)+dble(zt1)
        end do
      end do
      xj(:)=xl(:)+xs(:)
      write(50,'("   L : ",3G18.10)') xl
      write(50,'("   S : ",3G18.10)') xs
      write(50,'("   J : ",3G18.10)') xj
    end do
  end do
! end loop over k-point and state list
end do
close(50)
write(*,*)
write(*,'("Info(writelsj):")')
write(*,'(" L, S and J expectation values written to LSJ.OUT")')
write(*,*)
deallocate(done,apwalm,evecfv,evecsv)
deallocate(wfmt1,wfmt2,wfmt3,wfmt4,zlflm)
return
end subroutine

