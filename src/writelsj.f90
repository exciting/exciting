
! Copyright (C) 2002-2007 J. K. Dewhurst, S. Sharma, C. Ambrosch-Draxl and
! F. Cricchio. This file is distributed under the terms of the GNU General
! Public License. See the file COPYING for license details.

subroutine writelsj
use modmain
implicit none
! local variables
integer kst,ik,ist,lm
integer ispn,is,ia,ias
real(8) xl(3),xs(3),t1
! allocatable arrays
complex(8), allocatable :: apwalm(:,:,:,:,:)
complex(8), allocatable :: evecfv(:,:,:)
complex(8), allocatable :: evecsv(:,:)
complex(8), allocatable :: dmat1(:,:,:,:,:)
complex(8), allocatable :: dmat2(:,:,:,:,:)
complex(8), allocatable :: zlflm(:,:)
! initialise universal variables
call init0
call init1
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
allocate(evecfv(nmatmax,nstfv,nspnfv))
allocate(evecsv(nstsv,nstsv))
allocate(dmat1(lmmaxapw,lmmaxapw,nspinor,nspinor,natmtot))
allocate(dmat2(lmmaxapw,lmmaxapw,nspinor,nspinor,nstsv))
allocate(zlflm(lmmaxapw,3))
! read density and potentials from file
call readstate
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
if (task.eq.15) then
! compute total L, S and J
  dmat1(:,:,:,:,:)=0.d0
  do ik=1,nkpt
! get the eigenvectors and occupancies from file
    call getevecfv(vkl(:,ik),vgkl(:,:,:,ik),evecfv)
    call getevecsv(vkl(:,ik),evecsv)
    call getoccsv(vkl(:,ik),occsv(:,ik))
! find the matching coefficients
    do ispn=1,nspnfv
      call match(ngk(ispn,ik),gkc(:,ispn,ik),tpgkc(:,:,ispn,ik), &
       sfacgk(:,:,ispn,ik),apwalm(:,:,:,:,ispn))
    end do
! loop over species and atoms
    do is=1,nspecies
      do ia=1,natoms(is)
        ias=idxas(ia,is)
! generate the density matrix
        call gendmat(.false.,.false.,0,lmaxapw,is,ia,ngk(:,ik),apwalm,evecfv, &
         evecsv,lmmaxapw,dmat2)
        do ist=1,nstsv
          t1=wkpt(ik)*occsv(ist,ik)
          dmat1(:,:,:,:,ias)=dmat1(:,:,:,:,ias)+t1*dmat2(:,:,:,:,ist)
        end do
      end do
    end do
! end loop over k-points
  end do
! symmetrise the density matrix
  call symdmat(lmaxapw,lmmaxapw,dmat1)
  open(50,file='LSJ.OUT',action='WRITE',form='FORMATTED')
  write(50,*)
  write(50,'("Expectation values are computed only over the muffin-tin")')
! loop over species and atoms
  do is=1,nspecies
    write(50,*)
    write(50,'("Species : ",I4," (",A,")")') is,trim(spsymb(is))
    do ia=1,natoms(is)
      ias=idxas(ia,is)
! compute tr(LD)
      xl(:)=0.d0
      do ispn=1,nspinor
        do lm=1,lmmaxapw
          call lopzflm(lmaxapw,dmat1(:,lm,ispn,ispn,ias),lmmaxapw,zlflm)
          xl(:)=xl(:)+dble(zlflm(lm,:))
        end do
      end do
! compute tr(sigma D)
      xs(:)=0.d0
      if (spinpol) then
        do lm=1,lmmaxapw
          xs(1)=xs(1)+dble(dmat1(lm,lm,2,1,ias)+dmat1(lm,lm,1,2,ias))
          xs(2)=xs(2)+dble(-zi*dmat1(lm,lm,2,1,ias)+zi*dmat1(lm,lm,1,2,ias))
          xs(3)=xs(3)+dble(dmat1(lm,lm,1,1,ias)-dmat1(lm,lm,2,2,ias))
        end do
      end if
! S = 1/2 sigma
      xs(:)=0.5d0*xs(:)
      write(50,'(" atom : ",I4)') ia
      write(50,'("  L : ",3G18.10)') xl(:)
      write(50,'("  S : ",3G18.10)') xs(:)
      write(50,'("  J : ",3G18.10)') xl(:)+xs(:)
! end loop over atoms and species
    end do
  end do
  close(50)
  write(*,*)
  write(*,'("Info(writelsj):")')
  write(*,'(" total L, S and J expectation values written to LSJ.OUT")')
  write(*,*)
else
! compute L, S and J for all states in kstlist
  open(50,file='LSJ_KST.OUT',action='WRITE',form='FORMATTED')
  write(50,*)
  write(50,'("Expectation values are computed only over the muffin-tin")')
  do kst=1,nkstlist
    ik=kstlist(1,kst)
    ist=kstlist(2,kst)
    if ((ik.le.0).or.(ik.gt.nkpt)) then
      write(*,*)
      write(*,'("Error(writelsj): k-point out of range : ",I8)') ik
      write(*,*)
      stop
    end if
    if ((ist.le.0).or.(ist.gt.nstsv)) then
      write(*,*)
      write(*,'("Error(writelsj): state out of range : ",I8)') ist
      write(*,*)
      stop
    end if
! get the eigenvectors and occupancies from file
    call getevecfv(vkl(:,ik),vgkl(:,:,:,ik),evecfv)
    call getevecsv(vkl(:,ik),evecsv)
    call getoccsv(vkl(:,ik),occsv(:,ik))
! find the matching coefficients
    do ispn=1,nspnfv
      call match(ngk(ispn,ik),gkc(:,ispn,ik),tpgkc(:,:,ispn,ik), &
       sfacgk(:,:,ispn,ik),apwalm(:,:,:,:,ispn))
    end do
! loop over species and atoms
    do is=1,nspecies
      do ia=1,natoms(is)
        ias=idxas(ia,is)
! generate the density matrix
        call gendmat(.false.,.false.,0,lmaxapw,is,ia,ngk(:,ik),apwalm,evecfv, &
         evecsv,lmmaxapw,dmat2)
! compute tr(LD)
        xl(:)=0.d0
        do ispn=1,nspinor
          do lm=1,lmmaxapw
            call lopzflm(lmaxapw,dmat2(:,lm,ispn,ispn,ist),lmmaxapw,zlflm)
            xl(:)=xl(:)+dble(zlflm(lm,:))
          end do
        end do
! compute tr(sigma D)
        xs(:)=0.d0
        if (spinpol) then
          do lm=1,lmmaxapw
            xs(1)=xs(1)+dble(dmat2(lm,lm,2,1,ist)+dmat2(lm,lm,1,2,ist))
            xs(2)=xs(2)+dble(-zi*dmat2(lm,lm,2,1,ist)+zi*dmat2(lm,lm,1,2,ist))
            xs(3)=xs(3)+dble(dmat2(lm,lm,1,1,ist)-dmat2(lm,lm,2,2,ist))
          end do
        else
          xs(3)=1.d0
        end if
! S = 1/2 sigma
        xs(:)=0.5d0*xs(:)
        write(50,*)
        write(50,'("k-point : ",I6,3G18.10)') ik,vkl(:,ik)
        write(50,'("state : ",I6)') ist
        write(50,'("species : ",I4," (",A,"), atom : ",I4)') is, &
         trim(spsymb(is)),ia
        write(50,'(" L : ",3G18.10)') xl(:)
        write(50,'(" S : ",3G18.10)') xs(:)
        write(50,'(" J : ",3G18.10)') xl(:)+xs(:)
      end do
    end do
  end do
  close(50)
  write(*,*)
  write(*,'("Info(writelsj):")')
  write(*,'(" L, S and J expectation values for each k-point and state")')
  write(*,'("  in kstlist written to LSJ_KST.OUT")')
  write(*,*)
end if
deallocate(apwalm,evecfv,evecsv,dmat1,dmat2,zlflm)
return
end subroutine

