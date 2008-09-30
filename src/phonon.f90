
! Copyright (C) 2002-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine phonon
use modmain
implicit none
! local variables
integer is,js,ia,ja,ka,jas,kas
integer iq,ip,jp,nph,iph,i
real(8) dph,a,b,t1
real(8) ftp(3,maxatoms,maxspecies)
complex(8) zt1,zt2
complex(8) dyn(3,maxatoms,maxspecies)
! allocatable arrays
real(8), allocatable :: veffmtp(:,:,:)
real(8), allocatable :: veffirp(:)
complex(8), allocatable :: dveffmt(:,:,:)
complex(8), allocatable :: dveffir(:)
!------------------------!
!     initialisation     !
!------------------------!
! require forces
tforce=.true.
! no primitive cell determination
primcell=.false.
! initialise universal variables
call init0
! initialise q-point dependent variables
call init2
! read original effective potential from file and store in global arrays
call readstate
if (allocated(veffmt0)) deallocate(veffmt0)
allocate(veffmt0(lmmaxvr,nrmtmax,natmtot))
if (allocated(veffir0)) deallocate(veffir0)
allocate(veffir0(ngrtot))
veffmt0(:,:,:)=veffmt(:,:,:)
veffir0(:)=veffir(:)
! allocate local arrays
allocate(dveffmt(lmmaxvr,nrcmtmax,natmtot))
allocate(dveffir(ngrtot))
! switch off automatic determination of muffin-tin radii
autormt=.false.
! no shifting of atomic basis allowed
tshift=.false.
! determine k-point grid size from radkpt
autokpt=.true.
! store original parameters
natoms0(1:nspecies)=natoms(1:nspecies)
natmtot0=natmtot
avec0(:,:)=avec(:,:)
ainv0(:,:)=ainv(:,:)
atposc0(:,:,:)=0.d0
do is=1,nspecies
  do ia=1,natoms(is)
    atposc0(:,ia,is)=atposc(:,ia,is)
  end do
end do
ngrid0(:)=ngrid(:)
ngrtot0=ngrtot
!---------------------------------------!
!     compute dynamical matrix rows     !
!---------------------------------------!
10 continue
natoms(1:nspecies)=natoms0(1:nspecies)
! find a dynamical matrix to calculate
call dyntask(80,iq,is,ia,ip)
! phonon dry run
if (task.eq.201) goto 10
! check to see if mass is considered infinite
if (spmass(is).le.0.d0) then
  do ip=1,3
    do js=1,nspecies
      do ja=1,natoms0(js)
        do jp=1,3
          write(80,'(2G18.10," : is = ",I4,", ia = ",I4,", ip = ",I4)') 0.d0, &
           0.d0,js,ja,jp
        end do
      end do
    end do
  end do
  close(80)
  goto 10
end if
task=200
nph=1
if ((ivq(1,iq).eq.0).and.(ivq(2,iq).eq.0).and.(ivq(3,iq).eq.0)) nph=0
dyn(:,:,:)=0.d0
dveffmt(:,:,:)=0.d0
dveffir(:)=0.d0
! loop over phases (cos and sin displacements)
do iph=0,nph
! restore input values
  natoms(1:nspecies)=natoms0(1:nspecies)
  avec(:,:)=avec0(:,:)
  atposc(:,:,:)=atposc0(:,:,:)
! generate the supercell
  call phcell(iph,deltaph,iq,is,ia,ip)
! run the ground-state calculation
  call gndstate
! store the total force for the first displacement
  do js=1,nspecies
    do ja=1,natoms(js)
      jas=idxas(ja,js)
      ftp(:,ja,js)=forcetot(:,jas)
    end do
  end do
! store the effective potential for the first displacement
  allocate(veffmtp(lmmaxvr,nrmtmax,natmtot))
  allocate(veffirp(ngrtot))
  veffmtp(:,:,:)=veffmt(:,:,:)
  veffirp(:)=veffir(:)
! restore input values
  natoms(1:nspecies)=natoms0(1:nspecies)
  avec(:,:)=avec0(:,:)
  atposc(:,:,:)=atposc0(:,:,:)
! generate the supercell again with twice the displacement
  dph=deltaph+deltaph
  call phcell(iph,dph,iq,is,ia,ip)
! run the ground-state calculation again starting from the previous density
  task=1
  call gndstate
! compute the complex perturbing effective potential with implicit q-phase
  call phdveff(iph,iq,veffmtp,veffirp,dveffmt,dveffir)
  deallocate(veffmtp,veffirp)
! Fourier transform the force differences to obtain the dynamical matrix
  zt1=1.d0/(dble(nphcell)*deltaph)
! multiply by i for sin-like displacement
  if (iph.eq.1) zt1=zt1*zi
  kas=0
  do js=1,nspecies
    ka=0
    do ja=1,natoms0(js)
      do i=1,nphcell
        ka=ka+1
        kas=kas+1
        t1=-dot_product(vqc(:,iq),vphcell(:,i))
        zt2=zt1*cmplx(cos(t1),sin(t1),8)
        do jp=1,3
          t1=-(forcetot(jp,kas)-ftp(jp,ka,js))
          dyn(jp,ja,js)=dyn(jp,ja,js)+zt2*t1
        end do
      end do
    end do
  end do
end do
! write dynamical matrix row to file
do js=1,nspecies
  do ja=1,natoms0(js)
    do jp=1,3
      a=dble(dyn(jp,ja,js))
      b=aimag(dyn(jp,ja,js))
      if (abs(a).lt.1.d-12) a=0.d0
      if (abs(b).lt.1.d-12) b=0.d0
      write(80,'(2G18.10," : is = ",I4,", ia = ",I4,", ip = ",I4)') a,b,js,ja,jp
    end do
  end do
end do
close(80)
! write the complex perturbing effective potential to file
call writedveff(iq,is,ia,ip,dveffmt,dveffir)
! delete the non-essential files
call phdelete
goto 10
end subroutine

