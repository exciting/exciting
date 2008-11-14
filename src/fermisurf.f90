
! Copyright (C) 2002-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine fermisurf
use modmain
implicit none
! local variables
integer ik,jk,ist
integer ist0,ist1,nst
real(8) prd1,prd2
! allocatable arrays
real(8), allocatable :: evalfv(:,:)
complex(8), allocatable :: evecfv(:,:,:)
complex(8), allocatable :: evecsv(:,:)
! initialise universal variables
call init0
call init1
! read density and potentials from file
call readstate
! read Fermi energy from file
call readfermi
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
! compute the overlap radial integrals
call olprad
! compute the Hamiltonian radial integrals
call hmlrad
! begin parallel loop over reduced k-points set
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(evalfv,evecfv,evecsv) 
!$OMP DO
do ik=1,nkpt
  allocate(evalfv(nstfv,nspnfv))
  allocate(evecfv(nmatmax,nstfv,nspnfv))
  allocate(evecsv(nstsv,nstsv))
!$OMP CRITICAL
  write(*,'("Info(fermisurf): ",I6," of ",I6," k-points")') ik,nkpt
!$OMP END CRITICAL
! solve the first- and second-variational secular equations
  call seceqn(ik,evalfv,evecfv,evecsv)
  deallocate(evalfv,evecfv,evecsv)
! end loop over reduced k-points set
end do
!$OMP END DO
!$OMP END PARALLEL
if (ndmag.eq.1) then
! special case of collinear magnetism
  open(50,file='FERMISURF_UP.OUT',action='WRITE',form='FORMATTED')
  open(51,file='FERMISURF_DN.OUT',action='WRITE',form='FORMATTED')
  if (task.eq.100) then
! write product of eigenstates minus the Fermi energy
    write(50,'(3I6," : grid size")') np3d(:)
    write(51,'(3I6," : grid size")') np3d(:)
    do ik=1,nkptnr
      jk=ikmap(ivknr(1,ik),ivknr(2,ik),ivknr(3,ik))
      prd1=1.d0
      prd2=1.d0
      do ist=1,nstfv
        prd1=prd1*(evalsv(ist,jk)-efermi)
        prd2=prd2*(evalsv(nstfv+ist,jk)-efermi)
      end do
      write(50,'(4G18.10)') vkcnr(:,ik),prd1
      write(51,'(4G18.10)') vkcnr(:,ik),prd2
    end do
  else
! write the eigenvalues minus the Fermi energy separately
    ist=nstfv-nempty
    ist0=max(ist-nstfsp/2,1)
    ist1=min(ist+nstfsp/2,nstfv)
    nst=ist1-ist0+1
    write(50,'(4I6," : grid size, number of states")') np3d(:),nst
    write(51,'(4I6," : grid size, number of states")') np3d(:),nst
    do ik=1,nkptnr
      jk=ikmap(ivknr(1,ik),ivknr(2,ik),ivknr(3,ik))
      write(50,'(G18.10)',advance='NO') vkcnr(:,ik)
      write(51,'(G18.10)',advance='NO') vkcnr(:,ik)
      do ist=ist0,ist1
        write(50,'(F14.8)',advance='NO') evalsv(ist,jk)-efermi
        write(51,'(F14.8)',advance='NO') evalsv(nstfv+ist,jk)-efermi
      end do
      write(50,*)
      write(51,*)
    end do
  end if
  close(50)
  close(51)
else
! spin-unpolarised and non-collinear cases
  open(50,file='FERMISURF.OUT',action='WRITE',form='FORMATTED')
  if (task.eq.100) then
! write product of eigenstates minus the Fermi energy
    write(50,'(3I6," : grid size")') np3d(:)
    do ik=1,nkptnr
      jk=ikmap(ivknr(1,ik),ivknr(2,ik),ivknr(3,ik))
      prd1=1.d0
      do ist=1,nstsv
        prd1=prd1*(evalsv(ist,jk)-efermi)
      end do
      write(50,'(4G18.10)') vkcnr(:,ik),prd1
    end do
  else
! write the eigenvalues minus the Fermi energy separately
    ist=(nstfv-nempty)*nspinor
    ist0=max(ist-nstfsp/2,1)
    ist1=min(ist+nstfsp/2,nstsv)
    nst=ist1-ist0+1
    write(50,'(4I6," : grid size, number of states")') np3d(:),nst
    do ik=1,nkptnr
      jk=ikmap(ivknr(1,ik),ivknr(2,ik),ivknr(3,ik))
      write(50,'(3G18.10)',advance='NO') vkcnr(:,ik)
      do ist=ist0,ist1
        write(50,'(F14.8)',advance='NO') evalsv(ist,jk)-efermi
      end do
      write(50,*)
    end do
  end if
  close(50)
end if
write(*,*)
write(*,'("Info(fermisurf):")')
if (ndmag.eq.1) then
  write(*,'(" 3D Fermi surface data written to FERMISURF_UP.OUT and&
   & FERMISURF_DN.OUT")')
else
  write(*,'(" 3D Fermi surface data written to FERMISURF.OUT")')
end if
if (task.eq.100) then
  write(*,'(" in terms of the product of eigenvalues minus the Fermi energy")')
else
  write(*,'(" in terms of separate eigenvalues minus the Fermi energy")')
end if
write(*,*)
return
end subroutine

