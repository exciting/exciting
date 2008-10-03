
! Copyright (C) 2002-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: dos
! !INTERFACE:
subroutine dos
! !USES:
use modmain
! !DESCRIPTION:
!   Produces a total and partial density of states (DOS) for plotting. The total
!   DOS is written to the file {\tt TDOS.OUT} while the partial DOS is written
!   to the file {\tt PDOS\_Sss\_Aaaaa.OUT} for atom {\tt aaaa} of species
!   {\tt ss}. In the case of the partial DOS, each symmetrised
!   $(l,m)$-projection is written consecutively and separated by blank lines.
!   If the global variable {\tt lmirep} is {\tt .true.}, then the density matrix
!   from which the $(l,m)$-projections are obtained is first rotated into a
!   irreducible representation basis, i.e. one that block diagonalises all the
!   site symmetry matrices in the $Y_{lm}$ basis. Eigenvalues of a quasi-random
!   matrix in the $Y_{lm}$ basis which has been symmetrised with the site
!   symmetries are written to {\tt ELMIREP.OUT}. This allows for identification
!   of the irreducible representations of the site symmetries, for example $e_g$
!   or $t_{2g}$, by the degeneracies of the eigenvalues. In the plot, spin-up is
!   made positive and spin-down negative. See the routines {\tt gendmat} and
!   {\tt brzint}.
!
! !REVISION HISTORY:
!   Created January 2004 (JKD)
!EOP
!BOC
implicit none
! local variables
logical tsqaz
integer lmax,lmmax,l,m,lm
integer ispn,jspn,is,ia,ias
integer ik,nsk(3),ist,iw
real(8) dw,th,t1
real(8) v1(3),v2(3),v3(3)
complex(8) su2(2,2),dm1(2,2),dm2(2,2)
character(256) fname
! allocatable arrays
real(8), allocatable :: e(:,:,:)
real(8), allocatable :: f(:,:)
real(8), allocatable :: w(:)
real(8), allocatable :: g(:,:)
real(8), allocatable :: gp(:)
! low precision for band character array saves memory
real(4), allocatable :: bc(:,:,:,:,:)
real(8), allocatable :: elm(:,:)
complex(8), allocatable :: ulm(:,:,:)
complex(8), allocatable :: a(:,:)
complex(8), allocatable :: dmat(:,:,:,:,:)
complex(8), allocatable :: sdmat(:,:,:,:)
complex(8), allocatable :: apwalm(:,:,:,:,:)
complex(8), allocatable :: evecfv(:,:,:)
complex(8), allocatable :: evecsv(:,:)
! initialise universal variables
call init0
call init1
lmax=min(3,lmaxapw)
lmmax=(lmax+1)**2
! allocate local arrays
allocate(e(nstsv,nkpt,nspinor))
allocate(f(nstsv,nkpt))
allocate(w(nwdos))
allocate(g(nwdos,nspinor))
allocate(gp(nwdos))
allocate(bc(lmmax,nspinor,natmtot,nstsv,nkpt))
if (lmirep) then
  allocate(elm(lmmax,natmtot))
  allocate(ulm(lmmax,lmmax,natmtot))
  allocate(a(lmmax,lmmax))
end if
allocate(dmat(lmmax,lmmax,nspinor,nspinor,nstsv))
allocate(sdmat(nspinor,nspinor,nstsv,nkpt))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
allocate(evecfv(nmatmax,nstfv,nspnfv))
allocate(evecsv(nstsv,nstsv))
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
! generate unitary matrices which convert the (l,m) basis into the irreducible
! representation basis of the symmetry group at each atomic site
if (lmirep) then
  call genlmirep(lmax,lmmax,elm,ulm)
end if
! compute the SU(2) operator used for rotating the density matrix to the
! desired spin-quantisation axis
v1(:)=sqados(:)
t1=sqrt(v1(1)**2+v1(2)**2+v1(3)**2)
if (t1.le.epslat) then
  write(*,*)
  write(*,'("Error(dos): spin-quantisation axis (sqados) has zero length")')
  write(*,*)
  stop
end if
v1(:)=v1(:)/t1
if (v1(3).ge.1.d0-epslat) then
  tsqaz=.true.
else
  tsqaz=.false.
  v2(1:2)=0.d0
  v2(3)=1.d0
  call r3cross(v1,v2,v3)
! note that the spin-quantisation axis is rotated, so the density matrix should
! be rotated in the opposite direction
  th=-acos(v1(3))
  call axangsu2(v3,th,su2)
end if
! loop over k-points
do ik=1,nkpt
! get the eigenvalues/vectors from file
  call getevalsv(vkl(:,ik),evalsv(:,ik))
  call getevecfv(vkl(:,ik),vgkl(:,:,:,ik),evecfv)
  call getevecsv(vkl(:,ik),evecsv)
! find the matching coefficients
  do ispn=1,nspnfv
    call match(ngk(ispn,ik),gkc(:,ispn,ik),tpgkc(:,:,ispn,ik), &
     sfacgk(:,:,ispn,ik),apwalm(:,:,:,:,ispn))
  end do
  do is=1,nspecies
    do ia=1,natoms(is)
      ias=idxas(ia,is)
! generate the density matrix
      call gendmat(.false.,.false.,0,lmax,is,ia,ngk(:,ik),apwalm,evecfv, &
       evecsv,lmmax,dmat)
! convert (l,m) part to an irreducible representation if required
      if (lmirep) then
        do ist=1,nstsv
          do ispn=1,nspinor
            do jspn=1,nspinor
              call zgemm('N','N',lmmax,lmmax,lmmax,zone,ulm(:,:,ias),lmmax, &
               dmat(:,:,ispn,jspn,ist),lmmax,zzero,a,lmmax)
              call zgemm('N','C',lmmax,lmmax,lmmax,zone,a,lmmax,ulm(:,:,ias), &
               lmmax,zzero,dmat(:,:,ispn,jspn,ist),lmmax)
            end do
          end do
        end do
      end if
! spin rotate the density matrices to desired spin-quantisation axis
      if (spinpol.and.(.not.tsqaz)) then
        do ist=1,nstsv
          do lm=1,lmmax
            dm1(:,:)=dmat(lm,lm,:,:,ist)
            call z2mm(su2,dm1,dm2)
            call z2mmct(dm2,su2,dm1)
            dmat(lm,lm,:,:,ist)=dm1(:,:)
          end do
        end do
      end if
! determine the band characters from the density matrix
      do ist=1,nstsv
        do ispn=1,nspinor
          do lm=1,lmmax
            t1=dble(dmat(lm,lm,ispn,ispn,ist))
            bc(lm,ispn,ias,ist,ik)=real(t1)
          end do
        end do
      end do
    end do
  end do
! compute the spin density matrices of the second-variational states
  call gensdmat(evecsv,sdmat(:,:,:,ik))
! spin rotate the density matrices to desired spin-quantisation axis
  if (spinpol.and.(.not.tsqaz)) then
    do ist=1,nstsv
      call z2mm(su2,sdmat(:,:,ist,ik),dm1)
      call z2mmct(dm1,su2,sdmat(:,:,ist,ik))
    end do
  end if
end do
! generate energy grid
dw=(wdos(2)-wdos(1))/dble(nwdos)
do iw=1,nwdos
  w(iw)=dw*dble(iw-1)+wdos(1)
end do
! number of subdivisions used for interpolation
nsk(:)=max(ngrdos/ngridk(:),1)
!--------------------------!
!     output total DOS     !
!--------------------------!
open(50,file='TDOS.OUT',action='WRITE',form='FORMATTED')
do ispn=1,nspinor
  if (ispn.eq.1) then
    t1=1.d0
  else
    t1=-1.d0
  end if
  do ik=1,nkpt
    do ist=1,nstsv
! subtract the Fermi energy
      e(ist,ik,ispn)=evalsv(ist,ik)-efermi
! correction for scissors operator
      if (e(ist,ik,ispn).gt.0.d0) e(ist,ik,ispn)=e(ist,ik,ispn)+scissor
! use diagonal of spin density matrix for weight
      f(ist,ik)=dble(sdmat(ispn,ispn,ist,ik))
    end do
  end do
  call brzint(nsmdos,ngridk,nsk,ikmap,nwdos,wdos,nstsv,nstsv,e(:,:,ispn),f, &
   g(:,ispn))
! multiply by the maximum occupancy (spin-polarised: 1, unpolarised: 2)
  g(:,ispn)=occmax*g(:,ispn)
  do iw=1,nwdos
    write(50,'(2G18.10)') w(iw),t1*g(iw,ispn)
  end do
  write(50,'("     ")')
end do
close(50)
!----------------------------!
!     output partial DOS     !
!----------------------------!
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    write(fname,'("PDOS_S",I2.2,"_A",I4.4,".OUT")') is,ia
    open(50,file=trim(fname),action='WRITE',form='FORMATTED')
    do ispn=1,nspinor
      if (ispn.eq.1) then
        t1=1.d0
      else
        t1=-1.d0
      end if
      do l=0,lmax
        do m=-l,l
          lm=idxlm(l,m)
          do ik=1,nkpt
            do ist=1,nstsv
              f(ist,ik)=bc(lm,ispn,ias,ist,ik)
            end do
          end do
          call brzint(nsmdos,ngridk,nsk,ikmap,nwdos,wdos,nstsv,nstsv, &
           e(:,:,ispn),f,gp)
          gp(:)=occmax*gp(:)
          do iw=1,nwdos
            write(50,'(2G18.10)') w(iw),t1*gp(iw)
! interstitial DOS
            g(iw,ispn)=g(iw,ispn)-gp(iw)
          end do
          write(50,'("     ")')
        end do
      end do
    end do
    close(50)
  end do
end do
if (lmirep) then
  open(50,file='ELMIREP.OUT',action='WRITE',form='FORMATTED')
  do is=1,nspecies
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      write(50,*)
      write(50,'("Species : ",I4," (",A,"), atom : ",I4)') is, &
       trim(spsymb(is)),ia
      do l=0,lmax
        do m=-l,l
          lm=idxlm(l,m)
          write(50,'(" l = ",I2,", m = ",I2,", lm= ",I3," : ",G18.10)') l,m, &
           lm,elm(lm,ias)
        end do
      end do
    end do
  end do
  close(50)
end if
!---------------------------------!
!     output interstitial DOS     !
!---------------------------------!
open(50,file='IDOS.OUT',action='WRITE',form='FORMATTED')
do ispn=1,nspinor
  if (ispn.eq.1) then
    t1=1.d0
  else
    t1=-1.d0
  end if
  do iw=1,nwdos
    write(50,'(2G18.10)') w(iw),t1*g(iw,ispn)
  end do
end do
close(50)
write(*,*)
write(*,'("Info(dos):")')
write(*,'(" Total density of states written to TDOS.OUT")')
write(*,*)
write(*,'(" Partial density of states written to PDOS_Sss_Aaaaa.OUT")')
write(*,'(" for all species and atoms")')
if (lmirep) then
  write(*,*)
  write(*,'(" Eigenvalues of a random matrix in the (l,m) basis symmetrised")')
  write(*,'(" with the site symmetries written to ELMIREP.OUT for all")')
  write(*,'(" species and atoms. Degenerate eigenvalues correspond to")')
  write(*,'(" irreducible representations of each site symmetry group")')
end if
write(*,*)
write(*,'(" Interstitial density of states written to IDOS.OUT")')
write(*,*)
write(*,'(" Fermi energy is at zero in plot")')
write(*,*)
write(*,'(" DOS units are states/Hartree/unit cell")')
write(*,*)
deallocate(e,f,w,g,gp,bc)
if (lmirep) deallocate(elm,ulm,a)
deallocate(dmat,sdmat,apwalm,evecfv,evecsv)
return
end subroutine
!EOC
