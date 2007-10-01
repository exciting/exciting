
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
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
!   Eigenvalues of a matrix in the $Y_{lm}$ basis which has been symmetrised
!   with the site symmetries are written to {\tt ELMSYM.OUT}. This allows for
!   identification of the irreducible representations of the site symmetries,
!   for example $e_g$ or $t_{2g}$. Spin-up is made positive and spin-down
!   negative for the partial DOS in both the collinear and non-collinear cases
!   and for the total DOS in the collinear case. See the routines {\tt bandchar}
!   and {\tt brzint}.
!
! !REVISION HISTORY:
!   Created January 2004 (JKD)
!EOP
!BOC
implicit none
! local variables
integer lmax,lmmax,l,m,lm,nsk(3)
integer ik,ispn,is,ia,ias,ist,iw,i
real(8) t1
character(256) fname
! allocatable arrays
real(8), allocatable :: bndchr(:,:,:,:,:)
real(8), allocatable :: elmsym(:,:)
real(8), allocatable :: e(:,:,:)
real(8), allocatable :: f(:,:)
real(8), allocatable :: w(:)
real(8), allocatable :: g(:,:)
real(8), allocatable :: gp(:)
complex(8), allocatable :: evecfv(:,:,:)
complex(8), allocatable :: evecsv(:,:)
! initialise universal variables
call init0
call init1
allocate(evecfv(nmatmax,nstfv,nspnfv))
allocate(evecsv(nstsv,nstsv))
! allocate the band character array
lmax=min(3,lmaxapw)
lmmax=(lmax+1)**2
allocate(bndchr(lmmax,natmtot,nspinor,nstsv,nkpt))
allocate(elmsym(lmmax,natmtot))
! allocate local arrays
allocate(e(nstsv,nkpt,nspinor))
allocate(f(nstsv,nkpt))
allocate(w(nwdos))
allocate(g(nwdos,nspinor))
allocate(gp(nwdos))
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
do ik=1,nkpt
! get the eigenvalues/vectors from file
  call getevalsv(vkl(1,ik),evalsv(1,ik))
  call getevecfv(vkl(1,ik),vgkl(1,1,ik,1),evecfv)
  call getevecsv(vkl(1,ik),evecsv)
! compute the band character (appromximate for spin-spirals)
  call bandchar(lmax,ik,evecfv,evecsv,lmmax,bndchr(1,1,1,1,ik),elmsym)
! compute the spin characters
  call spinchar(ik,evecsv)
end do
! generate energy grid
t1=(wdos(2)-wdos(1))/dble(nwdos)
do iw=1,nwdos
  w(iw)=t1*dble(iw-1)+wdos(1)
end do
! number of subdivisions used for interpolation
do i=1,3
  nsk(i)=max(ngrdos/ngridk(i),1)
end do
! output total (spin-projected) DOS
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
! use spin character for weight
      f(ist,ik)=spnchr(ispn,ist,ik)
    end do
  end do
  call brzint(nsmdos,ngridk,nsk,ikmap,nwdos,wdos,nstsv,nstsv,e(1,1,ispn),f, &
   g(1,ispn))
  do iw=1,nwdos
    write(50,'(2G18.10)') w(iw),t1*g(iw,ispn)
  end do
  write(50,'("     ")')
end do
close(50)
! output partial DOS
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
              f(ist,ik)=bndchr(lm,ias,ispn,ist,ik)
            end do
          end do
          call brzint(nsmdos,ngridk,nsk,ikmap,nwdos,wdos,nstsv,nstsv, &
           e(1,1,ispn),f,gp)
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
! output eigenvalues of (l,m)-projection operator
open(50,file='ELMSYM.OUT',action='WRITE',form='FORMATTED')
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    write(50,*)
    write(50,'("Species : ",I4," (",A,"), atom : ",I4)') is,trim(spsymb(is)),ia
    do l=0,lmax
      do m=-l,l
        lm=idxlm(l,m)
        write(50,'(" l =",I2,", m =",I2,", lm=",I3," : ",G18.10)') l,m,lm, &
         elmsym(lm,ias)
      end do
    end do
  end do
end do
close(50)
! output interstitial DOS
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
write(*,'(" total density of states written to TDOS.OUT")')
write(*,'(" partial density of states written to PDOS_Sss_Aaaaa.OUT")')
write(*,'("  for all species and atoms")')
write(*,'(" eigenvalues of a matrix in the Y_lm basis which has been")')
write(*,'("  symmetrised with the site symmetries written to ELMSYM.OUT")')
write(*,'(" interstitial density of states written to IDOS.OUT")')
write(*,*)
write(*,'(" (DOS units are states/Hartree/spin/unit cell)")')
write(*,*)
deallocate(bndchr,elmsym,e,f,w,g,gp)
deallocate(evecfv,evecsv)
return
end subroutine
!EOC
