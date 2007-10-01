
! Copyright (C) 2002-2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine oepresk(ik,vnlcv,vnlvv,zvxmt,zvxir,zbxmt,zbxir,dvxmt,dvxir,dbxmt, &
 dbxir)
use modmain
implicit none
! arguments
integer, intent(in) :: ik
complex(8), intent(in) :: vnlcv(ncrmax,natmtot,nstsv,nkpt)
complex(8), intent(in) :: vnlvv(nstsv,nstsv,nkpt)
complex(8), intent(in) :: zvxmt(lmmaxvr,nrcmtmax,natmtot)
complex(8), intent(in) :: zvxir(ngrtot)
complex(8), intent(in) :: zbxmt(lmmaxvr,nrcmtmax,natmtot,ndmag)
complex(8), intent(in) :: zbxir(ngrtot,ndmag)
complex(8), intent(inout) :: dvxmt(lmmaxvr,nrcmtmax,natmtot)
complex(8), intent(inout) :: dvxir(ngrtot)
complex(8), intent(inout) :: dbxmt(lmmaxvr,nrcmtmax,natmtot,ndmag)
complex(8), intent(inout) :: dbxir(ngrtot,ndmag)
! local variables
integer is,ia,ias,ist1,ist2
integer nr,ic,m,idm
real(8) occ,de
complex(8) zde,zvnl,zrvx,zmbx,zt1,zt2
! allocatable arrays
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: evecfv(:,:)
complex(8), allocatable :: evecsv(:,:)
complex(8), allocatable :: wfmt(:,:,:,:,:)
complex(8), allocatable :: wfir(:,:,:)
complex(8), allocatable :: wfcr(:,:,:)
complex(8), allocatable :: zrhomt(:,:,:)
complex(8), allocatable :: zrhoir(:)
complex(8), allocatable :: zmagmt(:,:,:,:)
complex(8), allocatable :: zmagir(:,:)
complex(8), allocatable :: zvfmt(:,:,:)
complex(8), allocatable :: zfmt(:,:)
! external functions
complex(8) zfinp,zfmtinp
external zfinp,zfmtinp
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(evecfv(nmatmax,nstfv))
allocate(evecsv(nstsv,nstsv))
allocate(wfmt(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfir(ngrtot,nspinor,nstsv))
allocate(wfcr(lmmaxvr,nrcmtmax,2))
allocate(zrhomt(lmmaxvr,nrcmtmax,natmtot))
allocate(zrhoir(ngrtot))
if (spinpol) then
  allocate(zmagmt(lmmaxvr,nrcmtmax,natmtot,ndmag))
  allocate(zmagir(ngrtot,ndmag))
  allocate(zvfmt(lmmaxvr,nrcmtmax,ndmag))
  allocate(zfmt(lmmaxvr,nrcmtmax))
end if
if (spinpol) then
  occ=1.d0
else
  occ=2.d0
end if
! get the eigenvalues/vectors from file for input k-point
call getevalsv(vkl(1,ik),evalsv(1,ik))
call getevecfv(vkl(1,ik),vgkl(1,1,ik,1),evecfv)
call getevecsv(vkl(1,ik),evecsv)
! find the matching coefficients
call match(ngk(ik,1),gkc(1,ik,1),tpgkc(1,1,ik,1),sfacgk(1,1,ik,1),apwalm)
! calculate the wavefunctions for all states
call genwfsv(.false.,ngk(ik,1),igkig(1,ik,1),evalsv(1,ik),apwalm,evecfv, &
 evecsv,wfmt,wfir)
!-----------------------------------------------------------!
!     core-conduction overlap density and magnetisation     !
!-----------------------------------------------------------!
do is=1,nspecies
  nr=nrcmt(is)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    ic=0
    do ist1=1,spnst(is)
      if (spcore(ist1,is)) then
        do m=-spk(ist1,is),spk(ist1,is)-1
          ic=ic+1
! pass in m-1/2 to wavefcr
          call wavefcr(lradstp,is,ia,ist1,m,nrcmtmax,wfcr)
          do ist2=1,nstsv
            if (evalsv(ist2,ik).gt.efermi) then
              de=evalcr(ist1,ias)-evalsv(ist2,ik)
              zde=occ*wkpt(ik)/(de+zi*swidth)
! calculate the complex overlap density in the muffin-tin
              call vnlrhomt(is,wfcr(1,1,1),wfmt(1,1,ias,1,ist2),zrhomt(1,1,ias))
              if (spinpol) then
                call vnlrhomt(is,wfcr(1,1,2),wfmt(1,1,ias,2,ist2),zfmt)
                zrhomt(:,1:nr,ias)=zrhomt(:,1:nr,ias)+zfmt(:,1:nr)
              end if
              zvnl=conjg(vnlcv(ic,ias,ist2,ik))
              zrvx=zfmtinp(lmaxvr,nr,rcmt(1,is),lmmaxvr,zrhomt(1,1,ias), &
               zvxmt(1,1,ias))
              zt1=zvnl-zrvx
! spin-polarised case
              if (spinpol) then
                call oepmagmt(is,wfcr(1,1,1),wfcr(1,1,2),wfmt(1,1,ias,1,ist2), &
                 wfmt(1,1,ias,2,ist2),zvfmt)
! integral of magnetisation dot exchange field
                zmbx=0.d0
                do idm=1,ndmag
                  zmbx=zmbx+zfmtinp(lmaxvr,nr,rcmt(1,is),lmmaxvr, &
                   zvfmt(1,1,idm),zbxmt(1,1,ias,idm))
                end do
                zt1=zt1-zmbx
! end spin-polarised case
              end if
              zt2=zde*zt1
! residues for exchange potential and field
!$OMP CRITICAL
              dvxmt(:,1:nr,ias)=dvxmt(:,1:nr,ias)+zt2*zrhomt(:,1:nr,ias)
              do idm=1,ndmag
                dbxmt(:,1:nr,ias,idm)=dbxmt(:,1:nr,ias,idm) &
                 +zt2*zvfmt(:,1:nr,idm)
              end do
!$OMP END CRITICAL
! end loop over ist2
            end if
          end do
        end do
! end loop over ist1
      end if
    end do
! end loops over atoms and species
  end do
end do
!--------------------------------------------------------------!
!     valence-conduction overlap density and magnetisation     !
!--------------------------------------------------------------!
do ist1=1,nstsv
  if (evalsv(ist1,ik).lt.efermi) then
    do ist2=1,nstsv
      if (evalsv(ist2,ik).gt.efermi) then
! calculate the overlap density
        call vnlrho(wfmt(1,1,1,1,ist1),wfmt(1,1,1,1,ist2),wfir(1,1,ist1), &
         wfir(1,1,ist2),zrhomt,zrhoir)
        de=evalsv(ist1,ik)-evalsv(ist2,ik)
        zde=occ*wkpt(ik)/(de+zi*swidth)
        zvnl=conjg(vnlvv(ist1,ist2,ik))
        zrvx=zfinp(zrhomt,zvxmt,zrhoir,zvxir)
        zt1=zvnl-zrvx
        if (spinpol) then
          call oepmag(wfmt(1,1,1,1,ist1),wfmt(1,1,1,1,ist2),wfir(1,1,ist1), &
           wfir(1,1,ist2),zmagmt,zmagir)
! integral of magnetisation dot exchange field
          zmbx=0.d0
          do idm=1,ndmag
            zmbx=zmbx+zfinp(zmagmt(1,1,1,idm),zbxmt(1,1,1,idm),zmagir(1,idm), &
             zbxir(1,idm))
          end do
          zt1=zt1-zmbx
        end if
        zt2=zde*zt1
! residues for exchange potential and field
!$OMP CRITICAL
        do is=1,nspecies
          nr=nrcmt(is)
          do ia=1,natoms(is)
            ias=idxas(ia,is)
            dvxmt(:,1:nr,ias)=dvxmt(:,1:nr,ias)+zt2*zrhomt(:,1:nr,ias)
            do idm=1,ndmag
              dbxmt(:,1:nr,ias,idm)=dbxmt(:,1:nr,ias,idm) &
               +zt2*zmagmt(:,1:nr,ias,idm)
            end do
          end do
        end do
        dvxir(:)=dvxir(:)+zt2*zrhoir(:)
        do idm=1,ndmag
          dbxir(:,idm)=dbxir(:,idm)+zt2*zmagir(:,idm)
        end do
!$OMP END CRITICAL
! end loop over ist2
      end if
    end do
! end loop over ist1
  end if
end do
deallocate(apwalm,evecfv,evecsv)
deallocate(wfmt,wfir,wfcr,zrhomt,zrhoir)
if (spinpol) then
  deallocate(zmagmt,zmagir,zvfmt,zfmt)
end if
return
end subroutine

