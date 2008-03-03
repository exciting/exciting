
! Copyright (C) 2002-2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine oepresk(ik,vnlcv,vnlvv,dvxmt,dvxir,dbxmt,dbxir)
use modmain
implicit none
! arguments
integer, intent(in) :: ik
complex(8), intent(in) :: vnlcv(ncrmax,natmtot,nstsv,nkpt)
complex(8), intent(in) :: vnlvv(nstsv,nstsv,nkpt)
real(8), intent(inout) :: dvxmt(lmmaxvr,nrcmtmax,natmtot)
real(8), intent(inout) :: dvxir(ngrtot)
real(8), intent(inout) :: dbxmt(lmmaxvr,nrcmtmax,natmtot,ndmag)
real(8), intent(inout) :: dbxir(ngrtot,ndmag)
! local variables
integer is,ia,ias,ist,jst
integer nr,ic,m,idm
real(8) de
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
    do ist=1,spnst(is)
      if (spcore(ist,is)) then
        do m=-spk(ist,is),spk(ist,is)-1
          ic=ic+1
! pass in m-1/2 to wavefcr
          call wavefcr(lradstp,is,ia,ist,m,nrcmtmax,wfcr)
          do jst=1,nstsv
            if (evalsv(jst,ik).gt.efermi) then
              de=evalcr(ist,ias)-evalsv(jst,ik)
              zde=occmax*wkpt(ik)/(de+zi*swidth)
! calculate the complex overlap density in the muffin-tin
              call vnlrhomt(.false.,is,wfcr(1,1,1),wfmt(1,1,ias,1,jst), &
               zrhomt(1,1,ias))
              if (spinpol) then
                call vnlrhomt(.false.,is,wfcr(1,1,2),wfmt(1,1,ias,2,jst),zfmt)
                zrhomt(:,1:nr,ias)=zrhomt(:,1:nr,ias)+zfmt(:,1:nr)
              end if
              zvnl=conjg(vnlcv(ic,ias,jst,ik))
              zrvx=zfmtinp(.false.,lmaxvr,nr,rcmt(1,is),lmmaxvr, &
               zrhomt(1,1,ias),zvxmt(1,1,ias))
              zt1=zvnl-zrvx
! spin-polarised case
              if (spinpol) then
                call oepmagmt(.false.,is,wfcr(1,1,1),wfcr(1,1,2), &
                 wfmt(1,1,ias,1,jst),wfmt(1,1,ias,2,jst),zvfmt)
! integral of magnetisation dot exchange field
                zmbx=0.d0
                do idm=1,ndmag
                  zmbx=zmbx+zfmtinp(.false.,lmaxvr,nr,rcmt(1,is),lmmaxvr, &
                   zvfmt(1,1,idm),zbxmt(1,1,ias,idm))
                end do
                zt1=zt1-zmbx
! end spin-polarised case
              end if
              zt2=zde*zt1
! residues for exchange potential and field
!$OMP CRITICAL
              dvxmt(:,1:nr,ias)=dvxmt(:,1:nr,ias)+dble(zt2*zrhomt(:,1:nr,ias))
              do idm=1,ndmag
                dbxmt(:,1:nr,ias,idm)=dbxmt(:,1:nr,ias,idm) &
                 +dble(zt2*zvfmt(:,1:nr,idm))
              end do
!$OMP END CRITICAL
! end loop over jst
            end if
          end do
        end do
! end loop over ist
      end if
    end do
! end loops over atoms and species
  end do
end do
!--------------------------------------------------------------!
!     valence-conduction overlap density and magnetisation     !
!--------------------------------------------------------------!
do ist=1,nstsv
  if (evalsv(ist,ik).lt.efermi) then
    do jst=1,nstsv
      if (evalsv(jst,ik).gt.efermi) then
! calculate the overlap density
        call vnlrho(.false.,wfmt(1,1,1,1,ist),wfmt(1,1,1,1,jst),wfir(1,1,ist), &
         wfir(1,1,jst),zrhomt,zrhoir)
        de=evalsv(ist,ik)-evalsv(jst,ik)
        zde=occmax*wkpt(ik)/(de+zi*swidth)
        zvnl=conjg(vnlvv(ist,jst,ik))
        zrvx=zfinp(.false.,zrhomt,zvxmt,zrhoir,zvxir)
        zt1=zvnl-zrvx
        if (spinpol) then
          call oepmag(.false.,wfmt(1,1,1,1,ist),wfmt(1,1,1,1,jst), &
           wfir(1,1,ist),wfir(1,1,jst),zmagmt,zmagir)
! integral of magnetisation dot exchange field
          zmbx=0.d0
          do idm=1,ndmag
            zmbx=zmbx+zfinp(.false.,zmagmt(1,1,1,idm),zbxmt(1,1,1,idm), &
             zmagir(1,idm),zbxir(1,idm))
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
            dvxmt(:,1:nr,ias)=dvxmt(:,1:nr,ias)+dble(zt2*zrhomt(:,1:nr,ias))
            do idm=1,ndmag
              dbxmt(:,1:nr,ias,idm)=dbxmt(:,1:nr,ias,idm) &
               +dble(zt2*zmagmt(:,1:nr,ias,idm))
            end do
          end do
        end do
        dvxir(:)=dvxir(:)+dble(zt2*zrhoir(:))
        do idm=1,ndmag
          dbxir(:,idm)=dbxir(:,idm)+dble(zt2*zmagir(:,idm))
        end do
!$OMP END CRITICAL
! end loop over jst
      end if
    end do
! end loop over ist
  end if
end do
deallocate(apwalm,evecfv,evecsv)
deallocate(wfmt,wfir,wfcr,zrhomt,zrhoir)
if (spinpol) then
  deallocate(zmagmt,zmagir,zvfmt,zfmt)
end if
return
end subroutine

