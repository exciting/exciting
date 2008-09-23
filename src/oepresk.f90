
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
integer nrc,ic,m,idm
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
call getevalsv(vkl(:,ik),evalsv(:,ik))
call getevecfv(vkl(:,ik),vgkl(:,:,:,ik),evecfv)
call getevecsv(vkl(:,ik),evecsv)
! find the matching coefficients
call match(ngk(1,ik),gkc(:,1,ik),tpgkc(:,:,1,ik),sfacgk(:,:,1,ik),apwalm)
! calculate the wavefunctions for all states
call genwfsv(.false.,ngk(1,ik),igkig(:,1,ik),evalsv(:,ik),apwalm,evecfv, &
 evecsv,wfmt,wfir)
!-----------------------------------------------------------!
!     core-conduction overlap density and magnetisation     !
!-----------------------------------------------------------!
do is=1,nspecies
  nrc=nrcmt(is)
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
              call vnlrhomt(.false.,is,wfcr(:,:,1),wfmt(:,:,ias,1,jst), &
               zrhomt(:,:,ias))
              if (spinpol) then
                call vnlrhomt(.false.,is,wfcr(:,:,2),wfmt(:,:,ias,2,jst),zfmt)
                zrhomt(:,1:nrc,ias)=zrhomt(:,1:nrc,ias)+zfmt(:,1:nrc)
              end if
              zvnl=conjg(vnlcv(ic,ias,jst,ik))
              zrvx=zfmtinp(.false.,lmaxvr,nrc,rcmt(:,is),lmmaxvr, &
               zrhomt(:,:,ias),zvxmt(:,:,ias))
              zt1=zvnl-zrvx
! spin-polarised case
              if (spinpol) then
                call oepmagmt(.false.,is,wfcr(:,:,1),wfcr(:,:,2), &
                 wfmt(:,:,ias,1,jst),wfmt(:,:,ias,2,jst),zvfmt)
! integral of magnetisation dot exchange field
                zmbx=0.d0
                do idm=1,ndmag
                  zmbx=zmbx+zfmtinp(.false.,lmaxvr,nrc,rcmt(:,is),lmmaxvr, &
                   zvfmt(:,:,idm),zbxmt(:,:,ias,idm))
                end do
                zt1=zt1-zmbx
! end spin-polarised case
              end if
              zt2=zde*zt1
! residuals for exchange potential and field
!$OMP CRITICAL
              dvxmt(:,1:nrc,ias)=dvxmt(:,1:nrc,ias) &
               +dble(zt2*zrhomt(:,1:nrc,ias))
              do idm=1,ndmag
                dbxmt(:,1:nrc,ias,idm)=dbxmt(:,1:nrc,ias,idm) &
                 +dble(zt2*zvfmt(:,1:nrc,idm))
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
        call vnlrho(.false.,wfmt(:,:,:,:,ist),wfmt(:,:,:,:,jst),wfir(:,:,ist), &
         wfir(:,:,jst),zrhomt,zrhoir)
        de=evalsv(ist,ik)-evalsv(jst,ik)
        zde=occmax*wkpt(ik)/(de+zi*swidth)
        zvnl=conjg(vnlvv(ist,jst,ik))
        zrvx=zfinp(.false.,zrhomt,zvxmt,zrhoir,zvxir)
        zt1=zvnl-zrvx
! spin-polarised case
        if (spinpol) then
          call oepmag(.false.,wfmt(:,:,:,:,ist),wfmt(:,:,:,:,jst), &
           wfir(:,:,ist),wfir(:,:,jst),zmagmt,zmagir)
! integral of magnetisation dot exchange field
          zmbx=0.d0
          do idm=1,ndmag
            zmbx=zmbx+zfinp(.false.,zmagmt(:,:,:,idm),zbxmt(:,:,:,idm), &
             zmagir(:,idm),zbxir(:,idm))
          end do
          zt1=zt1-zmbx
        end if
        zt2=zde*zt1
! residuals for exchange potential and field
!$OMP CRITICAL
        do is=1,nspecies
          nrc=nrcmt(is)
          do ia=1,natoms(is)
            ias=idxas(ia,is)
            dvxmt(:,1:nrc,ias)=dvxmt(:,1:nrc,ias)+dble(zt2*zrhomt(:,1:nrc,ias))
            do idm=1,ndmag
              dbxmt(:,1:nrc,ias,idm)=dbxmt(:,1:nrc,ias,idm) &
               +dble(zt2*zmagmt(:,1:nrc,ias,idm))
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

