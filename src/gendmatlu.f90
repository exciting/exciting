
! Copyright (C) 2007 F. Bultmark, F. Cricchio, L. Nordstrom and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gendmatlu
use modmain
implicit none
! local variables
integer ik,ispn,jspn,i,j,n
integer is,ia,ias,ist,irc
integer l,m1,m2,lm1,lm2
real(8) wo,t1,t2
complex(8) zt1
! automatic arrays
real(8) fr1(nrcmtmax),fr2(nrcmtmax)
real(8) gr(nrcmtmax),cf(3,nrcmtmax)
! allocatable arrays
logical, allocatable :: done(:,:)
complex(8), allocatable :: apwalm(:,:,:,:,:)
complex(8), allocatable :: evecfv(:,:,:)
complex(8), allocatable :: evecsv(:,:)
complex(8), allocatable :: wfmt1(:,:,:,:)
complex(8), allocatable :: wfmt2(:,:,:)
! allocate local arrays
allocate(done(nstfv,nspnfv))
allocate(evecsv(nstsv,nstsv))
allocate(evecfv(nmatmax,nstfv,nspnfv))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
allocate(wfmt1(lmmaxvr,nrcmtmax,nstfv,nspnfv))
allocate(wfmt2(lmmaxvr,nrcmtmax,nspinor))
! zero the LDA+U density matrix
dmatlu(:,:,:,:,:)=0.d0
! begin loop over k-points
do ik=1,nkpt
! get the eigenvectors from file
  call getevecfv(vkl(1,ik),vgkl(1,1,ik,1),evecfv)
  call getevecsv(vkl(1,ik),evecsv)
! find the matching coefficients
  do ispn=1,nspnfv
    call match(ngk(ik,ispn),gkc(1,ik,ispn),tpgkc(1,1,ik,ispn), &
     sfacgk(1,1,ik,ispn),apwalm(1,1,1,1,ispn))
  end do
! begin loop over atoms and species
  do is=1,nspecies
    l=llu(is)
    if (l.lt.0) goto 10
    n=lmmaxvr*nrcmt(is)
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      done(:,:)=.false.
! begin loop over second-variational states
      do j=1,nstsv
        wo=wkpt(ik)*occsv(j,ik)
        if (abs(wo).gt.epsocc) then
          if (tevecsv) then
! generate spinor wavefunction from second-variational eigenvectors
            wfmt2(:,:,:)=0.d0
            i=0
            do ispn=1,nspinor
              if (spinsprl) then
                jspn=ispn
              else
                jspn=1
              end if
              do ist=1,nstfv
                i=i+1
                zt1=evecsv(i,j)
                if (abs(dble(zt1))+abs(aimag(zt1)).gt.epsocc) then
                  if (.not.done(ist,jspn)) then
                    call wavefmt(lradstp,lmaxvr,is,ia,ngk(ik,jspn), &
                     apwalm(1,1,1,1,jspn),evecfv(1,ist,jspn),lmmaxvr, &
                     wfmt1(1,1,ist,jspn))
                    done(ist,jspn)=.true.
                  end if
! add to spinor wavefunction
                  call zaxpy(n,zt1,wfmt1(1,1,ist,jspn),1,wfmt2(1,1,ispn),1)
                end if
              end do
            end do
          else
! spin-unpolarised wavefunction
            call wavefmt(lradstp,lmaxvr,is,ia,ngk(ik,1),apwalm,evecfv(1,j,1), &
             lmmaxvr,wfmt2)
          end if
        end if
        do ispn=1,nspinor
          do jspn=1,nspinor
            do m1=-l,l
              lm1=idxlm(l,m1)
              do m2=-l,l
                lm2=idxlm(l,m2)
                do irc=1,nrcmt(is)
                  zt1=wfmt2(lm1,irc,ispn)*conjg(wfmt2(lm2,irc,jspn))
                  t1=rcmt(irc,is)**2
                  fr1(irc)=dble(zt1)*t1
                  fr2(irc)=aimag(zt1)*t1
                end do
                call fderiv(-1,nrcmt(is),rcmt(1,is),fr1,gr,cf)
                t1=gr(nrcmt(is))
                call fderiv(-1,nrcmt(is),rcmt(1,is),fr2,gr,cf)
                t2=gr(nrcmt(is))
                dmatlu(lm1,lm2,ispn,jspn,ias)=dmatlu(lm1,lm2,ispn,jspn,ias) &
                 +wo*cmplx(t1,t2,8)
              end do
            end do
          end do
        end do
! end loop over second-variational states
      end do
! end loop over atoms
    end do
10 continue
! end loop over species
  end do
! end loop over k-points
end do
! symmetrise the density matrix
call symdmatlu(dmatlu)
deallocate(done,evecfv,evecsv,apwalm,wfmt1,wfmt2)
return
end subroutine

