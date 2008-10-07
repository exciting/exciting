
! Copyright (C) 2007 F. Bultmark, F. Cricchio, L. Nordstrom and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gendmat(tspndg,tlmdg,lmin,lmax,is,ia,ngp,apwalm,evecfv,evecsv,ld, &
 dmat)
use modmain
implicit none
! arguments
logical, intent(in) :: tspndg
logical, intent(in) :: tlmdg
integer, intent(in) :: lmin
integer, intent(in) :: lmax
integer, intent(in) :: is
integer, intent(in) :: ia
integer, intent(in) :: ngp(nspnfv)
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv)
complex(8), intent(in) :: evecfv(nmatmax,nstfv,nspnfv)
complex(8), intent(in) :: evecsv(nstsv,nstsv)
integer, intent(in) :: ld
complex(8), intent(out) :: dmat(ld,ld,nspinor,nspinor,nstsv)
! local variables
integer ispn,jspn,lmmax
integer l,m1,m2,lm1,lm2
integer i,j,n,ist,irc
real(8) t1,t2
complex(8) zt1
! automatic arrays
real(8) fr1(nrcmtmax),fr2(nrcmtmax)
real(8) gr(nrcmtmax),cf(3,nrcmtmax)
! allocatable arrays
logical, allocatable :: done(:,:)
complex(8), allocatable :: wfmt1(:,:,:,:)
complex(8), allocatable :: wfmt2(:,:,:)
if (lmin.lt.0) then
  write(*,*)
  write(*,'("Error(gendmat): lmin < 0 : ",I8)') lmin
  write(*,*)
  stop
end if
if (lmax.gt.lmaxapw) then
  write(*,*)
  write(*,'("Error(gendmat): lmax > lmaxapw : ",2I8)') lmax,lmaxapw
  write(*,*)
  stop
end if
lmmax=(lmax+1)**2
! allocate local arrays
allocate(done(nstfv,nspnfv))
allocate(wfmt1(lmmax,nrcmtmax,nstfv,nspnfv))
allocate(wfmt2(lmmax,nrcmtmax,nspinor))
! zero the density matrix
dmat(:,:,:,:,:)=0.d0
n=lmmax*nrcmt(is)
done(:,:)=.false.
! begin loop over second-variational states
do j=1,nstsv
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
            call wavefmt(lradstp,lmax,is,ia,ngp(jspn),apwalm(:,:,:,:,jspn), &
             evecfv(:,ist,jspn),lmmax,wfmt1(:,:,ist,jspn))
            done(ist,jspn)=.true.
          end if
! add to spinor wavefunction
          call zaxpy(n,zt1,wfmt1(:,:,ist,jspn),1,wfmt2(:,:,ispn),1)
        end if
      end do
    end do
  else
! spin-unpolarised wavefunction
    call wavefmt(lradstp,lmax,is,ia,ngp,apwalm,evecfv(:,j,1),lmmax,wfmt2)
  end if
  do ispn=1,nspinor
    do jspn=1,nspinor
      if (tspndg.and.(ispn.ne.jspn)) goto 20
      do l=lmin,lmax
        do m1=-l,l
          lm1=idxlm(l,m1)
          do m2=-l,l
            lm2=idxlm(l,m2)
            if (tlmdg.and.(lm1.ne.lm2)) goto 10
            do irc=1,nrcmt(is)
              zt1=wfmt2(lm1,irc,ispn)*conjg(wfmt2(lm2,irc,jspn))
              t1=rcmt(irc,is)**2
              fr1(irc)=dble(zt1)*t1
              fr2(irc)=aimag(zt1)*t1
            end do
            call fderiv(-1,nrcmt(is),rcmt(:,is),fr1,gr,cf)
            t1=gr(nrcmt(is))
            call fderiv(-1,nrcmt(is),rcmt(:,is),fr2,gr,cf)
            t2=gr(nrcmt(is))
            dmat(lm1,lm2,ispn,jspn,j)=cmplx(t1,t2,8)
10 continue
          end do
        end do
      end do
20 continue
    end do
  end do
! end loop over second-variational states
end do
deallocate(done,wfmt1,wfmt2)
return
end subroutine

