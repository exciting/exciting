
! Copyright (C) 2007 F. Bultmark, F. Cricchio and L. Nordstrom.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genvmatlu
use modmain
! local variables
integer ias
real(8) t1,t2
! allocatable arrays
real(8), allocatable :: enfll(:)
complex(8), allocatable :: vmfll(:,:,:,:,:)
! fully localised limit (FLL) or around mean field (AFM)
if ((ldapu.eq.1).or.(ldapu.eq.2)) then
  call genvmatlu_12
  return
end if
! interpolation between the two (PRB 67, 153106 (2003))
if (ldapu.eq.3) then
  allocate(enfll(natmtot))
  allocate(vmfll(lmmaxlu,lmmaxlu,nspinor,nspinor,natmtot))
  ldapu=1
  call genvmatlu_12
  vmfll(:,:,:,:,:)=vmatlu(:,:,:,:,:)
  enfll(:)=engyalu(:)
  ldapu=2
  call genvmatlu_12
! reset ldapu value
  ldapu=3
  do ias=1,natmtot
    t1=alphalu(ias)
    t2=1.d0-t1
    vmatlu(:,:,:,:,ias)=t1*vmfll(:,:,:,:,ias)+t2*vmatlu(:,:,:,:,ias)
    engyalu(ias)=t1*enfll(ias)+t2*engyalu(ias)
  end do
  deallocate(enfll,vmfll)
  return
end if
write(*,*)
write(*,'("Error(genvmatlu): invalid ldapu : ",I8)') ldapu
write(*,*)
stop
end subroutine

subroutine genvmatlu_12
use modmain
implicit none
! local variables
integer is,ia,ias,ispn,jspn
integer l,m1,m2,m3,m4,nm
integer lm1,lm2,lm3,lm4
real(8) u,j,n,n0
real(8) mg(3),mg0(3),mg2
real(8) edc,sum1,sum2
complex(8) zt1,zt2
! automatic arrays
real(8) vee(-lmaxlu:lmaxlu,-lmaxlu:lmaxlu,-lmaxlu:lmaxlu,-lmaxlu:lmaxlu)
complex(8) dm(lmmaxlu,lmmaxlu,nspinor,nspinor)
complex(8) dmt(nspinor,nspinor)
! zero the LDA+U potential for each atom
vmatlu(:,:,:,:,:)=0.d0
! zero the LDA+U energy for each atom
engyalu(:)=0.d0
! zero the interpolation constants
alphalu(:)=0.d0
! begin loop over species
do is=1,nspecies
! define LDA+U parameters
  l=llu(is)
  if (l.lt.0) goto 10
  nm=2*l+1
  u=ujlu(1,is)
  j=ujlu(2,is)
  if ((abs(u).lt.1.d-10).and.(abs(j).lt.1.d-10)) goto 10
! calculate the Coulomb matrix elements
  call genveelu(l,u,j,lmaxlu,vee)
! begin loop over atoms
  do ia=1,natoms(is)
    ias=idxas(ia,is)
! copy the density matrix for this atom
    dm(:,:,:,:)=dmatlu(:,:,:,:,ias)
! trace of density matrix for each spin
    dmt(:,:)=0.d0
    do ispn=1,nspinor
      do jspn=1,nspinor
        do m1=-l,l
          lm1=idxlm(l,m1)
          dmt(ispn,jspn)=dmt(ispn,jspn)+dm(lm1,lm1,ispn,jspn)
        end do
      end do
    end do
! trace over spin
    n=dble(dmt(1,1))
    if (spinpol) n=n+dble(dmt(2,2))
    n0=n/dble(nspinor*nm)
! magnetisation
    if (spinpol) then
      mg(:)=0.d0
      mg(3)=dble(dmt(1,1)-dmt(2,2))
! non-collinear terms
      if (ncmag) then
        mg(1)=dble(dmt(1,2)+dmt(2,1))
        mg(2)=dble(zi*(dmt(1,2)-dmt(2,1)))
      end if
      mg0(:)=mg(:)/dble(nspinor*nm)
    end if
! around mean field (AFM) approach
    if (ldapu.eq.2) then
! modify density matrices
      do m1=-l,l
        lm1=idxlm(l,m1)
        if (spinpol) then
          dm(lm1,lm1,1,1)=dm(lm1,lm1,1,1)-(n0+mg0(3))
          dm(lm1,lm1,2,2)=dm(lm1,lm1,2,2)-(n0-mg0(3))
! non-collinear terms
          if (ncmag) then
            dm(lm1,lm1,1,2)=dm(lm1,lm1,1,2)-(mg0(1)-zi*mg0(2))
            dm(lm1,lm1,2,1)=dm(lm1,lm1,2,1)-(mg0(1)+zi*mg0(2))
          end if
        else
! spin-unpolarised case
          dm(lm1,lm1,1,1)=dm(lm1,lm1,1,1)-n0
        end if
      end do
! determine alpha (PRB 67,153106 (2003))
      sum1=0.d0
      do ispn=1,nspinor
        do m1=-l,l
          lm1=idxlm(l,m1)
          do jspn=1,nspinor
            do m2=-l,l
              lm2=idxlm(l,m2)
              sum1=sum1+dble(dm(lm1,lm2,ispn,jspn)*dm(lm2,lm1,jspn,ispn))
            end do
          end do
        end do
      end do
      if (spinpol) then
        mg2=mg(3)**2
        if (ncmag) mg2=mg2+mg(1)**2+mg(2)**2
      else
        mg2=0.d0
      end if
      sum2=n*(1.d0-0.5d0*n/dble(nm))-0.5d0*mg2/dble(nm)
      if (abs(sum2).gt.1.d-14) then
        alphalu(ias)=sum1/sum2
      else
        alphalu(ias)=0.d0
      end if
    end if
! calculation of LDA+U potential and energy
! begin loops over m1 and m2
    do m1=-l,l
      lm1=idxlm(l,m1)
      do m2=-l,l
        lm2=idxlm(l,m2)
! begin loops over m3 and m4
        do m3=-l,l
          lm3=idxlm(l,m3)
          do m4=-l,l
            lm4=idxlm(l,m4)
            do ispn=1,nspinor
              do jspn=1,nspinor
                zt1=dm(lm2,lm1,ispn,ispn)*dm(lm4,lm3,jspn,jspn)
                zt2=dm(lm4,lm1,jspn,ispn)*dm(lm2,lm3,ispn,jspn)
                engyalu(ias)=engyalu(ias)+dble(zt1-zt2)*vee(m1,m3,m2,m4)
                vmatlu(lm1,lm2,ispn,ispn,ias)=vmatlu(lm1,lm2,ispn,ispn,ias) &
                 +dm(lm4,lm3,jspn,jspn)*vee(m1,m3,m2,m4)
                vmatlu(lm1,lm4,ispn,jspn,ias)=vmatlu(lm1,lm4,ispn,jspn,ias) &
                 -dm(lm2,lm3,ispn,jspn)*vee(m1,m3,m2,m4)
              end do
            end do
! end loops over m3 and m4
          end do
        end do
! end loops over m1 and m2
      end do
    end do
! multiply energy by factor 1/2
    engyalu(ias)=0.5d0*engyalu(ias)
! fully localised limit (FLL) approach: double counting corrections
    if (ldapu.eq.1) then
      if (spinpol) then
! spin-polarised
        if (ncmag) then
! non-collinear case
! correction to the energy 
          edc=0.5d0*u*n*(n-1.d0)
          edc=edc-0.5d0*j*dble(dmt(1,1)*(dmt(1,1)-1.d0))
          edc=edc-0.5d0*j*dble(dmt(2,2)*(dmt(2,2)-1.d0))
          edc=edc-0.5d0*j*dble(dmt(1,2)*dmt(2,1))
          edc=edc-0.5d0*j*dble(dmt(2,1)*dmt(1,2))
! correction to the potential
          do m1=-l,l
            lm1=idxlm(l,m1)
            vmatlu(lm1,lm1,1,1,ias)=vmatlu(lm1,lm1,1,1,ias) &
             -u*(n-0.5d0)+j*(dmt(1,1)-0.5d0)
            vmatlu(lm1,lm1,2,2,ias)=vmatlu(lm1,lm1,2,2,ias) &
             -u*(n-0.5d0)+j*(dmt(2,2)-0.5d0)
            vmatlu(lm1,lm1,1,2,ias)=vmatlu(lm1,lm1,1,2,ias)+j*dmt(1,2)
            vmatlu(lm1,lm1,2,1,ias)=vmatlu(lm1,lm1,2,1,ias)+j*dmt(2,1)
          end do
        else
! collinear case
! correction to the energy 
          edc=0.5d0*u*n*(n-1.d0)
          edc=edc-0.5d0*j*dble(dmt(1,1)*(dmt(1,1)-1.d0))
          edc=edc-0.5d0*j*dble(dmt(2,2)*(dmt(2,2)-1.d0))
! correction to the potential
          do m1=-l,l
            lm1=idxlm(l,m1)
            vmatlu(lm1,lm1,1,1,ias)=vmatlu(lm1,lm1,1,1,ias) &
             -u*(n-0.5d0)+j*(dmt(1,1)-0.5d0)
            vmatlu(lm1,lm1,2,2,ias)=vmatlu(lm1,lm1,2,2,ias) &
             -u*(n-0.5d0)+j*(dmt(2,2)-0.5d0)
          end do
        end if
      else
!spin-unpolarised
! correction to the energy
        edc=0.5d0*u*n*(n-1.d0)
        edc=edc-0.5d0*j*n*(n-1.d0)
! correction to the potential
        do m1=-l,l
          lm1=idxlm(l,m1)
          vmatlu(lm1,lm1,1,1,ias)=vmatlu(lm1,lm1,1,1,ias)-u*(n-0.5d0) &
           +j*(n-0.5d0)
        end do
      end if
      engyalu(ias)=engyalu(ias)-edc
    end if
! trace of dmatlu times vmatlu
    sum1=0.d0
    do ispn=1,nspinor
      do m1=-l,l
        lm1=idxlm(l,m1)
        do jspn=1,nspinor
          do m2=-l,l
            lm2=idxlm(l,m2)
            sum1=sum1+dble(dm(lm1,lm2,ispn,jspn)*vmatlu(lm2,lm1,jspn,ispn,ias))
          end do
        end do
      end do
    end do
! subtract contribution to the energy of LDA+U potential
    engyalu(ias)=engyalu(ias)-sum1
! end loop over atoms
  end do
10 continue
! end loop over species
end do
! symmetrise the potential
call symdmat(lmaxlu,lmmaxlu,vmatlu)
return
end subroutine

