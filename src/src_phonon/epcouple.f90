
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine epcouple
use modmain
implicit none
! local variables
integer is,ia,ias,ip
integer js,ja,jas
integer i,j,n,iv(3),isym
integer iq,ik,jk,ikq
integer ist,jst,ir,irc
real(8) vkql(3),x
real(8) t1,t2,t3,t4
complex(8) zt1
! allocatable arrays
real(8), allocatable :: wq(:,:)
real(8), allocatable :: gq(:,:)
complex(8), allocatable :: dynq(:,:,:)
complex(8), allocatable :: ev(:,:)
complex(8), allocatable :: dveffmt(:,:,:,:)
complex(8), allocatable :: dveffir(:,:)
complex(8), allocatable :: zfmt(:,:,:)
complex(8), allocatable :: gzfmt(:,:,:,:)
complex(8), allocatable :: zflm(:)
complex(8), allocatable :: zfir(:)
complex(8), allocatable :: epmat(:,:,:)
complex(8), allocatable :: gmq(:,:,:)
complex(8), allocatable :: b(:,:)
! external functions
real(8) sdelta,gaunt
external sdelta,gaunt
! initialise universal variables
call init0
call init1
call init2
! check k-point grid is commensurate with q-point grid
iv(:)=mod(ngridk(:),ngridq(:))
if ((iv(1).ne.0).or.(iv(2).ne.0).or.(iv(3).ne.0)) then
  write(*,*)
  write(*,'("Error(epcouple): k-point grid incommensurate with q-point grid")')
  write(*,'(" ngridk : ",3I6)') ngridk
  write(*,'(" ngridq : ",3I6)') ngridq
  write(*,*)
  stop
end if
n=3*natmtot
! allocate local arrays
allocate(wq(n,nqpt))
allocate(gq(n,nqpt))
allocate(dynq(n,n,nqpt))
allocate(ev(n,n))
allocate(dveffmt(lmmaxapw,nrcmtmax,natmtot,n))
allocate(dveffir(ngrtot,n))
allocate(zfmt(lmmaxvr,nrcmtmax,natmtot))
allocate(gzfmt(lmmaxvr,nrcmtmax,3,natmtot))
allocate(zflm(lmmaxapw))
allocate(zfir(ngrtot))
allocate(gmq(n,n,nqpt))
allocate(b(n,n))
! read in the density and potentials from file
call readstate
! read Fermi energy from file
call readfermi
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
! get the eigenvalues from file
do ik=1,nkpt
  call getevalsv(vkl(:,ik),evalsv(:,ik))
end do
! compute the occupancies and density of states at the Fermi energy
call occupy
! read in the dynamical matrices
call readdyn(dynq)
! apply the acoustic sum rule
call sumrule(dynq)
! compute the gradients of the effective potential for the rigid-ion term
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
! convert potential to complex spherical harmonic expansion
    irc=0
    do ir=1,nrmt(is),lradstp
      irc=irc+1
      call rtozflm(lmaxvr,veffmt(:,ir,ias),zfmt(:,irc,ias))
    end do
    call gradzfmt(lmaxvr,nrcmt(is),rcmt(:,is),lmmaxvr,nrcmtmax,zfmt(:,:,ias), &
     gzfmt(:,:,:,ias))
  end do
end do
! loop over phonon q-points
do iq=1,nqpt
  write(*,'("Info(epcouple): ",I6," of ",I6," q-points")') iq,nqpt
! diagonalise the dynamical matrix
  call dyndiag(dynq(:,:,iq),wq(:,iq),ev)
! loop over phonon branches
  do j=1,n
! find change effective potential for mode j
    dveffmt(:,:,:,j)=0.d0
    dveffir(:,j)=0.d0
    i=0
    do is=1,nspecies
! prefactor
      t1=2.d0*spmass(is)*wq(j,iq)
      if (t1.gt.1.d-8) then
        t1=1.d0/sqrt(t1)
      else
        t1=0.d0
      end if
      do ia=1,natoms(is)
        ias=idxas(ia,is)
        do ip=1,3
          i=i+1
! read in the Cartesian change in effective potential
          call readdveff(iq,is,ia,ip,zfmt,zfir)
! add the rigid-ion term
          do irc=1,nrcmt(is)
            zfmt(:,irc,ias)=zfmt(:,irc,ias)-gzfmt(:,irc,ip,ias)
          end do
! multiply with eigenvector component and add to total
          zt1=t1*ev(i,j)
          do js=1,nspecies
            do ja=1,natoms(js)
              jas=idxas(ja,js)
              do irc=1,nrcmt(js)
                dveffmt(1:lmmaxvr,irc,jas,j)=dveffmt(1:lmmaxvr,irc,jas,j) &
                 +zt1*zfmt(1:lmmaxvr,irc,jas)
              end do
            end do
          end do
          dveffir(:,j)=dveffir(:,j)+zt1*zfir(:)
        end do
      end do
    end do
! multiply the interstitial potential with the characteristic function
    dveffir(:,j)=dveffir(:,j)*cfunir(:)
! convert muffin-tin potential to spherical coordinates on the lmaxapw covering
    zflm(:)=0.d0
    do is=1,nspecies
      do ia=1,natoms(is)
        ias=idxas(ia,is)
        do irc=1,nrcmt(is)
          zflm(1:lmmaxvr)=dveffmt(1:lmmaxvr,irc,ias,j)
          call zgemv('N',lmmaxapw,lmmaxapw,zone,zbshtapw,lmmaxapw,zflm,1, &
           zzero,dveffmt(:,irc,ias,j),1)
        end do
      end do
    end do
  end do
! zero the phonon linewidths array
  gq(:,iq)=0.d0
! begin parallel loop over non-reduced k-points
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(epmat,jk,vkql,isym) &
!$OMP PRIVATE(ikq,ist,jst,i) &
!$OMP PRIVATE(x,t1,t2,t3,t4)
!$OMP DO
  do ik=1,nkptnr
    allocate(epmat(nstsv,nstsv,n))
! equivalent reduced k-point
    jk=ikmap(ivknr(1,ik),ivknr(2,ik),ivknr(3,ik))
! compute the electron-phonon coupling matrix elements
    call genepmat(iq,vklnr(:,ik),dveffmt,dveffir,epmat)
! k+q-vector in lattice coordinates
    vkql(:)=vklnr(:,ik)+vql(:,iq)
! index to k+q-vector
    call findkpt(vkql,isym,ikq)
    t1=twopi*wkptnr(ik)*(occmax/2.d0)
! loop over second-variational states
    do ist=1,nstsv
      x=(evalsv(ist,ikq)-efermi)/swidth
      t2=sdelta(stype,x)/swidth
! loop over phonon branches
      do i=1,n
        do jst=1,nstsv
          x=(evalsv(jst,jk)-efermi)/swidth
          t3=sdelta(stype,x)/swidth
          t4=dble(epmat(ist,jst,i))**2+aimag(epmat(ist,jst,i))**2
!$OMP CRITICAL
          gq(i,iq)=gq(i,iq)+wq(i,iq)*t1*t2*t3*t4
!$OMP END CRITICAL
        end do
      end do
    end do
    deallocate(epmat)
! end loop over k-points
  end do
!$OMP END DO
!$OMP END PARALLEL
! end loop over phonon q-points
end do
! write the phonon linewidths to file
call writegamma(gq)
! write electron-phonon coupling constants to file
call writelambda(wq,gq)
deallocate(wq,gq,dynq,ev,dveffmt,dveffir)
deallocate(zfmt,gzfmt,zflm,zfir,gmq,b)
return
end subroutine

