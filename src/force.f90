
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: force
! !INTERFACE:
subroutine force
! !USES:
use modmain
! !DESCRIPTION:
!   Computes the various contributions to the atomic forces. In principle, the
!   force acting on a nucleus is simply the gradient at that site of the
!   classical electrostatic potential from the other nuclei and the electronic
!   density. This is a result of the Hellmann-Feynman theorem. However because
!   the basis set is dependent on the nuclear coordinates and is not complete,
!   the Hellman-Feynman force is inacurate and corrections to it are required.
!   The first is the core correction which arises because the core wavefunctions
!   were determined by neglecting the non-spherical parts of the effective
!   potential $v_s$. Explicitly this is given by
!   $$ {\bf F}_{\rm core}^{\alpha}=\int_{\rm MT_{\alpha}} v_{\rm s}({\bf r})
!    \nabla\rho_{\rm core}^{\alpha}({\bf r})\,d{\bf r} $$
!   for atom $\alpha$. The second is due to the position dependence of the APW
!   functions, and is derived by considering the change in total energy if the
!   eigenvector coefficients were fixed and the APW functions themselves were
!   changed. This so-called incomplete basis set (IBS) correction to the $i$th
!   first-variational eigenvalue, $\epsilon^i$, is given by
!   \begin{align*}
!    {\bf F}^{i\alpha}_{\rm FV}=
!    &-\sum_{\bf G,G'}\left\{\frac{1}{2}({\bf G+k})\cdot({\bf G'+k})
!    -\epsilon_i\right\}\Phi^{i*}_{\bf G+k}\Phi^{i}_{\bf G'+k}i({\bf G-G'})
!    \tilde{\Theta}_{\alpha}({\bf G-G'})e^{-i({\bf G-G'})\cdot
!    {\bf r}_{\alpha}}\\
!    &-\sum_{\bf G}{\sum_{\bf G'}}'\left\{\left(
!    H^{{\rm MT}_{\alpha}}_{\bf G+k,G'+k}-\epsilon_i
!    O^{{\rm MT}_{\alpha}}_{\bf G+k,G'+k}\right)\Phi^{i*}_{\bf G+k}
!    \Phi^i_{\bf G'+k}i({\bf G+k})+{\rm c.c.}\right\},
!   \end{align*}
!   where $\Phi^i$ is the first-variational eigenvector;
!   $\tilde{\Theta}_{\alpha}$ is the form factor of the smooth step function
!   for muffin-tin $\alpha$; $H^{{\rm MT}_{\alpha}}$ and $O^{{\rm MT}_{\alpha}}$
!   are the muffin-tin Hamiltonian and overlap matrices, respectively; and the
!   primed sum is over both the ${\bf G}$-vectors and local-orbital indices.
!   Having obtained the derivative of the first-variational eigenvalues with
!   respect to the basis set, the derivatives of the second-variational
!   eigenvalues can be obtained with
!   $$ {\bf F}^{i\alpha}_{\rm SV}=\sum_j \left|\psi^i_j\right|^2
!    {\bf F}^{j\alpha}_{\rm FV}, $$
!   where $\psi^i$ is the second-variational eigenvector. Finally the IBS
!   correction is
!   $$ {\bf F}^{\alpha}_{\rm IBS}=\sum_i n_i {\bf F}^{i\alpha}_{\rm SV}
!    +\int_{\rm MT_{\alpha}}v_{\rm s}({\bf r})\nabla\left[\rho({\bf r})
!    -\rho^{\alpha}_{\rm core}({\bf r})\right]\,d{\bf r}, $$
!   where $n_i$ are the state occupancies. See routines {\tt hmlaa},
!   {\tt olpaa}, {\tt hmlalo}, {\tt olpalo}, {\tt energy}, {\tt seceqn} and
!   {\tt gencfun}.
!
! !REVISION HISTORY:
!   Created January 2004 (JKD)
!EOP
!BOC
implicit none
! local variables
integer ik,is,ia,ias,nr,i,ig
real(8) sum,t1,t2,t3,t4
real(8) cpu0,cpu1
! allocatable arrays
real(8), allocatable :: rfmt(:,:)
real(8), allocatable :: grfmt(:,:,:)
real(8), allocatable :: ff(:,:)
! external functions
real(8) rfmtinp
external rfmtinp
call cpu_time(cpu0)
allocate(rfmt(lmmaxvr,nrmtmax))
allocate(grfmt(lmmaxvr,nrmtmax,3))
!--------------------------------!
!     Hellmann-Feynman force     !
!--------------------------------!
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
! compute the gradient of the Coulomb potential
    call gradrfmt(1,nrmt(is),spr(1,is),lmmaxvr,nrmtmax,vclmt(1,1,ias),grfmt)
    forcehf(:,ias)=-spzn(is)*grfmt(1,1,:)*y00
  end do
end do
! symmetrise Hellmann-Feynman force
call symvect(forcehf)
!--------------------------------------!
!     core correction to the force     !
!--------------------------------------!
rfmt(:,:)=0.d0
do is=1,nspecies
  nr=nrmt(is)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
! compute the gradient of the core density
    rfmt(1,1:nr)=rhocr(1:nr,ias)/y00
    call gradrfmt(1,nr,spr(1,is),lmmaxvr,nrmtmax,rfmt,grfmt)
    do i=1,3
      forcecr(i,ias)=rfmtinp(1,1,nr,spr(1,is),lmmaxvr,veffmt(1,1,ias), &
       grfmt(1,1,i))
    end do
  end do
end do
! symmetrise core correction force
call symvect(forcecr)
!-------------------------------------!
!     IBS correction to the force     !
!-------------------------------------!
! set the IBS forces to zero
forceibs(:,:)=0.d0
if (tfibs) then
  allocate(ff(ngvec,nspecies))
! integral of effective potential with gradient of valence density
  do is=1,nspecies
    nr=nrmt(is)
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      rfmt(:,1:nr)=rhomt(:,1:nr,ias)
      rfmt(1,1:nr)=rfmt(1,1:nr)-rhocr(1:nr,ias)/y00
      call gradrfmt(lmaxvr,nr,spr(1,is),lmmaxvr,nrmtmax,rfmt,grfmt)
      do i=1,3
        t1=rfmtinp(1,lmaxvr,nr,spr(1,is),lmmaxvr,veffmt(1,1,ias),grfmt(1,1,i))
        forceibs(i,ias)=forceibs(i,ias)+t1
      end do
    end do
  end do
! step function form factors
  t1=fourpi/omega
  t2=cfdamp/gmaxvr
  do is=1,nspecies
    do ig=1,ngvec
      if (gc(ig).gt.epslat) then
        if (cfdamp.ne.0.d0) then
! use damping if required
          t3=exp(-(t2*gc(ig))**2)
        else
          t3=1.d0
        end if
        t4=gc(ig)*rmt(is)
        ff(ig,is)=t1*t3*(sin(t4)-t4*cos(t4))/(gc(ig)**3)
      else
        ff(ig,is)=t1*(rmt(is)**3)/3.d0
      end if
    end do
  end do
! compute k-point dependent contribution to the IBS force
!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO
  do ik=1,nkpt
    call forcek(ik,ff)
  end do
!$OMP END DO
!$OMP END PARALLEL
! symmetrise IBS force
  call symvect(forceibs)
  deallocate(ff)
end if
! total force
do ias=1,natmtot
  forcetot(:,ias)=forcehf(:,ias)+forcecr(:,ias)+forceibs(:,ias)
end do
! symmetrise total force
call symvect(forcetot)
! remove net total force (center of mass should not move)
do i=1,3
  sum=0.d0
  do ias=1,natmtot
    sum=sum+forcetot(i,ias)
  end do
  sum=sum/dble(natmtot)
  forcetot(i,:)=forcetot(i,:)-sum
end do
! compute maximum force magnitude over all atoms
forcemax=0.d0
do ias=1,natmtot
  t1=sqrt(forcetot(1,ias)**2+forcetot(2,ias)**2+forcetot(3,ias)**2)
  if (t1.gt.forcemax) forcemax=t1
end do
deallocate(rfmt,grfmt)
call cpu_time(cpu1)
timefor=timefor+cpu1-cpu0
return
end subroutine
!EOC

