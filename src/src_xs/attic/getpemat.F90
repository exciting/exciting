
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_getpemat
  implicit none
contains

  subroutine getpemat(iq,ik,pfilnam,efilnam,nstv,nstc,xou,xuo,pou,puo)
    use modmain
    use modxs
    use modtetra
    use m_getpmat
    use m_getemat
    implicit none
    ! arguments
    integer, intent(in) :: iq,ik,nstv,nstc
    character(*), intent(in) :: pfilnam,efilnam
    complex(8), optional, intent(out) :: xou(:,:,:), xuo(:,:,:)
    complex(8), optional, intent(out) :: pou(:,:,:), puo(:,:,:)
    ! local variables
    character(*), parameter :: thisnam='getpemat'
    complex(8), allocatable :: pm(:,:,:)
    real(8) :: fourpisqt
    integer :: n,igq,j
    logical :: tq0
    logical, external :: tqgamma
    tq0=tqgamma(iq)
    n=ngq(iq)
    fourpisqt=sqrt(fourpi)
    if (tq0) then
       ! Gamma q-point
       nstsv=nstv+nstc
       allocate(pm(3,nstsv,nstsv))
       ! read momentum matrix elements
       call getpmat(ik,vkl0,.true.,trim(pfilnam),pm)
       pou(:,:,:)=pm(:,1:nstv,nstv+1:nstsv)
       puo(:,:,:)=pm(:,nstv+1:nstsv,1:nstv)
       deallocate(pm)
       ! consider symmetric gauge wrt. Coulomb potential
       ! (multiply with v^(1/2))
       ! and normalization wrt. KS eigenvalues (no scissors correction!)
       forall (j=1:3)
          pou(j,:,:)=-pou(j,:,:)/deou(:,:)*fourpisqt
          puo(j,:,:)=-puo(j,:,:)/deuo(:,:)*fourpisqt
       end forall
    end if
    if ((.not.tq0).or.(n.gt.1)) then
       ! read matrix elemets of exponential expression
       call getemat(iq,ik,.true.,trim(efilnam),xou,xuo)
       ! consider symmetric gauge wrt. Coulomb potential (multiply with v^(1/2))
       if (.not.tq0) then
             xou(:,:,1)=xou(:,:,1)/gqc(1,iq)*fourpisqt
             xuo(:,:,1)=xuo(:,:,1)/gqc(1,iq)*fourpisqt
       end if
       if (n.gt.1) then
          forall (igq=2:n)
             xou(:,:,igq)=xou(:,:,igq)/gqc(igq,iq)*fourpisqt
             xuo(:,:,igq)=xuo(:,:,igq)/gqc(igq,iq)*fourpisqt
          end forall
       end if
    end if
  end subroutine getpemat
end module m_getpemat
