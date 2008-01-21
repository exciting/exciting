
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_getpemat2
  implicit none
contains

  subroutine getpemat2(iq,ik,pfilnam,efilnam,m12,m34,p12,p34)
    use modmain
    use modxs
    use modtetra
    use m_getpmat
    use m_getemat2
    implicit none
    ! arguments
    integer, intent(in) :: iq,ik
    character(*), intent(in) :: pfilnam,efilnam
    complex(8), optional, intent(out) :: m12(:,:,:),p34(:,:,:)
    complex(8), optional, intent(out) :: p12(:,:,:),m34(:,:,:)
    ! local variables
    character(*), parameter :: thisnam = 'getpemat2'
    complex(8), allocatable :: pm(:,:,:)
    real(8) :: fourpisqt
    integer :: n,igq,j
    logical :: tq0
    logical, external :: tqgamma
    tq0=tqgamma(iq)
    n=ngq(iq)
    fourpisqt=sqrt(fourpi)
    if (tq0.and.(.not.present(p12))) then
       write(*,*)
       write(*,'("Error(",a,"): Gamma q-point but momentum matrix elements not &
            &requested.")') thisnam
       write(*,*)
       call terminate
    end if
    if (tq0) then
       ! Gamma q-point
       allocate(pm(3,nstsv,nstsv))
       ! read momentum matrix elements
       call getpmat(ik,vkl0,.true.,trim(pfilnam),pm)
       p12(:,:,:)=pm(:,istlo1:isthi1,istlo2:isthi2)
       p34(:,:,:)=pm(:,istlo3:isthi3,istlo4:isthi4)
       deallocate(pm)
       ! consider symmetric gauge wrt. Coulomb potential
       ! (multiply with v^(1/2))
       ! and normalization wrt. KS eigenvalues (no scissors correction!)
       do j=1,3
!!$          where(abs(docc12).gt.epsocc)
!!$             p12(j,:,:)=-p12(j,:,:)/deou(:,:)*fourpisqt
!!$          elsewhere
!!$             p12(j,:,:)=zzero
!!$          end where
!!$          where(abs(docc21).gt.epsocc)
!!$             p34(j,:,:)=-p34(j,:,:)/deuo(:,:)*fourpisqt
!!$          elsewhere
!!$             p34(j,:,:)=zzero
!!$          end where
          where(abs(deou).gt.epsocc)
             p12(j,:,:)=-p12(j,:,:)/deou(:,:)*fourpisqt
          elsewhere
             p12(j,:,:)=zzero
          end where
          where(abs(deuo).gt.epsocc)
             p34(j,:,:)=-p34(j,:,:)/deuo(:,:)*fourpisqt
          elsewhere
             p34(j,:,:)=zzero
          end where
       end do
    end if
    if ((.not.tq0).or.(n.gt.1)) then
       ! read matrix elemets of exponential expression
       call getemat2(iq,ik,.true.,trim(efilnam),m12,m34)
       ! consider symmetric gauge wrt. Coulomb potential (multiply with v^(1/2))
       if (.not.tq0) then
             m12(:,:,1)=m12(:,:,1)/gqc(1,iq)*fourpisqt
             m34(:,:,1)=m34(:,:,1)/gqc(1,iq)*fourpisqt
       end if
       if (n.gt.1) then
          forall (igq=2:n)
             m12(:,:,igq)=m12(:,:,igq)/gqc(igq,iq)*fourpisqt
             m34(:,:,igq)=m34(:,:,igq)/gqc(igq,iq)*fourpisqt
          end forall
       end if
    end if
  end subroutine getpemat2

end module m_getpemat2
