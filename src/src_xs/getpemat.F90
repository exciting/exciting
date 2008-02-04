
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_getpemat
  implicit none
contains

  subroutine getpemat(iq,ik,pfilnam,efilnam,m12,m34,p12,p34)
    use modmain
    use modxs
    use modtetra
    use m_getpmat
    use m_getemat
    implicit none
    ! arguments
    integer, intent(in) :: iq,ik
    character(*), intent(in) :: pfilnam,efilnam
    complex(8), optional, intent(out) :: m12(:,:,:),p34(:,:,:)
    complex(8), optional, intent(out) :: p12(:,:,:),m34(:,:,:)
    ! local variables
    character(*), parameter :: thisnam='getpemat'
    real(8), parameter :: eps=1.d-8
    complex(8), allocatable :: pm(:,:,:)
    real(8) :: fourpisqt
    integer :: n,igq,j,i1,i2
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
          do i1=1,nst1
             do i2=1,nst2
                if (abs(deou(i1,i2)).ge.epsdfde) then
                   p12(j,i1,i2)=-p12(j,i1,i2)/deou(i1,i2)*fourpisqt
                else
                   p12(j,i1,i2)=zzero
                   if (abs(docc12(i1,i2)).gt.epsocc) then
                      write(*,'("Warning(",a,"): divergent energy denominator: &
                           &q-point, k-point, band indices 1-2:",4i6,g18.10)')&
                           thisnam,iq,ik,i1+istlo1-1,i2+istlo2-1,deou(i1,i2)
                   end if
                end if
                if (abs(deuo(i2,i1)).ge.epsdfde) then
                   p34(j,i2,i1)=-p34(j,i2,i1)/deuo(i2,i1)*fourpisqt
                else
                   p34(j,i2,i1)=zzero
                   if (abs(docc21(i2,i1)).gt.epsocc) then
                      write(*,'("Warning(",a,"): divergent energy denominator: &
                           &q-point, k-point, band indices 3-4:",4i6,g18.10)')&
                           thisnam,iq,ik,i1+istlo1-1,i2+istlo2-1,deuo(i2,i1)
                   end if
                end if
             end do
          end do
       end do
    end if
    if ((.not.tq0).or.(n.gt.1)) then
       ! for BSE(-kernel) matrix elements are calculated on the fly
       if ((task.ge.400).and.(task.le.499)) then
          m12(:,:,:)=xiou(:,:,:)
          m34(:,:,:)=xiuo(:,:,:)
       else
          ! read matrix elemets of exponential expression
          call getemat(iq,ik,.true.,trim(efilnam),m12,m34)
       end if
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
  end subroutine getpemat

end module m_getpemat
