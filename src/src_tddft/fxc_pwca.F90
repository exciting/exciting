
! Copyright (C) 2002-2007 S. Sagmeister, J. K. Dewhurst, S. Sharma and 
! C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

module m_fxc_alda
  implicit none
contains

  subroutine fxc_alda(iq,msiz,fxcg)
    use modmain
    use modtddft
    use m_ftfun
    implicit none
    ! arguments
    integer, intent(in) :: iq,msiz
    complex(8), intent(out) :: fxcg(:,:)
    ! local variables
    character(*), parameter :: thisnam = 'fxc_alda'
    complex(8), allocatable :: fxcft1(:)
    integer :: sh(2),ig,igq1,igq2,iv1(3),iv(3)

    sh=shape(fxcg)
    if ((sh(1).lt.msiz).or.(sh(2).lt.msiz)) then
       write(unitout,'(a,2i9,a,i9,a)') 'Error('//trim(thisnam)//'): size of &
            &fxc is to small (required)', sh, '(', msiz, ')'
       stop
    end if

    if (allocated(fxcmt)) deallocate(fxcmt)
    if (allocated(fxcir)) deallocate(fxcir)
    allocate(fxcmt(lmmaxvr,nrmtmax,natmtot))
    allocate(fxcir(ngrtot))

    ! calculate exchange-correlation kernel in real space
    call kernxc

    ! Fourier transform of muffin-tin and interstitial kernel
    allocate(fxcft1(ngvec))
    call ftfun(ngvec,.true.,.true.,fxcir,fxcmt,fxcft1)

    ! transform G''=G-G' to G and G'
    do igq1=1,msiz
       iv1(:)=ivg(:,igqig(igq1,iq))
       do igq2=1,msiz
          iv(:)=iv1(:)- ivg(:,igqig(igq2,iq))
          ig = ivgig(iv(1),iv(2),iv(3))
          fxcg(igq1,igq2)=fxcft1(ig)
          ! renormalization to symmetrized quantity wrt. G-space
          fxcg(igq1,igq2)=fxcg(igq1,igq2)*(gqc(igq1,iq)*gqc(igq2,iq))/fourpi
       end do
    end do

    ! deallocate
    deallocate(fxcmt,fxcir,fxcft1)

  end subroutine fxc_alda

end module m_fxc_alda

!BOP
! !ROUTINE: xcd_pwca
! !INTERFACE:
subroutine xcd_pwca(n,rho,dvx,dvc)
! !INPUT/OUTPUT PARAMETERS:
!   n     : number of density points (in,integer)
!   rho   : charge density (in,real(n))
!   dvx   : exchange potential derivative (out,real(n))
!   dvc   : correlation potential derivative (out,real(n))
! !DESCRIPTION:
!   Spin-unpolarised exchange-correlation potential derivative of the 
!   Perdew-Wang
!   parameterisation of the Ceperley-Alder electron gas,
!   Phys. Rev. B 45, 13244
!   (1992) and Phys. Rev. Lett. 45, 566 (1980). Based upon the routine
!   {\tt xc_pwca}.
!
! !REVISION HISTORY:
!   Created February 2007 (Sagmeister)
!EOP
!BOC
  implicit none
  ! arguments
  integer, intent(in) :: n
  real(8), intent(in) :: rho(n)
  real(8), intent(out) :: dvx(n)
  real(8), intent(out) :: dvc(n)
  ! local variables
  integer i
  real(8), parameter :: pi=3.1415926535897932385d0
  real(8), parameter :: thrd=1.d0/3.d0
  real(8), parameter :: thrd2=2.d0/3.d0
  ! beyond RPA
  real(8), parameter :: p=1.d0
  real(8) a(3),a1(3),b1(3),b2(3),b3(3),b4(3)
  data a  / 0.0310907d0, 0.01554535d0, 0.0168869d0 /
  data a1 / 0.21370d0,   0.20548d0,    0.11125d0   /
  data b1 / 7.5957d0,   14.1189d0,    10.357d0     /
  data b2 / 3.5876d0,    6.1977d0,     3.6231d0    /
  data b3 / 1.6382d0,    3.3662d0,     0.88026d0   /
  data b4 / 0.49294d0,   0.62517d0,    0.49671d0   /
  real(8) r,rs,srs
  real(8) drec0,ddrec0
  real(8) q0(3),q1(3),q1p,q1pp,q11t
  if (n.le.0) then
     write(*,*)
     write(*,'("Error(fxc_pwca): invalid n : ",I8)') n
     write(*,*)
     stop
  end if
  do i=1,n
     if (rho(i).gt.1.d-12) then ! *** check if derivative vanishes in this case
        r=rho(i)
        rs=(3.d0/(4.d0*pi*r))**thrd
        srs=sqrt(rs)
        ! exchange potential derivative
        dvx(i)=-thrd*(3/pi)**thrd*r**(-thrd2)
        ! correlation potential derivative
        q0(1)=-2.d0*a(1)*(1.d0+a1(1)*rs)
        q1(1)=2.d0*a(1)*(b1(1)*srs+b2(1)*rs+b3(1)*(srs**3)+b4(1)*rs**(p+1.d0))
        q1p=a(1)*(b1(1)/srs+2.d0*b2(1)+3.d0*b3(1)*srs+2.d0*(p+1.d0)*b4(1)* &
             rs**p)
        q1pp=a(1)*(-b1(1)/(srs**3)+3.d0*b3(1)/srs+4.d0*p*(p+1)*b4(1)* &
             rs**(p-1))/2.d0
        drec0=-2.d0*a(1)*a1(1)*log(1.d0+1.d0/q1(1))-q0(1)*q1p/(q1(1)**2+q1(1))
        q11t=1.d0/(q1(1)*(1.d0+q1(1)))
        ddrec0=-q11t*(-4.d0*a(1)*a1(1)*q1p+q0(1)*(1+2.d0*q1(1))*q11t*(q1p)**2 &
             +q0(1)*q1pp)
        dvc(i)=(1/9.d0)*(rs/r)*(-2.d0*drec0+rs*ddrec0)
     else
        dvx(i)=0.d0
        dvc(i)=0.d0
     end if
  end do
end subroutine xcd_pwca
!EOC
