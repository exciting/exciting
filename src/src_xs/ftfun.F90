
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_ftfun
  implicit none
contains

  subroutine ftfun(ng,tir,tmt,gir,gmt,ftg)
    use modmain
    implicit none
    ! arguments
    integer, intent(in) :: ng
    logical, intent(in) :: tir
    logical, intent(in) :: tmt
    complex(8), intent(in) :: gir(:)
    complex(8), intent(in) :: gmt(:,:,:)
    complex(8), intent(out) :: ftg(:)
    ! local variales
    complex(8), allocatable :: zfft(:)
    complex(8) :: zt1,zt2
    real(8), allocatable :: jbesslh(:),jbessl(:,:)
    real(8), allocatable :: r1(:),r2(:),fr(:),fr2(:),gr(:),cf(:,:)
    real(8) :: t1
    integer :: ig,ifg,is,ia,ias,ir,nr,l,m,lm

    ftg(:)=zzero

    ! interstitial part
    if (tir) then
       allocate(zfft(ngrtot))
       ! multiply effective potential with smooth characteristic function
       zfft(:)=gir(:)*cfunir(:)
       ! Fourier transform to G-space
       call zfftifc(3,ngrid,-1,zfft)
       do ig=1,ng
          ifg=igfft(ig)
          ftg(ig)=ftg(ig)+zfft(ifg)
       end do
       deallocate(zfft)
    end if

    ! muffin-tin part
    if (tmt) then
       ! allocate array for Bessel functions
       allocate(jbessl(0:lmaxvr,nrmtmax))
       allocate(jbesslh(0:lmaxvr))
       allocate(r1(nrmtmax),r2(nrmtmax),fr(nrmtmax),fr2(nrmtmax),gr(nrmtmax), &
            cf(3,nrmtmax))
       ! loop over G vectors
       do ig=1,ng
          ! loop over species
          do is=1,nspecies
             nr=nrmt(is)
             do ir=1,nr
                r1(ir)=spr(ir,is)
                r2(ir)=r1(ir)**2
             end do
             ! calculate bessel functions j_l(|G||r|)
             do ir=1,nr
                call sbessel(lmaxvr,gc(ig)*r1(ir),jbesslh)
                jbessl(:,ir)=jbesslh(:)
             end do
             ! loop over atoms
             do ia=1,natoms(is)
                ias=idxas(ia,is)
                zt1=zzero
                zt2=zzero
                do l=0,lmaxvr
                   do m=-l,l
                      lm=idxlm(l,m)
                      do ir=1,nr
                         ! integrand
                         fr(ir)= r2(ir)*jbessl(l,ir)*dble(gmt(lm,ir,ias))
                         fr2(ir)= r2(ir)*jbessl(l,ir)*aimag(gmt(lm,ir,ias))
                      end do
                      if (any(fr.ne.0.d0)) then
                         ! integration
                         call fderiv(-1,nr,spr(1,is),fr,gr,cf)
                         zt1=zt1+gr(nr)*conjg(zil(l))*ylmg(lm,ig)
                      end if
                      if (any(fr2.ne.0.d0)) then
                         ! integration
                         call fderiv(-1,nr,spr(1,is),fr2,gr,cf)
                         zt2=zt2+gr(nr)*conjg(zil(l))*ylmg(lm,ig)
                      end if
                   end do ! m
                end do ! l
                ! form factor summation
                ftg(ig)=ftg(ig)+(fourpi/omega)*conjg(sfacg(ig,ias))*(zt1+zi*zt2)
             end do ! ia
          end do ! is
       end do ! ig
       deallocate(jbessl,jbesslh,r1,r2,fr,fr2,gr,cf)
    end if

  end subroutine ftfun

end module m_ftfun


!*** make routine more general ***

! pass the following arguments to routine:

! * p-vector
! * |G+p| values
! * structure factors for -(G+q) !!!!! (G+p) is calculated in routine
! * spherical harmonics Y_lm(^G+q^)




! INTERFACE: genjlgpr(lmax,gpc,jlgpr)

