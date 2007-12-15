
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_dfqoscwg
  implicit none
contains

  subroutine dfqoscwg(iq,ik,sw,n,pou,puo,xou,xuo,you,yuo)
    use modmain
    use modxs
    implicit none
    ! arguments
    integer, intent(in) :: iq,ik,sw,n
    complex(8), intent(in) :: pou(3),puo(3),xou(:),xuo(:)
    complex(8), intent(out) :: you(:),yuo(:)
    ! local variables
    complex(8) :: pout,puot
    real(8) :: s(3,3)
    integer :: oc,i,k,nsym
    complex(8) :: zt1
    real(8) :: t1,t2,vgi(3),v1(3),v(3)
    integer :: igqi,ivi(3),isym,lspl
    ! *** check if these symmetry considerations are correct ***
    if (symwings) then
       s(:,:)=0.5d0*(symdfq0(:,:)+transpose(symdfq0))
       pout=zzero
       puot=zzero
       do i=1,3
          pout=pout+s(i,i)*pou(i)
          puot=pout+s(i,i)*puo(i)
       end do
    else
       oc=optcomp(1,1)
       pout=pou(oc)
       puot=puo(oc)
    end if

    ! <= version 0.9.114
!!$    if (sw.eq.1) then
!!$       ! first wing (G=0,G/=0)
!!$       you(:)=pout*conjg(xou(:))
!!$       yuo(:)=puot*conjg(xuo(:))
!!$    else
!!$       ! second wing (G/=0,G=0)
!!$       you(:)=xou(:)*conjg(pout)
!!$       yuo(:)=xuo(:)*conjg(puot)
!!$    end if


    ! consider symmetries
    if (nkpt.eq.nkptnr) then
       ! first wing (G=0,G/=0)
       you(:)=pout*conjg(xou(:))
       yuo(:)=puot*conjg(xuo(:))
    else
       nsym=nsymcrysstr(ik)
       ! *** simple loops for the moment ***
       do i=1,n
          ! first G-vector
          vgi(:)=vgql(:,i,iq)
          ! difference of G-vectors
          v1(:)=0.d0-vgi(:)
          do k=1,nsym
             ! symmetry element
             isym=scmapstr(k,ik)
             ! point group element
             lspl=lsplsymc(isym)
             ! rotation matrix in lattice coordinates
             s(:,:)=dble(symlat(:,:,lspl))
             ! Note: apply symmetry to difference from left side
             call r3mv(s,v1,v)
             t1=dot_product(v,vtlsymc(:,isym))
             ! phase factor
             t2=cmplx(cos(t1),sin(t1),8)
             ! index for first G-vector
             ivi(:)=ivg(:,igqig(i,iq))
             ivi=matmul(symlat(:,:,lspl),ivi)
             ivi(:)=ivi(:)+ivscwrapq(:,isym,iq)
             igqi=ivgigq(ivi(1),ivi(2),ivi(3),iq)
             ! update oscillators
             you(i)=you(i)+t2*pout*conjg(xou(igqi))
             yuo(i)=yuo(i)+t2*puot*conjg(xuo(igqi))                
          end do
       end do
    end if
    ! second wing (G/=0,G=0)
    if (sw.ne.1) then
       you(:)=conjg(you(:))
       yuo(:)=conjg(yuo(:))
    end if


  end subroutine dfqoscwg

end module m_dfqoscwg
