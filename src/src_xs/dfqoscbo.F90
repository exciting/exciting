
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_dfqoscbo
  implicit none
contains

  subroutine dfqoscbo(iq,ik,n,xou,xuo,you,yuo)
    use modmain
    use modxs
    use m_tdzoutpr2
    implicit none
    ! arguments
    integer, intent(in) :: iq,ik,n
    complex(8), intent(in) :: xou(:),xuo(:)
    complex(8), intent(out) :: you(:,:),yuo(:,:)
    ! local variables
    complex(8) :: zt
    integer :: i,j,k,nsym
    complex(8) :: zt1
    real(8) :: t1,t2,s(3,3),vgi(3),vgj(3),v1(3),v(3)
    integer :: igqi,igqj,ivi(3),ivj(3),isym,lspl
    zt=(1.d0,0.d0)
    ! consider symmetries
    if (nkpt.eq.nkptnr) then
       call tdzoutpr2(n,n,zt,xou,xou,you)
       call tdzoutpr2(n,n,zt,xuo,xuo,yuo)
    else
       nsym=nsymcrysstr(ik)
       ! *** simple loops for the moment ***
       do i=1,n
          ! first G-vector
          vgi(:)=vgql(:,i,iq)
          do j=1,n
             ! second G-vector
             vgj(:)=vgql(:,j,iq)
             ! difference of G-vectors
             v1(:)=vgi(:)-vgj(:)
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
                ! index for second G-vector
                ivj(:)=ivg(:,igqig(j,iq))
                ivj=matmul(symlat(:,:,lspl),ivj)
                ivj(:)=ivj(:)+ivscwrapq(:,isym,iq)
write(*,*) 'iq,ik,i,j,k',iq,ik,i,j,k
write(*,*) 'v1',v1
write(*,*) 'ivi',ivi
write(*,*) 'ivj',ivj
write(*,*) 'symlat',symlat(1,:,lspl)
write(*,*) 'symlat',symlat(2,:,lspl)
write(*,*) 'symlat',symlat(3,:,lspl)
write(*,*) 'v',v
write(*,*)

                igqj=ivgigq(ivj(1),ivj(2),ivj(3),iq)
                ! update oscillators
                you(i,j)=you(i,j)+t2*xou(igqi)*conjg(xou(igqj))
                yuo(i,j)=yuo(i,j)+t2*xuo(igqi)*conjg(xuo(igqj))                
             end do
          end do
       end do
    end if
  end subroutine dfqoscbo

end module m_dfqoscbo
