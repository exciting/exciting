
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
    ! *** draft ***
    complex(8) :: zt1
    real(8) :: t1,t2,s(3,3),vgi(3),vgj(3),v1(3),v(3),sn,cn
    integer :: igqi,igqj,ivi(3),ivj(3),isym,lspl
    zt=(1.d0,0.d0)
    ! consider symmetries
    if (nkpt.eq.nkptnr) then
       call tdzoutpr2(n,n,zt,xou,xou,you)
       call tdzoutpr2(n,n,zt,xuo,xuo,yuo)
    else
       nsym=nsymcrysstr(iqcu)
       ! *** simple loops for the moment ***
       do i=1,n
          do j=1,n
             do k=1,nsym
                ! symmetry element
                isym=scmapstr(k,ik)
                ! point group element
                lspl=lsplsymc(isym)
                ! rotation matrix in lattice coordinates
                s(:,:)=dble(symlat(:,:,lspl))
                ! first G-vector
                vgi(:)=vgql(:,i,iq)
                ! second G-vector
                vgj(:)=vgql(:,j,iq)
                v1(:)=vgi(:)-vgj(:)
                ! Note: apply symmetry to difference from left side
                call r3mv(s,v1,v)
                t1=dot_product(v,vtlsymc(:,isym))
                ! phase factor
                t2=cmplx(cos(t1),-sin(t1),8)
                ! index for first G-vector
                ivi(:)=nint(vgi(:)-vql(:,iq))
                ivi=matmul(symlat(:,:,lspl),ivi)
                ivi(:)=ivi(:)+ivscwrapq(:,isym,iq))
                igqi=
                ! index for second G-vector
                ivj(:)=nint(vgj(:)-vql(:,iq))
                ivj=matmul(symlat(:,:,lspl),ivj)
                ivj(:)=ivj(:)+ivscwrapq(:,isym,iq))
                
             end do
          end do
       end do
    end if
  end subroutine dfqoscbo

end module m_dfqoscbo
