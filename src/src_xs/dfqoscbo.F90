
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_dfqoscbo
  implicit none
contains

!BOP
! !ROUTINE: gengvec
! !INTERFACE:
  subroutine dfqoscbo(iq,ik,ni,n,xou,xuo,you,yuo)
! !USES:
    use modmain
    use modxs
    use m_tdzoutpr2
! !DESCRIPTION:
!   Generates the oscillators for the body of the dielectric function.
!   Note that the symmetries are applied from the left side in Cartesian
!   coordinates, so that the resulting operation in lattice coordinates is
!   an application of the inverse operation from the right side:
!   $(sG)^R = (s^{-1,L})^T G^R$ where the $R$ superscript denotes a vector
!   containing the coordinates with respect to the reciprocal frame and a $L$
!   superscript denotes a reference to the frame of the real lattice.
!
! !REVISION HISTORY:
!   Created December 2007 (Sagmeister)
!EOP
!BOC
    implicit none
    ! arguments
    integer, intent(in) :: iq,ik,ni,n
    complex(8), intent(in) :: xou(:),xuo(:)
    complex(8), intent(out) :: you(:,:),yuo(:,:)
    ! local variables
    complex(8) :: zt
    integer :: i,j,k,nsym
    complex(8) :: zt1
    real(8) :: t1,t2,s(3,3),vgi(3),vgj(3),v1(3),v(3),vr(3),vrs(3)
    integer :: igqi,igqj,ivi(3),ivj(3),isym,jsym,lspl,   lsplj
    zt=(1.d0,0.d0)
    ! consider symmetries
    if (nkpt.eq.nkptnr) then
       call tdzoutpr2(ni,n,n,zt,xou,xou,you)
       call tdzoutpr2(ni,n,n,zt,xuo,xuo,yuo)
    else
       nsym=nsymcrysstr(ik)
       ! *** simple loops for the moment ***
       do i=ni,n
          ! first G-vector
          vgi(:)=vgql(:,i,iq)
          do j=ni,n
             ! second G-vector
             vgj(:)=vgql(:,j,iq)
             ! difference of G-vectors
             v1(:)=vgi(:)-vgj(:)
             do k=1,nsym
                ! symmetry element
                isym=scmapstr(k,ik)
                ! inverse of symmetry element
                jsym=scimap(isym)
                ! point group element
                lspl=lsplsymc(jsym)
                lsplj=lsplsymc(isym)
!!$write(*,*) 'iq,ik,i,j,k',iq,ik,i,j,k
!!$write(*,*) 'jsym,lspl,symlat/symlatc',jsym,lspl
!!$write(*,*) symlat(1,:,lspl)
!!$write(*,*) symlat(2,:,lspl)
!!$write(*,*) symlat(3,:,lspl)
!!$write(*,*) symlatc(1,:,lspl)
!!$write(*,*) symlatc(2,:,lspl)
!!$write(*,*) symlatc(3,:,lspl)
                ! rotation matrix in lattice coordinates
                s(:,:)=dble(symlat(:,:,lspl))
                ! apply symmetry to difference
                call r3mtv(s,v1,v)
                t1=twopi*dot_product(v,vtlsymc(:,jsym))
                ! phase factor
                zt1=cmplx(cos(t1),sin(t1),8)
                ! index for first G-vector
                ivi(:)=ivg(:,igqig(i,iq))
!!$write(*,*) 'ivi             ',ivi
!!$vr(:)=ivi(1)*bvec(:,1)+ivi(2)*bvec(:,2)+ivi(3)*bvec(:,3)
!!$write(*,'(a,3f12.4)') 'vi              ',vr
                ivi=matmul(ivi,symlat(:,:,lspl))
!!$write(*,*) 'symlat(ivi)     ',ivi
!!$vr(:)=ivi(1)*bvec(:,1)+ivi(2)*bvec(:,2)+ivi(3)*bvec(:,3)
!!$write(*,'(a,3f12.4)') 'symlat(vi)      ',vr
                ivi(:)=ivi(:)+ivscwrapq(:,jsym,iq)
!!$write(*,*) 'symlat(ivi)+wrap',ivi
!!$vr(:)=ivi(1)*bvec(:,1)+ivi(2)*bvec(:,2)+ivi(3)*bvec(:,3)
!!$write(*,'(a,3f12.4)') 'symlat(vi)+wrap ',vr
                igqi=ivgigq(ivi(1),ivi(2),ivi(3),iq)
                ! index for second G-vector
                ivj(:)=ivg(:,igqig(j,iq))
!!$write(*,*) 'ivj             ',ivj
!!$vr(:)=ivj(1)*bvec(:,1)+ivj(2)*bvec(:,2)+ivj(3)*bvec(:,3)
!!$vrs=matmul(symlatc(:,:,lsplj),vr)
!!$write(*,'(a,3f12.4)') 'vj              ',vr
                ivj=matmul(ivj,symlat(:,:,lspl))
!!$write(*,*) 'symlat(ivj)     ',ivj
!!$vr(:)=ivj(1)*bvec(:,1)+ivj(2)*bvec(:,2)+ivj(3)*bvec(:,3)
!!$write(*,'(a,3f12.4,3x,3f12.4)') 'symlat(vj)/Cart. ',vr,vrs
                ivj(:)=ivj(:)+ivscwrapq(:,jsym,iq)
!!$write(*,*) 'symlat(ivj)+wrap',ivj
!!$vr(:)=ivj(1)*bvec(:,1)+ivj(2)*bvec(:,2)+ivj(3)*bvec(:,3)
!!$write(*,'(a,3f12.4)') 'symlat(vj)+wrap ',vr
                igqj=ivgigq(ivj(1),ivj(2),ivj(3),iq)
                ! update oscillators
                you(i,j)=you(i,j)+zt1*xou(igqi)*conjg(xou(igqj))
                yuo(i,j)=yuo(i,j)+zt1*xuo(igqi)*conjg(xuo(igqj))
!!$write(*,*) '==================================================================='
             end do
          end do
       end do
    end if
  end subroutine dfqoscbo
!EOC

end module m_dfqoscbo
