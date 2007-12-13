
! Copyright (C) 2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine findsymi(epslat,maxsymcrys,nsymcrys,symlat,lsplsymc,vtlsymc,scimap)
  implicit none
  ! arguments
  real(8), intent(in) :: epslat
  integer, intent(in) :: symlat(3,3,48)
  integer, intent(in) :: maxsymcrys,nsymcrys
  integer, intent(in) :: lsplsymc(nsymcrys)
  real(8), intent(in) :: vtlsymc(3,maxsymcrys)
  integer, intent(out) :: scimap(maxsymcrys)
  ! local variables
  real(8) :: c(3,3),si(3,3),sj(3,3),vtl(3)
  integer :: i,isym,jsym,lspli,lsplj,iv(3)

!!$  real(8) :: s1(3,3)

  scimap(:)=0
  do isym=1,nsymcrys
     lspli=lsplsymc(isym)
     si(:,:)=dble(symlat(:,:,lspli))
     do jsym=1,nsymcrys
        lsplj=lsplsymc(jsym)
        sj(:,:)=dble(symlat(:,:,lsplj))
        ! translation
        vtl(:)=vtlsymc(:,jsym)
        call r3mv(si,vtl,vtl)
        vtl(:)=vtl(:)+vtlsymc(:,isym)
        call r3frac(epslat,vtl,iv)
        ! rotation
        call r3mm(si,sj,c)
        ! subract unit matrix
        forall (i=1:3)
           c(i,i)=c(i,i)-1.d0
        end forall
        if ((sum(vtl).lt.epslat).and.(sum(abs(c)).lt.epslat)) then
           ! isym is inverse of jsym
           scimap(isym)=jsym
           goto 10
        end if
     end do
10   continue
  end do

!!$  ! Test to see if ((alpha^-1)_latt == (alpha_latt)^-1
!!$  do isym=1,nsymcrys
!!$     jsym=scimap(isym)
!!$     lspli=lsplsymc(isym)
!!$     lsplj=lsplsymc(jsym)
!!$     si(:,:)=dble(symlat(:,:,lspli))
!!$     sj(:,:)=dble(symlat(:,:,lsplj))
!!$     call r3minv(si,s1)
!!$     write(*,*) 'isym,jsym,lspli,lsplj,diff',isym,jsym,lspli,lsplj,&
!!$          sum(abs(s1-sj))
!!$  end do



end subroutine findsymi
