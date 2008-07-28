
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gensymcmut(eps,maxsymcrys,nsymcrys,symlatc,lsplsymc,vtlsymc,scmut, &
     tabel)
! !DESCRIPTION:
!   Sets up the group multiplication table. The table is checked for consistency
!   in a way that it is required that every elements occurrs once and only once
!   in each row and column of the table. The first row and colmuns must consist
!   of the indentity since the first symmetry element is the identity by
!   convention.
!
! !REVISION HISTORY:
!   Created July 2008 (Sagmeister)
!EOP
!BOC
  implicit none
  ! arguments
  real(8), intent(in) :: eps
  integer, intent(in) :: maxsymcrys,nsymcrys
  real(8), intent(in) :: symlatc(3,3,48)
  integer, intent(in) :: lsplsymc(nsymcrys)
  real(8), intent(in) :: vtlsymc(3,maxsymcrys)
  integer, intent(out) :: scmut(maxsymcrys,maxsymcrys)
  logical, intent(out) :: tabel
  ! local variables
  integer, parameter :: maxsymlat=48
  integer :: i,asymlat,isymlat,jsymlat,isym,jsym,lspli,lsplj,iv(3)
  integer :: symlatmut(maxsymlat,maxsymlat),doner(maxsymlat),donec(maxsymlat)
  real(8) :: c(3,3),ct(3,3),s(3,3),si(3,3),sj(3,3),vtli(3),vtlj(3)
  scmut(:,:)=0
  do isymlat=1,maxsymlat
     si(:,:)=dble(symlatc(:,:,isymlat))
     do jsymlat=1,maxsymlat
        sj(:,:)=dble(symlatc(:,:,jsymlat))
!        call r3mm(si,sj,c)
        c=matmul(si,sj)
        do asymlat=1,maxsymlat
           s(:,:)=dble(symlatc(:,:,asymlat))
           ct(:,:)=c(:,:)-s(:,:)
           if (sum(abs(ct)).lt.eps) then
              ! add element to multiplication table
              scmut(isymlat,jsymlat)=asymlat
              exit
           end if
        end do
        write(*,*) isymlat,jsymlat,scmut(isymlat,jsymlat)
     end do
  end do
  ! check multiplication table for consistency
  do isymlat=1,maxsymlat
     donec(:)=0
     doner(:)=0
     do jsymlat=1,maxsymlat
        doner(scmut(isymlat,jsymlat))=doner(scmut(isymlat,jsymlat))+1
        donec(scmut(jsymlat,isymlat))=donec(scmut(jsymlat,isymlat))+1
     end do
     do jsymlat=1,maxsymlat
        if (doner(jsymlat).ne.1) then  
           write(*,*)
           write(*,'("Error(gensymcmut): error in multiplication table")')
           write(*,'(" row number    : ",i6)') isymlat
           write(*,'(" column number : ",i6)') jsymlat
           write(*,'(" multiple occurrence : ",i6)') doner(jsymlat)
           write(*,*)
           stop
        end if
        if (donec(jsymlat).ne.1) then  
           write(*,*)
           write(*,'("Error(gensymcmut): error in multiplication table")')
           write(*,'(" row number    : ",i6)') jsymlat
           write(*,'(" column number : ",i6)') isymlat
           write(*,'(" * multiple occurrence : ",i6)') donec(jsymlat)
           write(*,*)
           stop
        end if
     end do
  end do
  if (any(scmut(:,1).ne.1).or.any(scmut(1,:).ne.1)) then
     write(*,*)
     write(*,'("Error(gensymcmut): error in multiplication table")')
     write(*,'(" first row or column does not consist of identity")')
     write(*,*)
write(*,*) 'col:',scmut(
     stop
  end if
  ! check if group is Abelian
  tabel=.false.
  if (all(scmut-transpose(scmut).eq.0)) tabel=.true.
end subroutine gensymcmut
