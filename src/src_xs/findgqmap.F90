
! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine findgqmap(iq,iqr,nsc,sc,ivgsc,nmax,n,isc,isci,ivgu,igqmap)
  use modmain
  use modxs
  implicit none
  ! arguments
  integer, intent(in) :: iq,iqr,nsc,sc(maxsymcrys),ivgsc(3,maxsymcrys),nmax,n
  integer, intent(out) :: isc,isci,ivgu(3),igqmap(nmax)
  ! local variables
  real(8) :: vqr(3),v2(3),t1
  integer :: iqrnr,j,isym,isymi,lspl,lspli,iv(3),ivg1(3),igq1
  integer, external :: iplocnr
  ! find map from G-vectors to rotated G-vectors
  iqrnr=iplocnr(ivqr(1,iqr),ngridq)
  vqr(:)=vqlr(:,iqr)
  do j=1,nsc
     isym=sc(j)
     lspl=lsplsymc(isym)
     isymi=scimap(isym)
     lspli=lsplsymc(isymi)
     do igq1=1,n
        ivg1(:)=ivg(:,igqig(igq1,iq))
        ! G1 = si^-1 * ( G + G_s ) , where si is the inverse of s
        iv=matmul(transpose(symlat(:,:,lspli)),ivg1+ivgsc(:,j))
        ! |G1 + q|
        v2=matmul(bvec,iv+vqr)
        t1=sqrt(sum(v2**2))
!!$write(*,'(a,5i6,3g18.10,3x,3g18.10)') 'reduce:',iq,j,isym,&
!!$     igq1,ivgigq(iv(1),iv(2),iv(3),iqrnr),vgql(:,igq1,iq)- &
!!$     matmul(transpose(symlat(:,:,lspl)),vqr+iv)
        if ((n.gt.1).and.(t1.gt.gqmax)) then
           write(*,*) '*** need one more symmetry operation'
           goto 10
        end if
        ! locate G1 + q in G+q-vector set
        igqmap(igq1)=ivgigq(iv(1),iv(2),iv(3),iqrnr)
        if (igqmap(igq1).le.0) then
           write(*,*)
           write(*,'("Error(findgqmap): failed to map rotated G-vector")')
           write(*,'(" non-reduced q-point                    :",i8)') iq
           write(*,'(" reduced q-point                        :",i8)') iqr
           write(*,'(" reduced q-point in non-reduced set     :",i8)') iqrnr
           write(*,'(" G+q-vector index (non-reduced q-point) :",i8)') igq1
           write(*,'(" rotated G-vector                       :",3i8)') iv
           write(*,*)
           call terminate
        end if
        ! end loop over G
     end do
     ! store G1 vector
     ivgu(:)=ivgsc(:,j)
     isc=isym
     isci=isymi
     goto 20
10   continue
     ! end loop over symmetry operations
  end do
  write(*,*)
  write(*,'("Error(findgqmap): failed to reduce q-point: ",i8)') iq
  write(*,*)
  call terminate
20 continue
end subroutine findgqmap
