


! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_sccli
  implicit none
contains


subroutine gensccli(s, lb, ub, emat13, emat24, sci, scimat)
    use indices
    implicit none
    ! arguments
    real(8), intent(in) :: s
    integer, intent(in) :: lb(4), ub(4)
    complex(8), intent(in) :: emat13(:, :), emat24(:, :), sci(:, :)
    complex(8), intent(out) :: scimat(:, :)
    ! local variables    
    integer :: siz(4), j12, j34, j13, j24, n12, n34, n13, n24
    complex(8), allocatable :: scm(:, :)
    siz(:)=ub(:)-lb(:)+1
    n12=siz(1)*siz(2); n34=siz(3)*siz(4); n13=siz(1)*siz(3); n24=siz(2)*siz(4)
    allocate(scm(n13, n24))
    ! carry out the double summation
    scm(:, :)=s*matmul(transpose(conjg(emat13)), matmul(sci, emat24))
    do j13=1, n13
       do j24=1, n24
	  call idxt1324(j12, j34, j13, j24, lb, ub, inv=.true.)
	  scimat(j12, j34)=scm(j13, j24)
       end do
    end do
    deallocate(scm)
  end subroutine gensccli

  ! tests ----------------------------------------------------------------------


subroutine test_sccli
    implicit none
    integer, parameter :: n1=3, n2=4, n3=3, n4=4
    integer, parameter :: lb(4)=(/1, 1, 1, 1/), ub(4)=(/n1, n2, n3, n4/)
    integer, parameter :: n=9
    real(8), parameter :: s=1.d0
    complex(8), parameter :: zi=(0.d0, 1.d0), zzero=(0.d0, 0.d0)
    integer :: g1, g2, i, j, i1, i2, i3, i4, i12, i34, i13, i24
    complex(8) :: emat13(n, n1*n3), emat24(n, n2*n4)
    complex(8) :: em13(n, n1, n3), em24(n, n2, n4)
    complex(8) :: sci(n, n), scimat(n1*n2, n3*n4), scimat2(n1*n2, n3*n4)
    complex(8) :: sc(n1, n2, n3, n4), zt, zt2

!!$    ! version 1 of assigning arrays
!!$    emat13(:,:)=zzero; emat24(:,:)=zzero; sci(:,:)=zzero
!!$    emat13(3:7,2:5)=7.d0+zi*2.d0
!!$    emat24(1:4,7:15)=-3.d0+zi*1.d0
!!$    sci(2:5,3:8)=-0.9d0+zi*0.4d0

    ! version 2 of assigning arrays
    call random_seed
    call zrandom_number(emat13)
    call zrandom_number(emat24)
    call zrandom_number(sci)

    call gensccli(s, (/1, 1, 1, 1/), (/n1, n2, n3, n4/), emat13, emat24, sci, scimat)

    do j=1, n3*n4
       do i=1, n1*n2
	  write(*, '(2i5, 2g18.10)') i, j, scimat(i, j)
       end do
    end do

    !*************************************************************************

    i13=0
    do i3=1, n3
       do i1=1, n1
	  i13=i13+1
	  em13(:, i1, i3)=emat13(:, i13)	  
       end do
    end do
    i24=0
    do i4=1, n4
       do i2=1, n2
	  i24=i24+1
	  em24(:, i2, i4)=emat24(:, i24)	  
       end do
    end do

    sc(:, :, :, :)=zzero
    i34=0
    do i4=1, n4
       do i3=1, n3
	  i34=i34+1
	  i12=0
	  do i2=1, n2
	     do i1=1, n1
		i12=i12+1
		do g1=1, n
		   do g2=1, n
		      sc(i1, i2, i3, i4) = sc(i1, i2, i3, i4)+ &
			   conjg(em13(g1, i1, i3)) * sci(g1, g2) * em24(g2, i2, i4)
		   end do
		end do
		scimat2(i12, i34)=sc(i1, i2, i3, i4)
	     end do
	  end do
       end do
    end do

    do j=1, n3*n4
       do i=1, n1*n2
	  zt=scimat(i, j); zt2=scimat2(i, j)
	  write(*, '(2i5, 2g18.8, 3x, 2g18.8, 3x, g18.8)') i, j, zt, zt2, abs(zt-zt2)
       end do
    end do

    write(*, *) 'deviation:', sum(abs(scimat-scimat2))

  end subroutine test_sccli


subroutine zrandom_number(z)
    implicit none
    complex(8), intent(out) :: z(:, :)
    real(8) :: x(size(z, 1), size(z, 2)), y(size(z, 1), size(z, 2))
    call random_number(x)
    call random_number(y)
    z=cmplx(x, y, 8)
  end subroutine zrandom_number

end module m_sccli
