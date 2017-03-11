subroutine energy2index(n, ns, evals, e1, e2, i1, i2)
  use modmpi
  integer(4), intent(in) :: n, ns
  real(8), intent(in) :: evals(n), e1, e2
  integer(4), intent(out) :: i1, i2

  integer(4) :: i
  logical :: f1, f2

  f1=.false.
  f2=.false.

  if(e2-e1 < 0.0d0) then 
    if(rank == 0) then 
      write(*,'("Error(energy2index): e2 < e1", 2E12.3)') e1, e2
    end if
    i1=1
    i2=1
    return
  end if

  if(evals(1) > e1) then 
    !write(*,'("Waring(energy2index): e1 smaller than smalles eigenvalue")')
    !write(*,*) "e1=",e1
    !write(*,*) "evals(1)=", evals(1)
    i1 = 1
    f1 = .true.
  end if
  if(evals(ns) < e2) then 
    !write(*,'("Waring(energy2index): e2 larger than largest eigenvalue")')
    !write(*,*) "e2=",e2
    !write(*,*) "evals(ns)=", evals(ns)
    i2 = ns
    f2 = .true.
  end if
  if(f2 .and. f1) return

  i1=0
  i2=1
  do i=1,ns
    if(evals(i) < e1) i1=i
    if(evals(i) < e2) i2=i
  end do
  i1=i1+1

end subroutine energy2index
