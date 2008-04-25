subroutine solvediislin(m,Pmatrix,Qmatrix,c)
  use diisinterfaces
  use sclcontroll,only:recalculate_preconditioner
  implicit none

  integer, intent(in)::m

  real(8), intent(inout)::Pmatrix(m+1,m+1),Qmatrix(m+1,m+1)
  real(8), intent(out)::c(m+1)
  integer:: ipiv(m+1),info=0,i,RANK,lwork=-1
  real(8),allocatable ::WORK (:)
  real(8)::tmp(1)
  Pmatrix(m+1,m+1)=0.0
  c(m+1)=1
  do i=1,m
     Pmatrix(i,m+1)=1.0
     Pmatrix(m+1,i)=1.0
     c(i)=0
  end do
  call  DGESV( m+1, 1,Pmatrix , m+1, IPIV, c, m+1, INFO )
  !call DGELSY( m+1, m+1, 1, Pmatrix, m+1, c, m+1, IPIV,0.0, RANK,&
  !          tmp, LWORK, INFO )
  !  LWORK=tmp(1)
  !  allocate(WORK(LWORK))
  !  INFO=1
  !call DGELSY( m+1, m+1, 1, Pmatrix, m+1, c, m+1, IPIV,1.0, RANK,&
  !              WORK, LWORK, INFO )            

  if (info.ne.0) then
     write(*,*)
     write(*,'("Error(solvediis):  failed")')
     write(*,'(" DGESV returned INFO = ",I8)') info
#ifdef DEBUG               
     write(775,*)(Pmatrix)
     write(776,*)(Qmatrix)
     ! stop
     recalculate_preconditioner=.true.
#endif
  end if

end subroutine solvediislin
