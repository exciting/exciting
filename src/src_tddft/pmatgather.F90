
subroutine pmatgather()
  use modmain
  use modtddft
  use modpar
  use m_filedel
  use m_getpmat
  use m_putpmat
  implicit none
  ! local variables
  character(*), parameter :: thisnam = 'pmatgather'
  character(256) :: filext_save
  complex(8), allocatable :: pm(:,:,:)
  integer :: ik,ikr,iproc,recl

  ! allocate matrix elements array
  allocate(pm(3,nstsv,nstsv))
  filext_save=trim(filext)

  ! file extension for q-point
  write(filext,'(".OUT")')
  do iproc=1,nproc
     write(filextp,'("_par",i3.3,".OUT")') iproc
     call getrange(iproc,nproc,nkpt,kpari,kparf)
     do ik=kpari,kparf
        ! momentum matrix elements
        call getpmat(ik,vkl,.false.,trim(fnpmat)//trim(filextp),pm)
        call putpmat(ik,.true.,'PMAT_TD'//trim(filext),pm)
     end do
  end do

  ! delete partial files
  do iproc=1,nproc
     write(filextp,'("_par",i3.3,".OUT")') iproc
     call filedel(trim(fnpmat)//trim(filextp))
  end do

  write(unitout,'(a,i8)') 'Info('//thisnam//'): Momentum matrix &
       &elements gathered.'

  ! deallocate
  deallocate(pm)
  filext=trim(filext_save)

end subroutine pmatgather
