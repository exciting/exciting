
subroutine idf
  use modmain
  use modtddft
  use modfxcifc
  use modmpi
  use m_idfq
  use m_tdlinopt
  use m_getunit
  implicit none
  ! local variables
  character(*), parameter :: thisnam = 'idf'
  integer, save :: called
  data called /0/
  integer :: iq, un, wi,wf

  called=called+1
  call getunit(un)

  ! initialise universal variables
  if (called.eq.1) call init0
  call init1

  ! initialize q-point set
  call init2xs

  ! w-point interval for process
  wpari=firstofset(rank,nwdf)
  wparf=lastofset(rank,nwdf)

  write(unitout,'("Exchange-correlation kernel type :",i4)') fxctype
  write(unitout,'("  ",a)') trim(fxcdescr)

  ! loop over q-points
  do iq = 1, nqpt
     ! call for q-point
     call idfq(iq)
     write(unitout,'(a,i8)') 'Info('//thisnam//'): inverse dielectric &
          &function finished for q-point:',iq
  end do

  call barrier(rank=rank,procs=procs,un=un,async=0,string='.barrier')

  if ((procs > 1).and.(rank == 0)) then
     call idfgather
        write(unitout,'(a)') 'Info('//thisnam//'): inverse dielectric &
             &function gathered for q-point:'
  end if

  call barrier(rank=rank,procs=procs,un=un,async=0,string='.barrier')

  if (rank == 0) then
     do iq=1,nqpt
        ! call for q-point
        call tdlinopt(iq)
        write(unitout,'(a,i8)') 'Info('//thisnam//'): TDDFT linear optics &
             &finished for q-point:',iq
     end do
  end if

  call barrier(rank=rank,procs=procs,un=un,async=0,string='.barrier')

  write(unitout,'(a)') "Info("//trim(thisnam)//"): TDDFT linear optics &
       &finished"

end subroutine idf
