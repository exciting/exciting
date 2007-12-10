
subroutine x0toasc
  use modmain
  use modxs
  use m_getx0
  use m_putx0
  use m_getunit
  use m_genfilname
  implicit none
  ! local variables
  character(*), parameter :: thisnam = 'x0toasc'
  character(256) :: filnam, filnama
  integer :: n,iq,iw,iproc,un
  real(8) :: vkloff_save(3)
  complex(8), allocatable :: chi0(:,:),chi0wg(:,:,:),chi0hd(:)
  logical :: tq0

  if (calledxs.eq.1) call init0

  ! initialise universal variables
  call init1

  ! save Gamma-point variables
  call tdsave0

  ! initialize q-point set
  call init2xs

  ! save k-point offset
  vkloff_save = vkloff

  ! loop over q-points
  do iq = 1, nqpt
     tq0 = tq1gamma.and.(iq.eq.1)
     ! shift k-mesh by q-point
     vkloff(:)=qvkloff(:,iq)
     ! calculate k+q and G+k+q related variables
     call init1td
     ! size of local field effects
     n = ngq(iq)
     ! allocate
     allocate(chi0(n,n),chi0wg(n,2,3),chi0hd(3))

     ! filenames
     call genfilname(asc=.true.,basename='X0',bzsampl=bzsampl,acont=acont,&
          nar=.not.aresdf,iq=iq,filnam=filnama)
     call genfilname(basename='X0',bzsampl=bzsampl,acont=acont,&
          nar=.not.aresdf,iq=iq,filnam=filnam)

     ! open file to write ASCI
     call getunit(un)
     open(unit=un,file=trim(filnama),form='formatted', &
          action='write',status='replace')

     do iw=1,nwdf
        ! read from binary file
        call getx0(tq0,iq,iw,trim(filnam),'',chi0,chi0wg,chi0hd)
        ! write to ASCII file
        if (tq0) then
           write(un,*) ngq(iq),vql(:,iq),chi0,chi0wg,chi0hd
        else
           write(un,*) ngq(iq),vql(:,iq),chi0
        end if
     end do

     ! close file
     close(un)

     deallocate(chi0,chi0wg,chi0hd)
     write(unitout,'(a,i8)') 'Info('//thisnam//'): Kohn Sham response &
          &function gathered for q-point:',iq
  end do

  ! restore offset
  vkloff = vkloff_save
  call genfilname(setfilext=.true.)

end subroutine x0toasc
