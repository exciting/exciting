subroutine writederived(iqmt, eps, nw, w)
  use modmpi
  use modxs, only: unitout, fneps, fnloss, fnsigma, fnsumrules
  use modinput, only: input
  use m_genfilname
  use m_genloss
  use m_gensigma
  !use m_gensumrls
  use m_writeeps
  use m_writeloss
  use m_writesigma
  use m_write_hdf5
  use m_writesumrls

  implicit none

  ! Arguments
  integer(4), intent(in) :: nw, iqmt
  real(8), intent(in) :: w(nw)
  complex(8), intent(in) :: eps(3,3,nw)

  ! Local variables
  integer(4) :: o1, o2, io1, io2, ol, ou
  integer(4) :: optvec(3)
  real(8) :: loss(3, 3, nw)
  complex(8) :: sigma(nw)
  !real(8) :: sumrls(3)
  logical :: foff
  character(256) :: tdastring, bsetypestring, scrtypestring
  character(256) :: syscommand, epsilondir, lossdir, sigmadir
  character(256) :: gname

  character(*), parameter :: thisname = "writederived"

  write(unitout, '("Info(",a,"): Writing eps_m derived quantities.")') trim(thisname)

  epsilondir='EPSILON'
  lossdir='LOSS'
  sigmadir='SIGMA'
  syscommand = 'test ! -e '//trim(adjustl(epsilondir))//' && mkdir '//trim(adjustl(epsilondir))
  call system(trim(adjustl(syscommand)))
  syscommand = 'test ! -e '//trim(adjustl(lossdir))//' && mkdir '//trim(adjustl(lossdir))
  call system(trim(adjustl(syscommand)))
  syscommand = 'test ! -e '//trim(adjustl(sigmadir))//' && mkdir '//trim(adjustl(sigmadir))
  call system(trim(adjustl(syscommand)))

  if(input%xs%bse%coupling) then
    tdastring=''
  else
    tdastring="-TDA"
  end if

  bsetypestring = '-'//trim(input%xs%bse%bsetype)//trim(tdastring)
  scrtypestring = '-'//trim(input%xs%screening%screentype)

  ! Generate loss function as inverted dielectric tensor
  ! or 1/epsm for finite q
  loss = 0.0d0
  if(iqmt == 1) then 
    call genloss(eps, loss, 3)
  else
    call genloss(eps, loss, 1)
  end if

  ! Calculate derived quantities and write them to disk
  if(iqmt == 1) then 
    io1=1
    io2=3
    foff = input%xs%dfoffdiag
  else
    io1=1
    io2=1
    foff =.false.
  end if
  ! Generate filename for hdf5 datablock
  gname="spectra"//trim(bsetypestring)//trim(scrtypestring)

  !Write optical functions to hdf5
  call write_spectra_hdf5(iqmt, foff, w, eps, loss, sigma, gname)
  do o1 = io1, io2

    if(foff) then
      ol = 1
      ou = 3
    else
      ol = o1
      ou = o1
    end if

    do o2 = ol, ou

      optvec = (/ o1, o2, 0 /)

      ! Generate File names for resulting quantities
      if(iqmt == 1) then 
        call genfilname(basename='EPSILON', tq0=.true., oc1=o1, oc2=o2,&
          & bsetype=trim(bsetypestring), scrtype=trim(scrtypestring),&
          & nar= .not. input%xs%bse%aresbse, filnam=fneps)

        call genfilname(basename='LOSS', tq0=.true., oc1=o1, oc2=o2,&
          & bsetype=trim(bsetypestring), scrtype=trim(scrtypestring),&
          & nar= .not. input%xs%bse%aresbse, filnam=fnloss)

        call genfilname(basename='SIGMA', tq0=.true., oc1=o1, oc2=o2,&
          & bsetype=trim(bsetypestring), scrtype=trim(scrtypestring),&
          & nar= .not. input%xs%bse%aresbse, filnam=fnsigma)

       ! call genfilname(basename='SUMRULES', tq0=.true., oc1=o1, oc2=o2,&
       !   & bsetype=trim(bsetypestring), scrtype=trim(scrtypestring),&
       !   & nar= .not. input%xs%bse%aresbse, filnam=fnsumrules)
      else
        call genfilname(basename='EPSILON', iqmt=iqmt,&
          & bsetype=trim(bsetypestring), scrtype=trim(scrtypestring),&
          & nar= .not. input%xs%bse%aresbse, filnam=fneps)

        call genfilname(basename='LOSS', iqmt=iqmt,&
          & bsetype=trim(bsetypestring), scrtype=trim(scrtypestring),&
          & nar= .not. input%xs%bse%aresbse, filnam=fnloss)

        call genfilname(basename='SIGMA', iqmt=iqmt,&
          & bsetype=trim(bsetypestring), scrtype=trim(scrtypestring),&
          & nar= .not. input%xs%bse%aresbse, filnam=fnsigma)

       ! call genfilname(basename='SUMRULES', iqmt=iqmt,&
       !   & bsetype=trim(bsetypestring), scrtype=trim(scrtypestring),&
       !   & nar= .not. input%xs%bse%aresbse, filnam=fnsumrules)
      end if
      fneps=trim(epsilondir)//'/'//trim(fneps)
      fnloss=trim(lossdir)//'/'//trim(fnloss)
      fnsigma=trim(sigmadir)//'/'//trim(fnsigma)


      ! Generate optical functions
      call gensigma(w, eps(o1,o2,:), optvec(1:2), sigma)
      !call gensumrls(w, eps(o1,o2,:), sumrls)
      
      ! Write optical functions to file
      call writeeps(iqmt, o1, o2, w, eps(o1,o2,:), trim(fneps)) ! iqmt not used
      call writeloss(iqmt, w, loss(o1, o2, :), trim(fnloss))
      call writesigma(iqmt, w, sigma, trim(fnsigma))  ! iqmt not used
      !call writesumrls(iqmt, sumrls, trim(fnsumrules)) ! iqmt not used

    ! End loop over optical components
    end do
  end do

  write(unitout, '("  Derived quantities written.")')
end subroutine writederived
