subroutine writederived(iqmt, eps, nw, w)
  use modxs, only: fneps, fnloss, fnsigma, fnsumrules
  use modinput, only: input
  use m_genfilname
  use m_genloss
  use m_gensigma
  use m_gensumrls
  use m_writeeps
  use m_writeloss
  use m_writesigma
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
  real(8) :: sumrls(3)
  logical :: foff

  ! Generate loss function as inverted dielectric tensor
  if(iqmt == 0) then 
    call genloss(eps, loss, 3)
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
          & bsetype=input%xs%bse%bsetype, scrtype=input%xs%screening%screentype,&
          & nar= .not. input%xs%bse%aresbse, filnam=fneps)

        call genfilname(basename='LOSS', tq0=.true., oc1=o1, oc2=o2,&
          & bsetype=input%xs%bse%bsetype, scrtype=input%xs%screening%screentype,&
          & nar= .not. input%xs%bse%aresbse, filnam=fnloss)

        call genfilname(basename='SIGMA', tq0=.true., oc1=o1, oc2=o2,&
          & bsetype=input%xs%bse%bsetype, scrtype=input%xs%screening%screentype,&
          & nar= .not. input%xs%bse%aresbse, filnam=fnsigma)

        call genfilname(basename='SUMRULES', tq0=.true., oc1=o1, oc2=o2,&
          & bsetype=input%xs%bse%bsetype, scrtype=input%xs%screening%screentype,&
          & nar= .not. input%xs%bse%aresbse, filnam=fnsumrules)
      else
        call genfilname(basename='EPSILON', iqmt=iqmt,&
          & bsetype=input%xs%bse%bsetype, scrtype=input%xs%screening%screentype,&
          & nar= .not. input%xs%bse%aresbse, filnam=fneps)

        call genfilname(basename='SIGMA', iqmt=iqmt,&
          & bsetype=input%xs%bse%bsetype, scrtype=input%xs%screening%screentype,&
          & nar= .not. input%xs%bse%aresbse, filnam=fnsigma)

        call genfilname(basename='SUMRULES', iqmt=iqmt,&
          & bsetype=input%xs%bse%bsetype, scrtype=input%xs%screening%screentype,&
          & nar= .not. input%xs%bse%aresbse, filnam=fnsumrules)
      end if


      ! Generate optical functions
      call gensigma(w, eps(o1,o2,:), optvec(1:2), sigma)
      call gensumrls(w, eps(o1,o2,:), sumrls)

      ! Write optical functions to file
      call writeeps(iqmt, o1, o2, w, eps(o1,o2,:), trim(fneps)) ! iqmt not used
      ! For writeloss iqmt must be a G+q index.
      if(iqmt==1) then 
        call writeloss(iqmt, w, loss(o1, o2, :), trim(fnloss))
      end if
      call writesigma(iqmt, w, sigma, trim(fnsigma))  ! iqmt not used
      call writesumrls(iqmt, sumrls, trim(fnsumrules)) ! iqmt not used

    ! End loop over optical components
    end do
  end do
end subroutine writederived
