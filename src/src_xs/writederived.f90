subroutine writederived(iq, eps, nw, w)
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
  integer(4), intent(in) :: nw, iq
  real(8), intent(in) :: w(nw)
  complex(8), intent(in) :: eps(3,3,nw)

  ! Local variables
  integer(4) :: o1, o2, ol, ou
  integer(4) :: optvec(3)
  real(8) :: loss(3, 3, nw)
  complex(8) :: sigma(nw)
  real(8) :: sumrls(3)

  ! Generate loss function as inverted dielectric tensor
  call genloss(eps, loss, 3)

  ! Calculate derived quantities and write them to disk
  do o1 = 1, 3

    if(input%xs%dfoffdiag) then
      ol = 1
      ou = 3
    else
      ol = o1
      ou = o1
    end if

    do o2 = ol, ou

      optvec = (/ o1, o2, 0 /)

      ! Generate File names for resulting quantities
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


      ! Generate optical functions
      call gensigma(w, eps(o1,o2,:), optvec(1:2), sigma)
      call gensumrls(w, eps(o1,o2,:), sumrls)

      ! Write optical functions to file
      call writeeps(iq, o1, o2, w, eps(o1,o2,:), trim(fneps)) ! iq not used
      ! For writeloss iq must be a G+q index.
      call writeloss(iq, w, loss(o1, o2, :), trim(fnloss))
      call writesigma(iq, w, sigma, trim(fnsigma))  ! iq not used
      call writesumrls(iq, sumrls, trim(fnsumrules)) ! iq not used

    ! End loop over optical components
    end do
  end do
end subroutine writederived
