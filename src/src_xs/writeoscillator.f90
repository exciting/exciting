subroutine writeoscillator(hamsize, nexc, egap, beval, oszs)
  use modinput, only: input
  use modxs, only: bsed, escale
#ifdef DGRID
  use modxs, only: dgrid, iksubpt
#endif
  use m_genfilname
  use m_getunit
  implicit none

  ! I/O
  integer(4), intent(in) :: hamsize, nexc
  real(8), intent(in) :: egap, beval(hamsize)
  complex(8), intent(in) :: oszs(nexc,3)

  ! Local
  integer(4) :: o1, lambda, unexc
  character(256) :: fnexc, frmt
#ifdef DGRID
  character(256) :: dgrid_dotext
#endif

  ! Loop over optical components
  do o1=1,3

#ifdef DGRID
    ! Stk: Add case of double grid
    if(dgrid) then 

      write(dgrid_dotext, '("_SG", i3.3, ".OUT")') iksubpt

      call genfilname(basename='EXCITON', tq0=.true., oc1=o1, oc2=o1,&
        & bsetype=input%xs%bse%bsetype, scrtype=input%xs%screening%screentype,&
        & nar= .not. input%xs%bse%aresbse, dotext=dgrid_dotext, filnam=fnexc)

    else

      call genfilname(basename='EXCITON', tq0=.true., oc1=o1, oc2=o1,&
        & bsetype=input%xs%bse%bsetype, scrtype=input%xs%screening%screentype,&
        & nar= .not. input%xs%bse%aresbse, filnam=fnexc)
      
    endif
#else
    call genfilname(basename='EXCITON', tq0=.true., oc1=o1, oc2=o1,&
      & bsetype=input%xs%bse%bsetype, scrtype=input%xs%screening%screentype,&
      & nar= .not. input%xs%bse%aresbse, filnam=fnexc)
#endif

    ! Write out exciton energies and oscillator strengths
    call getunit(unexc)
    open(unexc, file=fnexc, form='formatted', action='write', status='replace')
    if(input%xs%tevout) write(unexc, '("# All energies are given in electron volts")')
    write(unexc, '("# E_bsegap : ", SP, E23.16)') egap * escale
    frmt='(a1,a7,1x,a23,1x,a23,1x,a23,1x,a23,1x,a23)'
    write(unexc, frmt) "#", "Nr.",&
      & "E",&
      & "E-E_bsegap",&
      & "|Osc. Str.|",&
      & "Re(Osc. Str.)",&
      & "Im(Osc. Str.)"
    frmt='(I8,1x,E23.16,1x,E23.16,1x,E23.16,1x,E23.16,1x,E23.16)'
    do lambda = 1, nexc
      write(unexc, frmt) lambda,&
        & (beval(lambda)+egap-dble(bsed))*escale,&
        & (beval(lambda)+dble(bsed))*escale,&
        & abs(oszs(lambda, o1)),&
        & dble(oszs(lambda, o1)),&
        & aimag(oszs(lambda, o1))
    end do
    close(unexc)

  end do

end subroutine writeoscillator
