module m_writeoscillator
  implicit none

  contains

    subroutine writeoscillator(hamsize, nexc, egap, evalre, oszstrr, evalim, oszstra)
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
      real(8), intent(in) :: egap, evalre(hamsize)
      complex(8), intent(in) :: oszstrr(nexc,3)
      real(8), intent(in), optional :: evalim(hamsize)
      complex(8), intent(in), optional :: oszstra(nexc,3)

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

        if(present(oszstra) .and. present(evalim)) then

          frmt='(a1,a7,8(1x,a23))'
          write(unexc, frmt) "#", "Nr.",&
            & "Re(E)",&
            & "Im(E)",&
            & "|Osc. Str. Res.|",&
            & "Re(Osc. Str. Res.)",&
            & "Im(Osc. Str. Res.)",&
            & "|Osc. Str. Ares.|",&
            & "Re(Osc. Str. Ares.)",&
            & "Im(Osc. Str. Ares.)"
          frmt='(I8,8(1x,E23.16))'
          do lambda = 1, nexc
            write(unexc, frmt) lambda,&
              & (evalre(lambda))*escale,&
              & evalim(lambda)*escale,&
              & abs(oszstrr(lambda, o1)),&
              & dble(oszstrr(lambda, o1)),&
              & aimag(oszstrr(lambda, o1)),&
              & abs(oszstra(lambda, o1)),&
              & dble(oszstra(lambda, o1)),&
              & aimag(oszstra(lambda, o1))
          end do

        else if(present(oszstra)) then

          frmt='(a1,a7,1x,a23,1x,a23,1x,a23,1x,a23,1x,a23,1x,a23,1x,a23)'
          write(unexc, frmt) "#", "Nr.",&
            & "Re(E)",&
            & "|Osc. Str. Res.|",&
            & "Re(Osc. Str. Res.)",&
            & "Im(Osc. Str. Res.)",&
            & "|Osc. Str. Ares.|",&
            & "Re(Osc. Str. Ares.)",&
            & "Im(Osc. Str. Ares.)"
          frmt='(I8,1x,E23.16,1x,E23.16,1x,E23.16,&
            &1x,E23.16,1x,E23.16,1x,E23.16,1x,E23.16)'
          do lambda = 1, nexc
            write(unexc, frmt) lambda,&
              & (evalre(lambda))*escale,&
              & abs(oszstrr(lambda, o1)),&
              & dble(oszstrr(lambda, o1)),&
              & aimag(oszstrr(lambda, o1)),&
              & abs(oszstra(lambda, o1)),&
              & dble(oszstra(lambda, o1)),&
              & aimag(oszstra(lambda, o1))
          end do

        else

          frmt='(a1,a7,1x,a23,1x,a23,1x,a23,1x,a23,1x,a23)'
          write(unexc, frmt) "#", "Nr.",&
            & "E",&
            & "E-E_bsegap",&
            & "|Osc. Str. Res.|",&
            & "Re(Osc. Str. Res.)",&
            & "Im(Osc. Str. Res.)"
          frmt='(I8,1x,E23.16,1x,E23.16,1x,E23.16,1x,E23.16,1x,E23.16)'
          do lambda = 1, nexc
            write(unexc, frmt) lambda,&
              & (evalre(lambda)+egap-dble(bsed))*escale,&
              & (evalre(lambda)+dble(bsed))*escale,&
              & abs(oszstrr(lambda, o1)),&
              & dble(oszstrr(lambda, o1)),&
              & aimag(oszstrr(lambda, o1))
          end do

        end if

        close(unexc)

      end do

    end subroutine writeoscillator

end module m_writeoscillator
