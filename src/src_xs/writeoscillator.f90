module m_writeoscillator
  use modinput, only: input
  use modxs, only: bsed, escale
#ifdef DGRID
  use modxs, only: dgrid, iksubpt
#endif
  use m_genfilname
  use m_getunit

  implicit none

  contains

    subroutine writeoscillator(hamsize, nexc, eshift, evalre, oszstrr,&
      & evalim, oszstra, sort)
      ! I/O
      integer(4), intent(in) :: hamsize, nexc
      real(8), intent(in) :: eshift
      real(8), intent(in) :: evalre(hamsize)
      complex(8), intent(in) :: oszstrr(nexc,3)
      real(8), intent(in), optional :: evalim(hamsize)
      complex(8), intent(in), optional :: oszstra(nexc,3)
      logical, intent(in), optional :: sort

      ! Local
      logical :: fsort
      integer(4) :: o1, lambda, unexc
      integer(4), allocatable :: idxsort(:)
      real(8) :: pm
      character(256) :: fnexc, frmt
#ifdef DGRID
      character(256) :: dgrid_dotext
#endif
      
      if(present(sort)) then 
        fsort = sort
      else
        fsort = .false.
      end if

      allocate(idxsort(nexc))
      if(fsort) then 
        call sortidx(nexc, evalre, idxsort)
        idxsort = cshift(idxsort, nexc/2)
      else
        do lambda = 1, nexc
          idxsort(lambda) = lambda
        end do
      end if

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
        write(unexc, '("# E_shift : ", SP, E23.16)') eshift * escale

        ! Full BSE
        if(present(oszstra) .and. present(evalim)) then

          frmt='(a1,a7,9(1x,a23))'
          write(unexc, frmt) "#", "Nr.",&
            & "Re(E)",&
            & "Re(E)-+shift",&
            & "Im(E)",&
            & "|Osc. Str. Res.|",&
            & "Re(Osc. Str. Res.)",&
            & "Im(Osc. Str. Res.)",&
            & "|Osc. Str. Ares.|",&
            & "Re(Osc. Str. Ares.)",&
            & "Im(Osc. Str. Ares.)"
          frmt='(I8,9(1x,E23.16))'
          do lambda = 1, nexc
            if(evalre(idxsort(lambda)) > 0.0d0) then
              pm = -1.0d0
            else
              pm = 1.0d0
            end if
            write(unexc, frmt) idxsort(lambda),&
              & evalre(idxsort(lambda))*escale,&
              & (evalre(idxsort(lambda))+eshift*pm)*escale,&
              & evalim(idxsort(lambda))*escale,&
              & abs(oszstrr(idxsort(lambda), o1)),&
              & dble(oszstrr(idxsort(lambda), o1)),&
              & aimag(oszstrr(idxsort(lambda), o1)),&
              & abs(oszstra(idxsort(lambda), o1)),&
              & dble(oszstra(idxsort(lambda), o1)),&
              & aimag(oszstra(idxsort(lambda), o1))
          end do

        else if(present(oszstra)) then

          frmt='(a1,a7,8(1x,a23))'
          write(unexc, frmt) "#", "Nr.",&
            & "Re(E)",&
            & "Re(E)-+shift",&
            & "|Osc. Str. Res.|",&
            & "Re(Osc. Str. Res.)",&
            & "Im(Osc. Str. Res.)",&
            & "|Osc. Str. Ares.|",&
            & "Re(Osc. Str. Ares.)",&
            & "Im(Osc. Str. Ares.)"
          frmt='(I8,8(1x,E23.16))'
          do lambda = 1, nexc
            if(evalre(idxsort(lambda)) > 0.0d0) then
              pm = -1.0d0
            else
              pm = 1.0d0
            end if
            write(unexc, frmt) idxsort(lambda),&
              & evalre(idxsort(lambda))*escale,&
              & (evalre(idxsort(lambda))+eshift*pm)*escale,&
              & abs(oszstrr(idxsort(lambda), o1)),&
              & dble(oszstrr(idxsort(lambda), o1)),&
              & aimag(oszstrr(idxsort(lambda), o1)),&
              & abs(oszstra(idxsort(lambda), o1)),&
              & dble(oszstra(idxsort(lambda), o1)),&
              & aimag(oszstra(idxsort(lambda), o1))
          end do

        ! TDA case
        else

          frmt='(a1,a7,5(1x,a23))'
          write(unexc, frmt) "#", "Nr.",&
            & "E",&
            & "E+E_shift",&
            & "|Osc. Str. Res.|",&
            & "Re(Osc. Str. Res.)",&
            & "Im(Osc. Str. Res.)"
          frmt='(I8,5(1x,E23.16))'
          do lambda = 1, nexc
            write(unexc, frmt) lambda,&
              & evalre(lambda)*escale,&
              & (evalre(lambda)+eshift)*escale,&
              & abs(oszstrr(lambda, o1)),&
              & dble(oszstrr(lambda, o1)),&
              & aimag(oszstrr(lambda, o1))
          end do

        end if

        close(unexc)

      end do

    end subroutine writeoscillator

end module m_writeoscillator
