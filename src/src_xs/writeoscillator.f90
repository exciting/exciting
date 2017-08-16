module m_writeoscillator
  use modmpi
  use modinput, only: input
  use modxs, only: bsed, escale, ivgmt, vqlmt, vgcmt, vqcmt
#ifdef DGRID
  use modxs, only: dgrid, iksubpt
#endif
  use m_genfilname
  use m_getunit
  use m_write_hdf5

  implicit none

  contains

    subroutine writeoscillator(hamsize, nexc, eshift, evalre, oscstrr, iqmt)
      ! I/O
      integer(4), intent(in) :: hamsize, nexc
      real(8), intent(in) :: eshift
      real(8), intent(in) :: evalre(hamsize)
      complex(8), intent(in) :: oscstrr(:,:)
      integer(4), intent(in), optional :: iqmt

      ! Local
      integer(4) :: o1, lambda, unexc, io1, io2, iq
      character(256) :: fnexc, frmt, tdastring, bsetypestring, scrtypestring
      character(256) :: syscommand, excitondir
      character(128) :: gname
#ifdef DGRID
      character(256) :: dgrid_dotext
#endif
     
      excitondir='EXCITON'
      syscommand = 'test ! -e '//trim(adjustl(excitondir))//' && mkdir '//trim(adjustl(excitondir))
      call system(trim(adjustl(syscommand)))

      if(input%xs%bse%coupling) then
        tdastring=''
      else
        tdastring="-TDA"
      end if

      bsetypestring = '-'//trim(input%xs%bse%bsetype)//trim(tdastring)
      scrtypestring = '-'//trim(input%xs%screening%screentype)

      !write hdf5 output
      gname="excitons"//trim(bsetypestring)//trim(scrtypestring)
      call write_excitons_hdf5(hamsize, nexc, eshift, evalre, oscstrr,&
      & gname, iqmt=iqmt)

      ! Loop over optical components
      if(present(iqmt)) then 
        iq = iqmt
      else
        iq = 1
      end if

      io1=1
      io2=3
      if(iq /= 1) then 
        io1=1
        io2=1
      end if

      do o1=io1,io2

#ifdef DGRID
        ! Stk: Add case of double grid
        if(dgrid) then 

          write(dgrid_dotext, '("_SG", i3.3, ".OUT")') iksubpt

          call genfilname(basename='EXCITON', tq0=.true., oc1=o1, oc2=o1,&
            & bsetype=trim(bsetypestring), scrtype=trim(scrtypestring),&
            & nar= .not. input%xs%bse%aresbse, dotext=dgrid_dotext, filnam=fnexc)

        else

          call genfilname(basename='EXCITON', tq0=.true., oc1=o1, oc2=o1,&
            & bsetype=trim(bsetypestring), scrtype=trim(scrtypestring),&
            & nar= .not. input%xs%bse%aresbse, filnam=fnexc)
          
        endif
#else
        if(iq /= 1) then 
          call genfilname(basename='EXCITON', iqmt=iq,&
            & bsetype=trim(bsetypestring), scrtype=trim(scrtypestring),&
            & nar= .not. input%xs%bse%aresbse, filnam=fnexc)
        else
          call genfilname(basename='EXCITON', tq0=.true., oc1=o1, oc2=o1,&
            & bsetype=trim(bsetypestring), scrtype=trim(scrtypestring),&
            & nar= .not. input%xs%bse%aresbse, filnam=fnexc)
        end if
#endif
        fnexc=trim(excitondir)//'/'//trim(fnexc)

        ! Write out exciton energies and oscillator strengths
        call getunit(unexc)
        open(unexc, file=fnexc, form='formatted', action='write', status='replace')
        write(unexc, '("#",1x,"Excitonic eigen energies and oscillator strengths")')
        write(unexc, '("#")')
        write(unexc, '("# Momentum transfer Q=G+q in lattice cooridnates")')
        write(unexc, '("# G:",3i4)') ivgmt(1:3,iq) 
        write(unexc, '("# q:",3f12.7)') vqlmt(1:3,iq) 
        write(unexc, '("# Momentum transfer Q=G+q in Cartesian cooridnates")')
        write(unexc, '("# G:",3f12.7)') vgcmt(1:3,iq) 
        write(unexc, '("# q:",3f12.7)') vqcmt(1:3,iq) 
        write(unexc, '("# Norm2(G+q)",f12.7)') norm2(vgcmt(:,iq)+vqcmt(:,iq))
        write(unexc, '("#")')
        write(unexc, '("# Energy scale",f12.7)') escale
        write(unexc, '("# E_shift : ", SP, E23.16)') eshift * escale
        write(unexc, '("#")')


        frmt='(a1,a7,5(1x,a23))'
        write(unexc, frmt) "#", "Nr.",&
          & "E",&
          & "E+E_shift",&
          & "|Osc. Str.|",&
          & "Re(Osc. Str. Res.)",&
          & "Im(Osc. Str. Res.)"
        frmt='(I8,5(1x,E23.16))'
        do lambda = 1, nexc
          write(unexc, frmt) lambda,&
            & evalre(lambda)*escale,&
            & (evalre(lambda)+eshift)*escale,&
            & abs(oscstrr(lambda, o1)),&
            & dble(oscstrr(lambda, o1)),&
            & aimag(oscstrr(lambda, o1))
        end do

        close(unexc)

      end do

    end subroutine writeoscillator

end module m_writeoscillator
