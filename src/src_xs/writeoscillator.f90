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

    subroutine writeoscillator(hamsize, nexc, nk, eshift, evalre, oscstrr, iqmt)
      use modxs, only: symt2, ivgigq, sptclg
      use mod_constants, only: zzero, pi
      use mod_lattice, only: omega

      ! I/O
      integer(4), intent(in) :: hamsize, nexc
      integer(4), intent(in) :: nk
      real(8), intent(in) :: eshift
      real(8), intent(in) :: evalre(hamsize)
      complex(8), intent(in) :: oscstrr(:,:)
      integer(4), intent(in), optional :: iqmt

      ! Local
      integer(4) :: o1, o2, lambda, unexc, io1, io2, io3, io4, iq, igqmt
      character(256) :: fnexc, frmt, tdastring, bsetypestring, scrtypestring
      character(256) :: syscommand, excitondir
      character(128) :: gname
      complex (8) :: buf(3,3,nexc), oscstrr_(3,3,nexc)
      real(8) :: pref
#ifdef DGRID
      character(256) :: dgrid_dotext
#endif
     
      excitondir='EXCITON'
      syscommand = 'test ! -e '//trim(adjustl(excitondir))//' && mkdir '//trim(adjustl(excitondir))
      call system(trim(adjustl(syscommand)))

      if(input%xs%bse%coupling) then
        tdastring=''
      else
        if(input%xs%bse%chibarq) then 
          tdastring="-TDA-BAR"
        else
          tdastring="-TDA"
        end if
      end if

      bsetypestring = '-'//trim(input%xs%bse%bsetype)//trim(tdastring)
      scrtypestring = '-'//trim(input%xs%screening%screentype)

      ! Loop over optical components
      if(present(iqmt)) then 
        iq = iqmt
      else
        iq = 1
      end if
      
      if(iqmt==1) then
        if (.not. input%groundstate%tevecsv) then 
          ! -1 * 2 * 4 pi * 1/V * 1/nk
          pref = -2.d0*4.d0*pi/omega/nk
        else
          ! -1 * 4 pi * 1/V * 1/nk
          pref = -4.d0*pi/omega/nk
        end if
      else
        if (.not. input%groundstate%tevecsv) then
          ! 2 * ( (4 pi/|Gmt+qmt|^2)^1/2 )^2 * 1/V * 1/nk
          igqmt = ivgigq(ivgmt(1,iqmt),ivgmt(2,iqmt),ivgmt(3,iqmt),iqmt)
          pref = 2.0d0*sptclg(igqmt,iqmt)**2/omega/nk
        else
          ! ( (4 pi/|Gmt+qmt|^2)^1/2 )^2 * 1/V * 1/nk
          igqmt = ivgigq(ivgmt(1,iqmt),ivgmt(2,iqmt),ivgmt(3,iqmt),iqmt)
          pref = sptclg(igqmt,iqmt)**2/omega/nk
        end if
      end if
      if (input%xs%bse%chibarq) then
        pref=-pref
      end if
      ! Symmetrize the oscillator strength for qmt=0
      buf(:,:,:)=zzero
      if (iq .eq. 1) then
        io1=1
        io2=3
        ! Set diagonal values
        do o1=io1,io2
          buf(o1,o1,:)=pref*abs(oscstrr(:,o1))**2
        end do
        ! Set off-diagonal values
        if (input%xs%dfoffdiag) then
          do o1=io1,io2
            do o2=io1,io2
              if (o1 .ne. o2) then
                buf(o1,o2,:)=pref*conjg(oscstrr(:,o1))*oscstrr(:,o2)
              end if
            end do !o2
          end do !o1
        end if
        ! symmetrize the oscillator strength
        do o1=io1,io2
          do o2=io1,io2
            call symt2app(o1, o2, nexc, symt2, buf, oscstrr_(o1,o2,:))
          end do !o2
        end do !o1
      else ! qmt != 0
        oscstrr_(1,1,:)=pref*abs(oscstrr(:,1))**2
      end if 
      !write hdf5 output
      gname="excitons"//trim(bsetypestring)//trim(scrtypestring)
      call write_excitons_hdf5(hamsize, nexc, eshift, evalre, oscstrr,&
      & gname, iqmt=iq)

      io1=1
      io2=3
      if(iq /= 1) then 
        io1=1
        io2=1
      end if

      do o1=io1,io2
       if ((iq .eq. 1) .and. (input%xs%dfoffdiag)) then
         io3=1
         io4=3
       else
         io3=o1
         io4=o1
       end if
       do o2=io3, io4
#ifdef DGRID
        ! Stk: Add case of double grid
        if(dgrid) then 

          write(dgrid_dotext, '("_SG", i3.3, ".OUT")') iksubpt

          call genfilname(basename='EXCITON', tq0=.true., oc1=o1, oc2=o2,&
            & bsetype=trim(bsetypestring), scrtype=trim(scrtypestring),&
            & nar= .not. input%xs%bse%aresbse, dotext=dgrid_dotext, filnam=fnexc)

        else

          call genfilname(basename='EXCITON', tq0=.true., oc1=o1, oc2=o2,&
            & bsetype=trim(bsetypestring), scrtype=trim(scrtypestring),&
            & nar= .not. input%xs%bse%aresbse, filnam=fnexc)
          
        endif
#else
        if(iq /= 1) then 
          call genfilname(basename='EXCITON', iqmt=iq,&
            & bsetype=trim(bsetypestring), scrtype=trim(scrtypestring),&
            & nar= .not. input%xs%bse%aresbse, filnam=fnexc)
        else
          call genfilname(basename='EXCITON', tq0=.true., oc1=o1, oc2=o2,&
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
            & abs(oscstrr_(o1,o2,lambda)),&
            & dble(oscstrr_(o1,o2,lambda)),&
            & aimag(oscstrr_(o1,o2,lambda))
        end do

        close(unexc)
        end do
      end do

    end subroutine writeoscillator

end module m_writeoscillator
