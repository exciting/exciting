
! Copyright (C) 2007-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_genfilname
  implicit none
contains

!BOP
! !ROUTINE: genfilname
! !INTERFACE:
  subroutine genfilname(nodotpar,basename,etype,asc,bzsampl,acont,&
       nar,tord,nlf,fxctype,scrtype,bsetype,markfxcbse,tq0,oc1,oc2,iq,iqmt,procs,rank,dotext, &
       setfilext,revertfilext,appfilext,filnam,fileext)
! !USES:
    use modmain, only: filext
    use modxs, only: filextrevert
! !DESCRIPTION:
!   Generates file name and extension according to optional input parameters (
!   see routine).
!   Interpret bzsampl variable as default (Lorentzian) for 0, as
!   tetrahedron method for 1. Trilinear method to be followed.
!
! !REVISION HISTORY:
!   Created October 2007 (Sagmeister)
!EOP
!BOC
    implicit none
    ! arguments
    integer, optional, intent(in) :: bzsampl,fxctype,oc1,oc2,iq,iqmt,procs,rank
    integer, optional, intent(in) :: etype
    logical, optional, intent(in) :: nodotpar,asc,acont,nar,tord,nlf,tq0,markfxcbse
    logical, optional, intent(in) :: revertfilext,setfilext,appfilext
    character(*), optional, intent(in) :: basename,dotext,scrtype,bsetype
    character(256), optional, intent(out) :: filnam,fileext
    ! local variables
    logical :: nodot0,revert,setfxt,appfxt,dotxt,oct,lnar
    character(*), parameter :: thisnam = 'genfilname'
    character(256) :: s,s1
    ! if file extension in "modmain" is to be reset to last value: reset
    ! else store current file extension
    revert=.false.
    setfxt=.false.
    if (present(revertfilext)) revert=revertfilext
    if (present(setfilext)) setfxt=setfilext
    if (revert) then
       filext=trim(filextrevert)
    else if (setfxt) then
       filextrevert=filext
    end if
    appfxt=.false.
    if (present(appfilext)) appfxt=appfilext
    dotxt=.false.
    if (present(dotext)) dotxt=.true.
    if ((appfxt.and.setfxt).or.(appfxt.and.dotxt)) then
       write(*,'(a)') 'Error('//trim(thisnam)//'): specified appfxt together &
            &with setfilext or dotext'
       call terminate
    end if
    oct=present(oc1).and.present(oc2)
    ! dot in front of filename in parallel output for rank eq. zero
    nodot0=.false.
    if (present(nodotpar)) nodot0=nodotpar
    ! start with empty string
    s=''
    ! type of band combinations for plane wave matrix elements
    if (present(etype)) then
       select case(etype)
       case(0)
          ! all band combinations
          s=trim(s)//'_FULL'
       case(1)
          ! v-c anc c-v combinations for response function
       case(2)
          ! v-v and c-c combinations for screened interaction
          s=trim(s)//'_SCRI'
       case default
          write(*,'(a)') 'Error('//trim(thisnam)//'): unknown etype: ', &
               etype
          call terminate
       end select
    end if
    ! ascii output identifier
    if (present(asc)) then
       if (asc) then
          s=trim(s)//'_ASC'
       end if
    end if
    ! sampling of Brillouine zone
    if (present(bzsampl)) then
       select case(bzsampl)
       case(0)
          ! do nothing (Lorentzian broadening)
       case(1)
          ! tetrahedron method
          s=trim(s)//'_TET'
       case default
          write(*,'(a)') 'Error('//trim(thisnam)//'): unknown bzsampl: ', &
               bzsampl
          call terminate
       end select
    end if
    ! analytic continuation
    if (present(acont)) then
       if (acont) then
          s=trim(s)//'_AC'
       end if
    end if
    ! exclusion of anti-resonant part
    lnar=.false.
    if (present(nar)) then
       if (nar) then
          lnar=.true.
          s=trim(s)//'_NAR'
       end if
    end if
    ! time-ordering
    if (present(tord)) then
       if (tord.and.(.not.lnar)) then
          s=trim(s)//'_TORD'
       end if
    end if
    ! no local field effects
    if (present(nlf)) then
       if (nlf) then
          s=trim(s)//'_NLF'
       end if
    end if
    ! xc-kernel type
    if (present(fxctype)) then
       write(s1,'("_FXC",i2.2)') fxctype
       s=trim(s)//trim(s1)
    end if
    ! BSE effective Hamiltonian type
    if (present(bsetype)) then
       write(s1,'("_BSE",a)') trim(adjustl(bsetype))
       s=trim(s)//trim(s1)
    end if
    ! screening type in screened Coulomb interaction
    if (present(scrtype)) then
       write(s1,'("_SCR",a)') trim(adjustl(scrtype))
       s=trim(s)//trim(s1)
    end if
    ! optical components
    if (present(tq0).and.oct) then
       if (tq0) then
          write(s1,'("_OC",2i1.1)') oc1,oc2
          s=trim(s)//trim(s1)
       end if
    end if
    ! mark file if it is generated in combination with fxc-BSE
    if (present(markfxcbse)) then
       if (markfxcbse) then
         s=trim(s)//"_FXCBSE"
       end if
    end if
    ! Q-point (finite momentum transfer)
    if (present(iqmt)) then
       write(s1,'("_QMT",i3.3)') iqmt
       s=trim(s)//trim(s1)
    end if
    ! q-point
    if (present(iq)) then
       write(s1,'("_Q",i5.5)') iq
       s=trim(s)//trim(s1)
    end if
    ! parallelization
    if (present(rank).and.present(procs)) then
       if ((procs > 1).and.((nodot0.and.(rank > 0)).or.(.not.nodot0))) then
          ! tag for rank
          write(s1,'("_par",i3.3)') rank+1
          s=trim(s)//trim(s1)
       end if
    end if
    ! extension (including the dot)
    if (present(dotext)) then
       s=trim(s)//trim(dotext)
    else if (appfxt) then
       s=trim(s)//trim(filext)
    else
       s=trim(s)//'.OUT'
    end if
    ! assign file extension if required
    if (present(fileext)) fileext=trim(s)
    ! assign global file extension if required
    if (setfxt) filext=trim(s)
    ! basename
    if (present(basename)) s=trim(basename)//trim(s)
    ! dot in front of filename determined by procs, rank and nodotpar
    if (present(rank).and.present(procs)) then
       if (((procs > 1).and.(rank > 0)).or. &
            ((.not.nodot0).and.(procs > 1).and.(rank == 0))) then
          s='.'//trim(s)
       end if
    end if
    if (present(filnam)) filnam=trim(s)
  end subroutine genfilname
!EOC

end module m_genfilname
