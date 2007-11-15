
! Copyright (C) 2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_genfilname
  implicit none
contains

  subroutine genfilname(nodotpar,basename,asc,bzsampl,acont,&
       nar,nlf,fxctype,tq0,oc,iq,procs,rank,dotext,setfilext,revertfilext,&
       filnam,fileext)
    use modmain, only: filext
    use modxs, only: filextrevert
    ! Generate file name and extension accoring to purpose and optional
    ! input parameters.
    ! Interpret bzsampl variable as default (Lorentzian) for 0, as
    ! Tetrahedron method for 1. Trilinear method to be followed.
    !
    ! Oktober 2007
    implicit none
    ! arguments
    integer, optional, intent(in) :: bzsampl,fxctype,oc,iq,procs,rank
    logical, optional, intent(in) :: nodotpar,asc,acont,nar,nlf,tq0
    logical, optional, intent(in) :: revertfilext,setfilext
    character(*), optional, intent(in) :: basename,dotext
    character(256), optional, intent(out) :: filnam,fileext
    ! local variables
    logical :: nodot0,revert,setfxt
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
    
    ! dot in front of filename in parallel output for rank eq. zero
    nodot0=.false.
    if (present(nodotpar)) nodot0=nodotpar
    s=''
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
    if (present(nar)) then
       if (nar) then
          s=trim(s)//'_NAR'
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
    ! optical components
    if (present(tq0).and.present(oc)) then
       if (tq0) then
          write(s1,'("_OC",i2.2)') oc*11
          s=trim(s)//trim(s1)
       end if
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

end module m_genfilname
