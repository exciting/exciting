
! Copyright (C) 2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine checkinput
  use modmain
  use modinput
  implicit none
  integer :: is,i
  character(1024) :: message
  if (input%structure%epslat.le.0.d0) then
    write(*,*)
    write(*,'("Error(checkinput): /input/structure/@epslat <= 0 : ",G18.10)') input%structure%epslat
    write(*,*)
    stop
  end if
  if (associated(input%groundstate))then
    if (input%groundstate%radkpt.le.0.d0) then
      write(*,*)
      write(*,'("Error(checkinput): /input/groundstate/@radkpt <= 0 : ",G18.10)') input%groundstate%radkpt
      write(*,*)
      stop
    end if
  end if
  if (associated(input%xs)) then
    if ((input%xs%ngridk(1).le.0).or.(input%xs%ngridk(2).le.0).or.(input%xs%ngridk(3).le.0)) then
      write(*,*)
      write(*,'("Error(checkinput): invalid /input/xs/@ngridk : ",3I8)') input%xs%ngridk
      write(*,*)
      stop
    end if
  end if
  if (associated(input%xs)) then
    if (associated(input%xs%screening)) then
      ! default: (0, 0, 0) (caught)
      if ((input%xs%screening%ngridk(1).lt.0).or.(input%xs%screening%ngridk(2).lt.0).or.(input%xs%screening%ngridk(3).lt.0)) then
        write(*,*)
        write(*,'("Error(checkinput): invalid /input/xs/screening/@ngridk : ",3I8)') input%xs%screening%ngridk
        write(*,*)
        stop
      end if
    end if
  end if
  if (associated(input%phonons)) then
    if ((input%phonons%ngridq(1).le.0).or.(input%phonons%ngridq(2).le.0).or.(input%phonons%ngridq(3).le.0)) then
      write(*,*)
      write(*,'("Error(checkinput): invalid /input/phonons/@ngridq : ",3I8)') input%phonons%ngridq
      write(*,*)
      stop
    end if
  end if
  if (associated(input%xs))then
    if ((input%xs%ngridq(1).le.0).or.(input%xs%ngridq(2).le.0).or.(input%xs%ngridq(3).le.0)) then
      write(*,*)
      write(*,'("Error(checkinput): invalid /input/xs/@ngridq : ",3I8)') input%xs%ngridq
      write(*,*)
      stop
    end if
  end if
  if (associated(input%groundstate))then
    if (input%groundstate%rgkmax.le.0.d0) then
      write(*,*)
      write(*,'("Error(checkinput): /input/groundstate/@rgkmax <= 0 : ",G18.10)') input%groundstate%rgkmax
      write(*,*)
      stop
    end if
  end if
  if (associated(input%xs))then
    if (input%xs%rgkmax.le.0.d0) then
      write(*,*)
      write(*,'("Error(checkinput): /input/xs/@rgkmax <= 0 : ",G18.10)') input%xs%rgkmax
      write(*,*)
      stop
    end if
  end if
  if (associated(input%xs))then
    if (associated(input%xs%screening))then
      ! default: 0.0 (caught)
      if (input%xs%screening%rgkmax.lt.0.d0) then
        write(*,*)
        write(*,'("Error(checkinput): /input/xs/screening/@rgkmax <= 0 : ",G18.10)') input%xs%screening%rgkmax
        write(*,*)
        stop
      end if
    end if
  end if
  if (associated(input%xs))then
    if (associated(input%xs%BSE))then
      ! default: 0.0 (caught)
      if (input%xs%BSE%rgkmax.lt.0.d0) then
        write(*,*)
        write(*,'("Error(checkinput): /input/xs/BSE/@rgkmax <= 0 : ",G18.10)') input%xs%BSE%rgkmax
        write(*,*)
        stop
      end if
    end if
  end if
  if (associated(input%groundstate))then
    if (input%groundstate%lmaxapw.lt.0) then
      write(*,*)
      write(*,'("Error(checkinput): /input/groundstate/@lmaxapw < 0 : ",I8)') input%groundstate%lmaxapw
      write(*,*)
      stop
    end if
    if (input%groundstate%lmaxapw.ge.maxlapw) then
      write(*,*)
      write(*,'("Error(checkinput): /input/groundstate/@lmaxapw too large : ",I8)') input%groundstate%lmaxapw
      write(*,'("Adjust maxlapw in modmain and recompile code")')
      write(*,*)
      stop
    end if
  end if
  if (associated(input%xs))then
    if (input%xs%lmaxapw.lt.0) then
      write(*,*)
      write(*,'("Error(checkinput): /input/xs/@lmaxapw < 0 : ",I8)') input%xs%lmaxapw
      write(*,*)
      stop
    end if
    if (input%xs%lmaxapw.ge.maxlapw) then
      write(*,*)
      write(*,'("Error(checkinput): /input/xs/@lmaxapw too large : ",I8)') input%xs%lmaxapw
      write(*,'("Adjust maxlapw in modmain and recompile code")')
      write(*,*)
      stop
    end if
  end if
  if (associated(input%groundstate))then
    if (input%groundstate%lmaxvr.lt.3) then
      write(*,*)
      write(*,'("Error(checkinput): /input/groundstate/@lmaxvr < 3 : ",I8)') input%groundstate%lmaxvr
      write(*,*)
      stop
    end if
  end if
  if (associated(input%groundstate))then
    if (input%groundstate%lmaxmat.lt.0) then
      write(*,*)
      write(*,'("Error(checkinput): /input/groundstate/@lmaxmat < 0 : ",I8)') input%groundstate%lmaxmat
      write(*,*)
      stop
    end if
  end if
  if (associated(input%xs))then
    if (input%xs%lmaxmat.lt.0) then
      write(*,*)
      write(*,'("Error(checkinput): /input/xs/@lmaxmat < 0 : ",I8)') input%xs%lmaxmat
      write(*,*)
      stop
    end if
  end if
  if (associated(input%groundstate))then
    if (input%groundstate%lmaxinr.lt.0) then
      write(*,*)
      write(*,'("Error(checkinput): /input/groundstate/@lmaxinr < 0 : ",I8)') input%groundstate%lmaxinr
      write(*,*)
      stop
    end if
  end if
  if (associated(input%groundstate))then
    if (input%groundstate%npsden.lt.2) then
      write(*,*)
      write(*,'("Error(checkinput): /input/groundstate/@npsden < 2 : ",I8)') input%groundstate%npsden
      write(*,*)
      stop
    end if
  end if
  if (associated(input%groundstate))then
    if (input%groundstate%swidth.lt.1.d-9) then
      write(*,*)
      write(*,'("Error(checkinput): /input/groundstate/@swidth too small or negative : ",G18.10)') &
       input%groundstate%swidth
      write(*,*)
      stop
    end if
  end if
  if (associated(input%xs))then
    if (input%xs%swidth.lt.1.d-9) then
      write(*,*)
      write(*,'("Error(checkinput): /input/xs/@swidth too small or negative : ",G18.10)') &
       input%xs%swidth
      write(*,*)
      stop
    end if
  end if
  if (associated(input%groundstate))then
    if (input%groundstate%epsocc.le.0.d0) then
      write(*,*)
      write(*,'("Error(checkinput): /input/groundstate/@epsocc <= 0 : ",G18.10)') input%groundstate%epsocc
      write(*,*)
      stop
    end if
  end if
  if (associated(input%groundstate))then
    if (input%groundstate%epschg.le.0.d0) then
      write(*,*)
      write(*,'("Error(checkinput): /input/groundstate/@epschg <= 0 : ",G18.10)') input%groundstate%epschg
      write(*,*)
      stop
    end if
  end if
  if (associated(input%groundstate))then
    if (input%groundstate%nempty.le.0) then
      write(*,*)
      write(*,'("Error(checkinput): /input/groundstate/@nempty <= 0 : ",I8)') input%groundstate%nempty
      write(*,*)
      stop
    end if
  end if
  if (associated(input%xs))then
    if (input%xs%nempty.le.0) then
      write(*,*)
      write(*,'("Error(checkinput): /input/xs/@nempty <= 0 : ",I8)') input%xs%nempty
      write(*,*)
      stop
    end if
  end if
  if (associated(input%xs)) then
    if (associated(input%xs%screening)) then
      ! default: 0 (caught)
      if (input%xs%screening%nempty.lt.0) then
        write(*,*)
        write(*,'("Error(checkinput): /input/xs/screening/@nempty <= 0 : ",I8)') input%xs%screening%nempty
        write(*,*)
        stop
      end if
    end if
  end if
  if (associated(input%groundstate))then
    if (input%groundstate%beta0.lt.0.d0) then
      write(*,*)
      write(*,'("Error(checkinput): /input/groundstate/@beta0 < 0 : ",G18.10)') input%groundstate%beta0
      write(*,*)
      stop
    end if
  end if
  if (associated(input%groundstate))then
    if (input%groundstate%betainc.lt.1.d0) then
      write(*,*)
      write(*,'("Error(checkinput): /input/groundstate/@betainc < 1 : ",G18.10)') input%groundstate%betainc
      write(*,*)
      stop
    end if
  end if
  if (associated(input%groundstate))then
    if ((input%groundstate%betadec.le.0.d0).or.(input%groundstate%betadec.gt.1.d0)) then
      write(*,*)
      write(*,'("Error(checkinput): /input/groundstate/@betadec should be in (0,1] : ",G18.10)') &
       input%groundstate%betadec
      write(*,*)
      stop
    end if
  end if
  if (associated(input%groundstate))then
    if (input%groundstate%maxscl.lt.0) then
      write(*,*)
      write(*,'("Error(checkinput): /input/groundstate/@maxscl < 0 : ",I8)') input%groundstate%maxscl
      write(*,*)
      stop
    end if
  end if
  if (associated(input%groundstate))then
    if (input%groundstate%cfdamp.lt.0.d0) then
      write(*,*)
      write(*,'("Error(checkinput): /input/groundstate/@cfdamp < 0 : ",G18.10)') input%groundstate%cfdamp
      write(*,*)
      stop
    end if
  end if
  if (associated(input%structure%speciesarray)) then
    if (size(input%structure%speciesarray).gt.maxspecies) then ! (nspecies > maxspecies)
      write(*,*)
      write(*,'("Error(checkinput): number of species too large : ",I8)') size(input%structure%speciesarray)
      write(*,'("Adjust maxspecies in modmain and recompile code")')
      write(*,*)
      stop
    end if
    do is=1,nspecies
      if (size(input%structure%speciesarray(is)%species%atomarray).gt.maxatoms) then ! (natoms(is) > maxatoms)
        write(*,*)
        write(*,'("Error(checkinput): number of atoms too large : ",I8)') size(input%structure%speciesarray(is)%species%atomarray)
        write(*,'(" for species ",I4)') is
        write(*,'("Adjust maxatoms in modmain and recompile code")')
        write(*,*)
        stop
      end if
    end do
  end if
  if (associated(input%properties)) then
    if (associated(input%properties%dos)) then
      if (input%properties%dos%nwdos.lt.2) then
       write(*,*)
       write(*,'("Error(checkinput): /input/properties/dos/@nwdos < 2 : ",I8)') input%properties%dos%nwdos
       write(*,*)
       stop
      end if
      if (input%properties%dos%ngrdos.lt.1) then
       write(*,*)
       write(*,'("Error(checkinput): /input/properties/dos/@ngrdos < 1 : ",I8)') input%properties%dos%ngrdos
       write(*,*)
       stop
      end if
      if (input%properties%dos%nsmdos.lt.0) then
       write(*,*)
       write(*,'("Error(checkinput): /input/properties/dos/@nsmdos < 0 : ",I8)') input%properties%dos%nsmdos
       write(*,*)
       stop
      end if
      if (input%properties%dos%winddos(1).ge.input%properties%dos%winddos(2)) then
       write(*,*)
       write(*,'("Error(checkinput): /input/properties/dos/@winddos(1) >= /input/properties/dos/@winddos(2) : ",&
         &2G18.10)') input%properties%dos%winddos
       write(*,*)
       stop
      end if
    end if
  end if
  if (associated(input%properties)) then
    if (associated(input%properties%fermisurfaceplot)) then
      if (input%properties%fermisurfaceplot%nstfsp.le.0) then
        write(*,*)
        write(*,'("Error(checkinput): /input/properties/fermisurfaceplot/@nstfsp <= 0 : ",I8)') &
          input%properties%fermisurfaceplot%nstfsp
        write(*,*)
        stop
      end if
    end if
  end if
  if (associated(input%groundstate)) then
    if (input%groundstate%lradstep.le.0) then
      write(*,*)
      write(*,'("Error(checkinput): /input/groundstate/@lradstep <= 0 : ",I8)') input%groundstate%lradstep
      write(*,*)
      stop
    end if
  end if
  if (associated(input%groundstate)) then
    if (input%groundstate%nprad.lt.2) then
      write(*,*)
      write(*,'("Error(checkinput): /input/groundstate/@nprad < 2 : ",I8)') input%groundstate%nprad
      write(*,*)
      stop
    end if
  end if
  if (associated(input%groundstate)) then
    if (associated(input%groundstate%solver)) then
      if (input%groundstate%solver%evaltol.le.0.d0) then
        write(*,*)
        write(*,'("Error(checkinput): /input/groundstate/solver/@evaltol <= 0 : ",G18.10)') input%groundstate%solver%evaltol
        write(*,*)
        stop
      end if
    end if
  end if
  if (associated(input%groundstate)) then
    if (input%groundstate%deband.lt.0.d0) then
      write(*,*)
      write(*,'("Error(checkinput): /input/groundstate/@deband < 0 : ",G18.10)') input%groundstate%deband
      write(*,*)
      stop
    end if
  end if
  if (associated(input%groundstate)) then
    if (associated(input%groundstate%spin)) then
      if (input%groundstate%spin%taufsm.lt.0.d0) then
        write(*,*)
        write(*,'("Error(checkinput): /input/groundstate/spin/@taufsm < 0 : ",G18.10)') input%groundstate%spin%taufsm
        write(*,*)
        stop
      end if
    end if
  end if
  if (associated(input%groundstate)) then
    if (input%structure%rmtapm(1).lt.0.d0) then
      write(*,*)
      write(*,'("Error(checkinput): /input/groundstate/@rmtapm(1) < 0 : ",G18.10)') input%structure%rmtapm(1)
      write(*,*)
      stop
    end if
    if ((input%structure%rmtapm(2).le.0.d0).or.(input%structure%rmtapm(2).gt.1.d0)) then
      write(*,*)
      write(*,'("Error(checkinput): /input/groundstate/@rmtapm(2) not in (0,1] : ",G18.10)') input%structure%rmtapm(2)
      write(*,*)
      stop
    end if
  end if
  if (associated(input%groundstate)) then
    if (associated(input%groundstate%OEP)) then
      if (input%groundstate%OEP%maxitoep.lt.1) then
        write(*,*)
        write(*,'("Error(checkinput): /input/groundstate/OEP/@maxitoep < 1 : ",I8)') input%groundstate%OEP%maxitoep
        write(*,*)
        stop
      end if
    end if
  end if
  if (associated(input%groundstate)) then
    if (associated(input%groundstate%OEP)) then
      if ((input%groundstate%OEP%tauoep(1).lt.0.d0).or. &
       (input%groundstate%OEP%tauoep(2).lt.0.d0).or. &
       (input%groundstate%OEP%tauoep(3).lt.0.d0)) then
        write(*,*)
        write(*,'("Error(checkinput): /input/groundstate/OEP/@tauoep < 0 : ",3G18.10)') input%groundstate%OEP%tauoep
        write(*,*)
        stop
      end if
    end if
  end if
  if (associated(input%properties)) then
    if (associated(input%properties%masstensor)) then
      if ((input%properties%masstensor%ndspem.lt.1).or.(input%properties%masstensor%ndspem.gt.3)) then
        write(*,*)
        write(*,'("Error(checkinput): /input/properties/masstensor/@ndspem out of range : ",I8)') &
          input%properties%masstensor%ndspem
        write(*,*)
        stop
      end if
    end if
  end if
  if (associated(input%structure%speciesarray)) then
    do is=1,size(input%structure%speciesarray)
      if (associated(input%structure%speciesarray(is)%species%LDAplusU)) then
        if (input%structure%speciesarray(is)%species%LDAplusU%l.gt.lmaxlu) then
          write(*,*)
          write(*,'("Error(checkinput): /input/structure/species/LDAplusU/@l > lmaxlu in lda+u block : ",2I8)') &
           input%structure%speciesarray(is)%species%LDAplusU%l,lmaxlu
          write(*,*)
          stop
        end if
      end if
    end do
  end if
  if (associated(input%groundstate)) then
    if (associated(input%groundstate%RDMFT)) then
      if (input%groundstate%RDMFT%rdmmaxscl.lt.0) then
        write(*,*)
        write(*,'("Error(checkinput): /input/groundstate/RDMFT/@rdmmaxscl < 0 : ",I8)') &
         input%groundstate%RDMFT%rdmmaxscl
        write(*,*)
      end if
    end if
  end if
  if (associated(input%groundstate)) then
    if (associated(input%groundstate%RDMFT)) then
      if (input%groundstate%RDMFT%taurdmn.lt.0.d0) then
        write(*,*)
        write(*,'("Error(checkinput): /input/groundstate/RDMFT/@taurdmn < 0 : ",G18.10)') &
         input%groundstate%RDMFT%taurdmn
        write(*,*)
        stop
      end if
    end if
  end if
  if (associated(input%groundstate)) then
    if (associated(input%groundstate%RDMFT)) then
      if (input%groundstate%RDMFT%taurdmc.lt.0.d0) then
        write(*,*)
        write(*,'("Error(checkinput): /input/groundstate/RDMFT/@taurdmc < 0 : ",G18.10)') &
         input%groundstate%RDMFT%taurdmc
        write(*,*)
        stop
      end if
    end if
  end if
  if (associated(input%groundstate)) then
    if (associated(input%groundstate%RDMFT)) then
      if ((input%groundstate%RDMFT%rdmalpha.le.0.d0).or.(input%groundstate%RDMFT%rdmalpha.ge.1.d0)) then
        write(*,*)
        write(*,'("Error(checkinput): /input/groundstate/RDMFT/@rdmalpha not in (0,1) : ",G18.10)') &
         input%groundstate%RDMFT%rdmalpha
        write(*,*)
        stop
      end if
    end if
  end if
  if (associated(input%groundstate)) then
      if (associated(input%groundstate%RDMFT)) then
      if (input%groundstate%RDMFT%rdmtemp.lt.0.d0) then
        write(*,*)
        write(*,'("Error(checkinput): /input/groundstate/RDMFT/@rdmtemp < 0 : ",G18.10)') &
         input%groundstate%RDMFT%rdmtemp
        write(*,*)
        stop
      end if
    end if
  end if
  if (associated(input%groundstate)) then
    if (associated(input%groundstate%spin)) then
      if ((input%groundstate%spin%reducebf.lt.0.d0).or.(input%groundstate%spin%reducebf.gt.1.d0)) then
        write(*,*)
        write(*,'("Error(checkinput): /input/groundstate/spin/@reducebf not in [0,1] : ",G18.10)') input%groundstate%spin%reducebf
        write(*,*)
        stop
      end if
    end if
  end if
! ngridk optional
  if (associated(input%groundstate)) then
    if(any(input%groundstate%ngridk .le. 0) .and. (.not.input%groundstate%autokpt) .and. (input%groundstate%nktot .eq. 0)) then
      write(*,*)
      write(*,'("Error(checkinput): components in /input/groundstate/@ngridk < 1 : ",3I10)') input%groundstate%ngridk
      write(*,'("specifiy either /input/groundstate/@ngridk or /input/groundstate/@autokpt or /input/groundstate/@nktot ")')
      write(*,*)
      stop
    end if
  end if

! excited states elements and attributes not already contained in remaining part
  if (associated(input%xs))then
    ! default: -1 (caught)
    if (input%xs%lmaxapwwf.lt.-1) then
      write(*,*)
      write(*,'("Error(checkinput): /input/xs/@lmaxapwwf < 0 : ",I8)') input%xs%lmaxapwwf
      write(*,*)
      stop
    end if
  end if
  if (associated(input%xs)) then
    if (associated(input%xs%tddft)) then
      if ((input%xs%tddft%mdfqtype.lt.0).or.(input%xs%tddft%mdfqtype.gt.1)) then
        write(*,*)
        write(*,'("Error(checkinput): /input/xs/tddft/@mdfqtype not in {0,1} : ",I8)') input%xs%tddft%mdfqtype
        write(*,*)
        stop
      end if
      if (input%xs%tddft%mdfqtype.eq.1) then
        write(*,*)
        write(*,'("Error(checkinput): /input/xs/tddft/@mdfqtype=1; not compatible with xcifc-&
         &routine, if local fields are neglected, needs to project out special &
         &G-vector - code limitation")')
        write(*,'(" use equivalent choice mdfqtype=0 in place")')
        write(*,*)
        stop
      end if
    end if
  end if
  if (associated(input%xs)) then
    if (input%xs%lmaxemat.lt.0) then
      write(*,*)
      write(*,'("Error(checkinput): /input/xs/@lmaxemat < 0 : ",I8)') input%xs%lmaxemat
      write(*,*)
      stop
    end if
  end if
  if (associated(input%xs)) then
    if (associated(input%xs%tddft)) then
      if (input%xs%tddft%lmaxalda.lt.0) then
        write(*,*)
        write(*,'("Error(checkinput): /input/xs/tddft/@lmaxalda < 0 : ",I8)') input%xs%tddft%lmaxalda
        write(*,*)
        stop
      end if
    end if
  end if
  if (associated(input%xs)) then
    if (associated(input%xs%BSE)) then
      if (input%xs%BSE%nleblaik.lt.0) then
        write(*,*)
        write(*,'("Error(checkinput): /input/xs/BSE/@nleblaik < 0 : ",I8)') input%xs%BSE%nleblaik
        write(*,*)
        stop
      end if
    end if
  end if
  if (associated(input%xs)) then
    if (associated(input%xs%BSE)) then
      if (input%xs%BSE%lmaxdielt.lt.0) then
        write(*,*)
        write(*,'("Error(checkinput): /input/xs/BSE/@lmaxdielt < 0 : ",I8)') input%xs%BSE%lmaxdielt
        write(*,*)
        stop
      end if
    end if
  end if
  if (associated(input%xs)) then
    if (associated(input%xs%tddft)) then
      if (input%xs%tddft%nwacont.lt.0) then
        write(*,*)
        write(*,'("Error(checkinput): /input/xs/tddft/@nwacont <= 0 : ",g18.10)') input%xs%tddft%nwacont
        write(*,*)
        stop
      end if
    end if
  end if
  if (associated(input%xs)) then
    if (input%xs%broad.le.0.d0) then
      write(*,*)
      write(*,'("Warning(checkinput): /input/xs/@broad <= 0 : ",g18.10)') input%xs%broad
      write(*,*)
    end if
  end if
  if (associated(input%xs)) then
    if (input%xs%epsdfde.le.0.d0) then
      write(*,*)
      write(*,'("Warning(checkinput): /input/xs/@epsdfde <= 0 : ",g18.10)') input%xs%epsdfde
      write(*,*)
    end if
  end if
  if (associated(input%xs)) then
    if (associated(input%xs%tddft)) then
      if (input%xs%tddft%fxcbsesplit.le.0) then
        write(*,*)
        write(*,'("Error(checkinput): /input/xs/tddft/@fxcbsesplit <= 0 : ",g18.10)') input%xs%tddft%fxcbsesplit
        write(*,*)
        stop
      end if
    end if
  end if
  if (associated(input%xs)) then
    if (associated(input%xs%BSE)) then
      if ((input%xs%BSE%nstlbsemat(1).lt.0).or.(input%xs%BSE%nstlbsemat(2).lt.0).or. &
        (input%xs%BSE%nstlbsemat(3).lt.0).or.(input%xs%BSE%nstlbsemat(4).lt.0)) then
        write(*,*)
        write(*,'("Error(checkinput): /input/xs/BSE/@nstlbsemat(1), /input/xs/BSE/@nstlbsemat(2), ")')
        write(*,'(" /input/xs/BSE/@nstlbsemat(3) or /input/xs/BSE/@nstlbsemat(4) <= 0 : ",4I8)') &
          input%xs%BSE%nstlbsemat
        write(*,*)
        stop
      end if
    end if
  end if
  if (associated(input%xs)) then
    if (associated(input%xs%BSE)) then
      if ((input%xs%BSE%nstlbse(1).lt.0).or.(input%xs%BSE%nstlbse(2).lt.0).or. &
        (input%xs%BSE%nstlbse(3).lt.0).or.(input%xs%BSE%nstlbse(4).lt.0)) then
        write(*,*)
        write(*,'("Error(checkinput): /input/xs/BSE/@nstlbse(1), /input/xs/BSE/@nstlbse(2), ")')
        write(*,'(" /input/xs/BSE/@nstlbse(3) or /input/xs/BSE/@nstlbse(4) <= 0 : ",4I8)') input%xs%BSE%nstlbse
        write(*,*)
        stop
      end if
    end if
  end if
end subroutine
