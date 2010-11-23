! Copyright (C) 2002-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine readspeciesxml
      Use modmain
      Use modinput, Only: input
      Use modsp
      Use FoX_dom
      Use modspdb
      Implicit None
! local variables
      Integer :: is, ist
      Integer :: io, nlx, ilx, lx, ilo,slash,errorcode
      Type (Node), Pointer :: speciesnp, speciesdbnp
      character(2048)::command
      character(256)::spfile_string
      mine0=0.d0
      if (allocated(speziesdeflist)) deallocate(speziesdeflist)
      Allocate (speziesdeflist(nspecies))
      config => newDOMConfig ()
      parseerror = .False.
! parse xml and create derived type for species definitions  speziesdeflist
      Do is = 1, nspecies
      spfile_string=""

         if (index(input%structure%speciespath,"http://").ge.1) then
         write(*,*) "index ",index(input%structure%speciespath,"http://")
             If(trim(input%structure%speciesarray(is)%species%href).eq."")then
                 write(input%structure%speciesarray(is)%species%href,*) trim(input%structure%speciespath),&
                 &trim(input%structure%speciesarray(is)%species%speciesfile)
             endif
         endif

         If(trim(input%structure%speciesarray(is)%species%href).ne."") then
             slash=index(input%structure%speciesarray(is)%species%href,"/",.true.)
             write(spfile_string,*)trim(input%structure%speciesarray(is)%species%href(slash+1:))
#ifdef CURL
             write(command,*)"curl ",trim(input%structure%speciesarray(is)%species%href)," > ",trim(spfile_string)
#endif
#ifndef CURL
             write(command,*)"wget -c ",trim(input%structure%speciesarray(is)%species%href)
#endif
!
             call system(command)
         else
             Write (spfile_string,*)trim (input%structure%speciespath) // "/" &
            & // trim (input%structure%speciesarray(is)%species%speciesfile)
	     endif
         doc => parseFile (ADJUSTL(trim(spfile_string)),config,iostat=errorcode)
          if(errorcode.ne.0) then
        	 write(*,*) "### Could not open ", ADJUSTL(trim(spfile_string))
      		 write(*,*) "### Check if file exists and if it is well-formed XML."

       	     stop
         endif
         speciesdbnp => getDocumentElement (doc)
         speciesnp => item (getElementsByTagname(speciesdbnp, "sp"), 0)
         parseerror = .False.
         speziesdeflist(is)%sp => getstructsp (speciesnp)
         Call destroy (doc)
      End Do
!
      Do is = 1, nspecies
         spsymb (is) = speziesdeflist(is)%sp%chemicalSymbol
         input%structure%speciesarray(is)%species%chemicalSymbol = &
        & spsymb (is)
         spname (is) = speziesdeflist(is)%sp%name
         spzn (is) = speziesdeflist(is)%sp%z
         spmass (is) = speziesdeflist(is)%sp%mass
         sprmin (is) = speziesdeflist(is)%sp%muffinTin%rmin
         rmt (is) = speziesdeflist(is)%sp%muffinTin%radius
         If (input%structure%speciesarray(is)%species%rmt .Gt. 0) Then
            rmt (is) = input%structure%speciesarray(is)%species%rmt
         End If
         sprmax (is) = speziesdeflist(is)%sp%muffinTin%rinf
         nrmt (is) = speziesdeflist(is)%sp%muffinTin%radialmeshPoints
         If (sprmin(is) .Le. 0.d0) Then
            Write (*,*)
            Write (*, '("Error(readinput): sprmin <= 0 : ", G18.10)') &
           & sprmin (is)
            Write (*, '(" for species ", I4)') is
            Write (*,*)
            Stop
         End If
         If (rmt(is) .Le. sprmin(is)) Then
            Write (*,*)
            Write (*, '("Error(readinput): rmt <= sprmin : ", 2G18.10)') rmt (is), sprmin (is)
            Write (*, '(" for species ", I4)') is
            Write (*,*)
            Stop
         End If
         If (sprmax(is) .Lt. rmt(is)) Then
            Write (*,*)
            Write (*, '("Error(readinput): sprmax < rmt : ", 2G18.10)') &
           & sprmax (is), rmt (is)
            Write (*,*)
            Stop
         End If
         If (nrmt(is) .Lt. 20) Then
            Write (*,*)
            Write (*, '("Error(readinput): nrmt too small : ", I8)') &
           & nrmt (is)
            Write (*, '(" for species ", I4)') is
            Write (*,*)
            Stop
         End If
         spnst (is) = size (speziesdeflist(is)%sp%AtomicStatearray)
         If (spnst(is) .Le. 0) Then
            Write (*,*)
            Write (*, '("Error(readinput): invalid spnst : ", I8)') &
           & spnst (is)
            Write (*, '(" for species ", I4)') is
            Write (*,*)
            Stop
         End If
         If (spnst(is) .Gt. maxspst) Then
            Write (*,*)
            Write (*, '("Error(readinput): too many states for species &
           &", I8)') is
            Write (*,*)
            Stop
         End If
         Do ist = 1, spnst (is)
!
            spn (ist, is) = &
           & speziesdeflist(is)%sp%AtomicStatearray(ist)%atomicState%n
            spl (ist, is) = &
           & speziesdeflist(is)%sp%AtomicStatearray(ist)%atomicState%l
            spk (ist, is) = speziesdeflist(is)%sp%AtomicStatearray(ist)%atomicState%kappa
            spocc (ist, is) = speziesdeflist(is)%sp%AtomicStatearray(ist)%atomicState%occ
            spcore (ist, is) = speziesdeflist(is)%sp%AtomicStatearray(ist)%atomicState%core
            If (spn(ist, is) .Lt. 1) Then
               Write (*,*)
               Write (*, '("Error(readinput): spn < 1 : ", I8)') spn &
              & (ist, is)
               Write (*, '(" for species ", I4)') is
               Write (*, '(" and state ", I4)') ist
               Write (*,*)
               Stop
            End If
            If (spl(ist, is) .Lt. 0) Then
               Write (*,*)
               Write (*, '("Error(readinput): spl < 0 : ", I8)') spl &
              & (ist, is)
               Write (*, '(" for species ", I4)') is
               Write (*, '(" and state ", I4)') ist
               Write (*,*)
               Stop
            End If
            If (spk(ist, is) .Lt. 1) Then
               Write (*,*)
               Write (*, '("Error(readinput): spk < 1 : ", I8)') spk &
              & (ist, is)
               Write (*, '(" for species ", I4)') is
               Write (*, '(" and state ", I4)') ist
               Write (*,*)
               Stop
            End If
            If (spocc(ist, is) .Lt. 0.d0) Then
               Write (*,*)
               Write (*, '("Error(readinput): spocc < 0 : ", G18.10)') &
              & spocc (ist, is)
               Write (*, '(" for species ", I4)') is
               Write (*, '(" and state ", I4)') ist
               Write (*,*)
               Stop
            End If
         End Do
         apword (0, is) = speziesdeflist(is)%sp%basis%order
         If (apword(0, is) .Le. 0) Then
            Write (*,*)
            Write (*, '("Error(readinput): apword <= 0 : ", I8)') &
           & apword (0, is)
            Write (*, '(" for species ", I4)') is
            Write (*,*)
            Stop
         End If
         If (apword(0, is) .Gt. maxapword) Then
            Write (*,*)
            Write (*, '("Error(readinput): apword too large : ", I8)') &
           & apword (0, is)
            Write (*, '(" for species ", I4)') is
            Write (*, '("Adjust maxapword in modmain and recompile code&
           &")')
            Write (*,*)
            Stop
         End If
! set the APW orders for l>0
         apword (1:input%groundstate%lmaxapw, is) = apword (0, is)
         Do io = 1, apword (0, is)
            apwe0 (io, 0, is) = &
           & speziesdeflist(is)%sp%basis%wfarray(io)%wf%trialEnergy
            apwdm (io, 0, is) = &
           & speziesdeflist(is)%sp%basis%wfarray(io)%wf%matchingOrder
            apwve (io, 0, is) = &
           & speziesdeflist(is)%sp%basis%wfarray(io)%wf%searchE
            If (apwdm(io, 0, is) .Lt. 0) Then
               Write (*,*)
               Write (*, '("Error(readinput): apwdm < 0 : ", I8)') &
              & apwdm (io, 0, is)
               Write (*, '(" for species ", I4)') is
               Write (*, '(" and order ", I4)') io
               Write (*,*)
               Stop
            End If
! set the APW linearisation energies, derivative orders and variability for l>0
            apwe0 (io, 1:input%groundstate%lmaxapw, is) = apwe0 (io, 0, &
           & is)
            apwdm (io, 1:input%groundstate%lmaxapw, is) = apwdm (io, 0, &
           & is)
            apwve (io, 1:input%groundstate%lmaxapw, is) = apwve (io, 0, &
           & is)
         End Do
         mine0=min(apwe0(io,0,is),mine0)
         nlx = size (speziesdeflist(is)%sp%basis%exceptionarray)
         If (nlx .Lt. 0 .And. &
        & associated(speziesdeflist(is)%sp%basis%exceptionarray)) Then
            Write (*,*)
            Write (*, '("Error(readinput): nlx < 0 : ", I8)') nlx
            Write (*, '(" for species ", I4)') is
            Write (*,*)
            Stop
         End If
         Do ilx = 1, nlx
            lx = speziesdeflist(is)%sp%basis%exceptionarray(ilx)%exception%l
            io = size (speziesdeflist(is)%sp%basis%exceptionarray(ilx)%exception%wfarray)
            If (lx .Lt. 0) Then
               Write (*,*)
               Write (*, '("Error(readinput): lx < 0 : ", I8)') lx
               Write (*, '(" for species ", I4)') is
               Write (*, '(" and exception number ", I4)') ilx
               Write (*,*)
               Stop
            End If
            If (lx .Gt. input%groundstate%lmaxapw) Then
               Write (*,*)
               Write (*, '("Error(readinput): lx > lmaxapw : ", I8)') &
              & lx
               Write (*, '(" for species ", I4)') is
               Write (*, '(" and exception number ", I4)') ilx
               Write (*,*)
               Stop
            End If
            apword (lx, is) = io
            If (apword(lx, is) .Le. 0) Then
               Write (*,*)
               Write (*, '("Error(readinput): apword <= 0 : ", I8)') &
              & apword (lx, is)
               Write (*, '(" for species ", I4)') is
               Write (*, '(" and exception number ", I4)') ilx
               Write (*,*)
               Stop
            End If
            If (apword(lx, is) .Gt. maxapword) Then
               Write (*,*)
               Write (*, '("Error(readinput): apword too large : ", I8)&
              &') apword (lx, is)
               Write (*, '(" for species ", I4)') is
               Write (*, '(" and exception number ", I4)') ilx
               Write (*, '("Adjust maxapword in modmain and recompile c&
              &ode")')
               Write (*,*)
               Stop
            End If
            Do io = 1, apword (lx, is)
!
               apwe0 (io, lx, is) = speziesdeflist(is)%sp%basis%exceptionarray(ilx)%exception%wfarray(io)%wf%trialEnergy
               apwdm (io, lx, is) = speziesdeflist(is)%sp%basis%exceptionarray(ilx)%exception%wfarray(io)%wf%matchingOrder
               apwve (io, lx, is) = speziesdeflist(is)%sp%basis%exceptionarray(ilx)%exception%wfarray(io)%wf%searchE
               If (apwdm(io, lx, is) .Lt. 0) Then
                  Write (*,*)
                  Write (*, '("Error(readinput): apwdm < 0 : ", I8)') &
                 & apwdm (io, lx, is)
                  Write (*, '(" for species ", I4)') is
                  Write (*, '(" exception number ", I4)') ilx
                  Write (*, '(" and order ", I4)') io
                  Write (*,*)
                  Stop
               End If
               mine0=min(apwe0(io,lx,is),mine0)
            End Do
         End Do
         nlorb (is) = size (speziesdeflist(is)%sp%lorbarray)
         If (nlorb(is) .Lt. 0) Then
            Write (*,*)
            Write (*, '("Error(readinput): nlorb < 0 : ", I8)') nlorb &
           & (is)
            Write (*, '(" for species ", I4)') is
            Write (*,*)
            Stop
         End If
         If (nlorb(is) .Gt. maxlorb) Then
            Write (*,*)
            Write (*, '("Error(readinput): nlorb too large : ", I8)') &
           & nlorb (is)
            Write (*, '(" for species ", I4)') is
            Write (*, '("Adjust maxlorb in modmain and recompile code")&
           &')
            Write (*,*)
            Stop
         End If
         Do ilo = 1, nlorb (is)
            lorbl (ilo, is) = &
           & speziesdeflist(is)%sp%lorbarray(ilo)%lorb%l
            lorbord (ilo, is) = size &
           & (speziesdeflist(is)%sp%lorbarray(ilo)%lorb%wfarray)
            If (lorbl(ilo, is) .Lt. 0) Then
               Write (*,*)
               Write (*, '("Error(readinput): lorbl < 0 : ", I8)') &
              & lorbl (ilo, is)
               Write (*, '(" for species ", I4)') is
               Write (*, '(" and local-orbital ", I4)') ilo
               Write (*,*)
               Stop
            End If
            If (lorbl(ilo, is) .Gt. input%groundstate%lmaxmat) Then
               Write (*,*)
               Write (*, '("Error(readinput): lorbl > lmaxmat : ", 2I8)&
              &') lorbl (ilo, is), input%groundstate%lmaxmat
               Write (*, '(" for species ", I4)') is
               Write (*, '(" and local-orbital ", I4)') ilo
               Write (*,*)
               Stop
            End If
            If (lorbord(ilo, is) .Lt. 2) Then
               Write (*,*)
               Write (*, '("Error(readinput): lorbord < 2 : ", I8)') &
              & lorbord (ilo, is)
               Write (*, '(" for species ", I4)') is
               Write (*, '(" and local-orbital ", I4)') ilo
               Write (*,*)
               Stop
            End If
            If (lorbord(ilo, is) .Gt. maxlorbord) Then
               Write (*,*)
               Write (*, '("Error(readinput): lorbord too large : ", I8&
              &)') lorbord (ilo, is)
               Write (*, '(" for species ", I4)') is
               Write (*, '(" and local-orbital ", I4)') ilo
               Write (*, '("Adjust maxlorbord in modmain and recompile &
              &code")')
               Write (*,*)
               Stop
            End If
            Do io = 1, lorbord (ilo, is)
               lorbe0 (io, ilo, is) = speziesdeflist(is)%sp%lorbarray(ilo)%lorb%wfarray(io)%wf%trialEnergy
               lorbdm (io, ilo, is) = speziesdeflist(is)%sp%lorbarray(ilo)%lorb%wfarray(io)%wf%matchingOrder
               lorbve (io, ilo, is) = speziesdeflist(is)%sp%lorbarray(ilo)%lorb%wfarray(io)%wf%searchE
               If (lorbdm(io, ilo, is) .Lt. 0) Then
                  Write (*,*)
                  Write (*, '("Error(readinput): lorbdm < 0 : ", I8)') &
                 & lorbdm (io, ilo, is)
                  Write (*, '(" for species ", I4)') is
                  Write (*, '(" local-orbital ", I4)') ilo
                  Write (*, '(" and order ", I4)') io
                  Write (*,*)
                  Stop
               End If
               mine0=min(lorbe0(io,ilo,is),mine0)
            End Do
         End Do
         Close (50)
      End Do
      Return
End Subroutine
