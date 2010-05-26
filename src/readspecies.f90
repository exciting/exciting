!
!
!
! Copyright (C) 2002-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine readspecies
      Use modmain
      Use modinput
      Implicit None
! local variables
      Integer :: is, ist, iostat
      Integer :: io, nlx, ilx, lx, ilo
      mine0=0.d0
      Do is = 1, nspecies
         Open (50, File=trim(input%structure%speciespath)//&
        & trim(input%structure%speciesarray(is)%species%speciesfile), &
        & Action='READ', Status='OLD', Form='FORMATTED', IoStat=IoStat)
         If (iostat .Ne. 0) Then
            Write (*,*)
            Write (*, '("Error(readinput): error opening species file "&
           &, A)') trim (input%structure%speciespath) // trim &
           & (input%structure%speciesarray(is)%species%speciesfile)
            Write (*,*)
            Stop
         End If
         Read (50,*) &
        & input%structure%speciesarray(is)%species%chemicalSymbol
         Read (50,*) spname (is)
         Read (50,*) spzn (is)
         Read (50,*) spmass (is)
         Read (50,*) sprmin (is), rmt (is), sprmax (is), nrmt (is)
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
         Read (50,*) spnst (is)
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
            Read (50,*) spn (ist, is), spl (ist, is), spk (ist, is), &
           & spocc (ist, is), spcore (ist, is)
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
         Read (50,*) apword (0, is)
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
            Read (50,*) apwe0 (io, 0, is), apwdm (io, 0, is), apwve &
           & (io, 0, is)
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
           mine0=min(apwe0(io,0,is),mine0)
         End Do
         Read (50,*) nlx
         If (nlx .Lt. 0) Then
            Write (*,*)
            Write (*, '("Error(readinput): nlx < 0 : ", I8)') nlx
            Write (*, '(" for species ", I4)') is
            Write (*,*)
            Stop
         End If
         Do ilx = 1, nlx
            Read (50,*) lx, io
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
               Read (50,*) apwe0 (io, lx, is), apwdm (io, lx, is), &
              & apwve (io, lx, is)
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
         Read (50,*) nlorb (is)
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
            Read (50,*) lorbl (ilo, is), lorbord (ilo, is)
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
               Read (50,*) lorbe0 (io, ilo, is), lorbdm (io, ilo, is), &
              & lorbve (io, ilo, is)
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
