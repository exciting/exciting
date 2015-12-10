!
Subroutine readspeciesxml
  
  Use modmain
  Use modinput, Only: input
  Use modspdeflist
  Use FoX_dom
  Use modspdb
  Use modmpi

  Implicit None
! local variables
  Integer :: is, ist
  Integer :: io, nlx, ilx, lx, nlo, ilo, errorcode
  Type (Node), Pointer :: speciesnp, speciesdbnp
  character(2048) :: command
  character(256)  :: spfile_string
  character(256)  :: string

  if (allocated(speziesdeflist)) deallocate(speziesdeflist)
  Allocate(speziesdeflist(nspecies))
  config => newDOMConfig ()
  parseerror = .False.

! parse xml and create derived type for species definitions  speziesdeflist
  Do is = 1, nspecies
    spfile_string=""
    
    if (index(input%structure%speciespath,"http://").ge.1) then

#ifdef CURL
      write(command,*)"curl ",trim(input%structure%speciespath),          &
        trim(input%structure%speciesarray(is)%species%speciesfile)," > ", &
        trim(input%structure%speciesarray(is)%species%speciesfile)
#endif
#ifndef CURL
      write(command,*)"wget -c ",trim(input%structure%speciespath), &
        trim(input%structure%speciesarray(is)%species%speciesfile)
#endif

      write(spfile_string,*)"./", trim(input%structure%speciesarray(is)%species%speciesfile)
      If (rank .Eq. 0) Then
        call system(command)
      endif
#ifdef MPI
      Call MPI_barrier (MPI_COMM_WORLD, ierr)
#endif

    else
    
      write(spfile_string,*) trim(input%structure%speciespath) // "/" &
      & // trim(input%structure%speciesarray(is)%species%speciesfile)
      
      !-------------------------------------------------
      ! Copy all species files to the working directory
      !-------------------------------------------------
      if (rank==0) then
        if ((trim(input%structure%speciespath)/='./').and. &
        &   (trim(input%structure%speciespath)/='.')) then
          write(command,*) "cp" // trim(spfile_string) // " ."
          !write(*,*) trim(command)
          call system(trim(command))
        end if
      end if
      call barrier
  
    endif

    doc => parseFile(ADJUSTL(trim(spfile_string)),config,iostat=errorcode)

    if(errorcode.ne.0) then
      write(*,*) "### Could not open ", ADJUSTL(trim(spfile_string))
      write(*,*) "### Check if file exists and if it is well-formed XML."
      stop
    endif
    
    speciesdbnp => getDocumentElement(doc)
    speciesnp => item(getElementsByTagname(speciesdbnp, "sp"), 0)
    parseerror = .False.
    speziesdeflist(is)%sp => getstructsp (speciesnp)
    Call destroy (doc)
  End Do
!
!
!
  Do is = 1, nspecies
     spsymb(is) = trim(speziesdeflist(is)%sp%chemicalSymbol)
     input%structure%speciesarray(is)%species%chemicalSymbol = spsymb (is)
     spname(is) = trim(speziesdeflist(is)%sp%name)
     spzn(is) = speziesdeflist(is)%sp%z
     spmass(is) = speziesdeflist(is)%sp%mass
     sprmin(is) = speziesdeflist(is)%sp%muffinTin%rmin
     rmt(is) = speziesdeflist(is)%sp%muffinTin%radius
     If (input%structure%speciesarray(is)%species%rmt .Gt. 0) Then
        rmt(is) = input%structure%speciesarray(is)%species%rmt
     End If
     sprmax(is) = speziesdeflist(is)%sp%muffinTin%rinf
     nrmt(is) = speziesdeflist(is)%sp%muffinTin%radialmeshPoints
     If (sprmin(is) .Le. 0.d0) Then
        Write (*,*)
        Write (*, '("Error(readinput): sprmin <= 0 : ", G18.10)') sprmin(is)
        Write (*, '(" for species ", I4)') is
        Write (*,*)
        Stop
     End If
     If (rmt(is) .Le. sprmin(is)) Then
        Write (*,*)
        Write (*, '("Error(readinput): rmt <= sprmin : ", 2G18.10)') rmt(is), sprmin(is)
        Write (*, '(" for species ", I4)') is
        Write (*,*)
        Stop
     End If
     If (sprmax(is) .Lt. rmt(is)) Then
        Write (*,*)
        Write (*, '("Error(readinput): sprmax < rmt : ", 2G18.10)') sprmax(is), rmt(is)
        Write (*,*)
        Stop
     End If
     If (nrmt(is) .Lt. 20) Then
        Write (*,*)
        Write (*, '("Error(readinput): nrmt too small : ", I8)') nrmt(is)
        Write (*, '(" for species ", I4)') is
        Write (*,*)
        Stop
     End If
     spnst(is) = size(speziesdeflist(is)%sp%AtomicStatearray)
     If (spnst(is) .Le. 0) Then
        Write (*,*)
        Write (*, '("Error(readinput): invalid spnst : ", I8)') spnst(is)
        Write (*, '(" for species ", I4)') is
        Write (*,*)
        Stop
     End If
     If (spnst(is) .Gt. maxspst) Then
        Write (*,*)
        Write (*, '("Error(readinput): too many states for species", I8)') is
        Write (*,*)
        Stop
     End If
     Do ist = 1, spnst (is)
        spn(ist, is) = &
       &  speziesdeflist(is)%sp%AtomicStatearray(ist)%atomicState%n
        spl(ist, is) = &
       &  speziesdeflist(is)%sp%AtomicStatearray(ist)%atomicState%l
        spk(ist, is) = &
       &  speziesdeflist(is)%sp%AtomicStatearray(ist)%atomicState%kappa
        spocc(ist, is) = &
       &  speziesdeflist(is)%sp%AtomicStatearray(ist)%atomicState%occ
        spcore(ist, is) = &
       &  speziesdeflist(is)%sp%AtomicStatearray(ist)%atomicState%core
        If (spn(ist, is) .Lt. 1) Then
           Write (*,*)
           Write (*, '("Error(readinput): spn < 1 : ", I8)') spn(ist, is)
           Write (*, '(" for species ", I4)') is
           Write (*, '(" and state ", I4)') ist
           Write (*,*)
           Stop
        End If
        If (spl(ist, is) .Lt. 0) Then
           Write (*,*)
           Write (*, '("Error(readinput): spl < 0 : ", I8)') spl(ist, is)
           Write (*, '(" for species ", I4)') is
           Write (*, '(" and state ", I4)') ist
           Write (*,*)
           Stop
        End If
        If (spk(ist, is) .Lt. 1) Then
           Write (*,*)
           Write (*, '("Error(readinput): spk < 1 : ", I8)') spk(ist, is)
           Write (*, '(" for species ", I4)') is
           Write (*, '(" and state ", I4)') ist
           Write (*,*)
           Stop
        End If
        If (spocc(ist, is) .Lt. 0.d0) Then
           Write (*,*)
           Write (*, '("Error(readinput): spocc < 0 : ", G18.10)') spocc (ist, is)
           Write (*, '(" for species ", I4)') is
           Write (*, '(" and state ", I4)') ist
           Write (*,*)
           Stop
        End If
     End Do

!----------------------------
!    Default definitions
!----------------------------

     if (size(speziesdeflist(is)%sp%basis%default%wfarray)>0) then
       
!      DEFAULT: Element wf is specified
!      The definition based on the augmentation type is (if present) ignored

       apword(0, is) = size(speziesdeflist(is)%sp%basis%default%wfarray)
       If (apword(0, is) .Le. 0) Then
          Write (*,*)
          Write (*, '("Error(readinput): apword <= 0 : ", I8)') apword(0, is)
          Write (*, '(" for species ", I4)') is
          Write (*,*)
          Stop
       End If
       If (apword(0, is) .Gt. maxapword) Then
          Write (*,*)
          Write (*, '("Error(readinput): apword too large : ", I8)') apword (0, is)
          Write (*, '(" for species ", I4)') is
          Write (*, '("Adjust maxapword in modmain and recompile code")')
          Write (*,*)
          Stop
       End If

!      set the APW orders for l>0
       apword(1:input%groundstate%lmaxapw, is) = apword (0, is)
       Do io = 1, apword (0, is)
          apwe0(io, 0, is) = &
         & speziesdeflist(is)%sp%basis%default%wfarray(io)%wf%trialEnergy
          apwdm(io, 0, is) = &
         & speziesdeflist(is)%sp%basis%default%wfarray(io)%wf%matchingOrder
          apwve(io, 0, is) = &
         & speziesdeflist(is)%sp%basis%default%wfarray(io)%wf%searchE
          If (apwdm(io, 0, is) .Lt. 0) Then
             Write (*,*)
             Write (*, '("Error(readinput): apwdm < 0 : ", I8)') apwdm (io, 0, is)
             Write (*, '(" for species ", I4)') is
             Write (*, '(" and order ", I4)') io
             Write (*,*)
             Stop
          End If
       End Do

     else
       
       string=trim(speziesdeflist(is)%sp%basis%default%type)
       if ((string .eq. 'apw+lo').or.(string .eq. 'apw')) then
         apword(0, is) = 1
       else if (string .eq. 'lapw') then
         apword(0, is) = 2
       else
         write (*, *) 'Error(readspeciesxml): Unknown APW basis type = ', adjustl(string)
         write (*, *) 'Available options: lapw, apw, apw+lo'
         stop
       end if
       do io = 1, apword(0, is)                                              
          apwe0(io, 0, is) = speziesdeflist(is)%sp%basis%default%trialEnergy 
          apwdm(io, 0, is) = io-1                                            
          apwve(io, 0, is) = speziesdeflist(is)%sp%basis%default%searchE
       end do    

     end if

!    the APW orders, linearisation energies, derivative orders and variability for l>0
     apword(1:input%groundstate%lmaxapw, is) = apword(0, is)
     do io = 1, apword(0, is)
        apwe0(io, 1:input%groundstate%lmaxapw, is) = apwe0(io, 0, is)
        apwdm(io, 1:input%groundstate%lmaxapw, is) = apwdm(io, 0, is)
        apwve(io, 1:input%groundstate%lmaxapw, is) = apwve(io, 0, is)
     end do

!----------------------------
!    Custom definitions
!----------------------------

     mine0 = 0.d0
     mine0 = min(apwe0(io,0,is),mine0)
     nlx = size(speziesdeflist(is)%sp%basis%customarray)
     If ( (nlx .Lt. 0) .And. (associated(speziesdeflist(is)%sp%basis%customarray)) ) Then
        Write (*,*)
        Write (*, '("Error(readinput): nlx < 0 : ", I8)') nlx
        Write (*, '(" for species ", I4)') is
        Write (*,*)
        Stop
     End If

     nlo = 0
     Do ilx = 1, nlx
        
        lx = speziesdeflist(is)%sp%basis%customarray(ilx)%custom%l
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
        
        if (size(speziesdeflist(is)%sp%basis%customarray(ilx)%custom%wfarray)>0) then

!         CUSTOM: Element wf is specified
!         The definition bases on the augmentation type is (if present) ignored
        
          io = size (speziesdeflist(is)%sp%basis%customarray(ilx)%custom%wfarray)
          apword(lx, is) = io
          If (apword(lx, is) .Le. 0) Then
             Write (*,*)
             Write (*, '("Error(readinput): apword <= 0 : ", I8)') apword (lx, is)
             Write (*, '(" for species ", I4)') is
             Write (*, '(" and exception number ", I4)') ilx
             Write (*,*)
             Stop
          End If
          If (apword(lx, is) .Gt. maxapword) Then
             Write (*,*)
             Write (*, '("Error(readinput): apword too large : ", I8)') apword (lx, is)
             Write (*, '(" for species ", I4)') is
             Write (*, '(" and exception number ", I4)') ilx
             Write (*, '("Adjust maxapword in modmain and recompile code")')
             Write (*,*)
             Stop
          End If
          Do io = 1, apword (lx, is)
             apwe0(io, lx, is) = speziesdeflist(is)%sp%basis%customarray(ilx)%custom%wfarray(io)%wf%trialEnergy
             apwdm(io, lx, is) = speziesdeflist(is)%sp%basis%customarray(ilx)%custom%wfarray(io)%wf%matchingOrder
             apwve(io, lx, is) = speziesdeflist(is)%sp%basis%customarray(ilx)%custom%wfarray(io)%wf%searchE
             If (apwdm(io, lx, is) .Lt. 0) Then
                Write (*,*)
                Write (*, '("Error(readinput): apwdm < 0 : ", I8)') apwdm (io, lx, is)
                Write (*, '(" for species ", I4)') is
                Write (*, '(" exception number ", I4)') ilx
                Write (*, '(" and order ", I4)') io
                Write (*,*)
                Stop
             End If
             mine0=min(apwe0(io,lx,is),mine0)
          End Do
        
        else
        
          string = trim(speziesdeflist(is)%sp%basis%customarray(ilx)%custom%type)
          if ((string .eq. 'apw+lo').or.(string .eq. 'apw')) then
            apword(lx, is) = 1
          else if (string .eq. 'lapw') then
            apword(lx, is) = 2
          else
            write (*, *) 'Error(readspeciesxml): Unknown APW basis type = ', adjustl(string)
            write (*, *) 'Available options: lapw, apw, apw+lo'
            stop
          end if
          Do io = 1, apword(lx, is)
             apwe0(io, lx, is) = speziesdeflist(is)%sp%basis%customarray(ilx)%custom%trialEnergy
             apwdm(io, lx, is) = io-1
             apwve(io, lx, is) = speziesdeflist(is)%sp%basis%customarray(ilx)%custom%searchE
             mine0=min(apwe0(io,lx,is),mine0)
          End Do
!        
!         lo in APW+lo method
!
          if (string .eq. 'apw+lo') then
            nlo = nlo+1
            lorbl(nlo, is) = lx
            if (lorbl(nlo, is) .Gt. input%groundstate%lmaxmat) then
              write (*,*)
              write (*, '("Error(readinput): lorbl > lmaxmat : ", 2I8)') &
             &  lorbl(nlo, is), input%groundstate%lmaxmat
              write (*, '(" for species ", I4)') is
              write (*,*)
              stop
            end if
            lorbord(nlo, is) = 2
            do io = 1, lorbord(nlo, is)
              lorbe0(io, nlo, is) = apwe0(1, lx, is)
              lorbdm(io, nlo, is) = io-1
              lorbve(io, nlo, is) = apwve(1, lx, is)
            end do
            mine0=min(lorbe0(io,nlo,is),mine0)
          end if
        
        end if

     End Do ! ilx

!-----------------------------------
!    Local orbitals definitions
!-----------------------------------

     nlorb(is) = nlo+size(speziesdeflist(is)%sp%basis%loarray)
     If (nlorb(is) .Gt. maxlorb) Then
        Write (*,*)
        Write (*, '("Error(readinput): nlorb too large : ", I8)') nlorb(is)
        Write (*, '(" for species ", I4)') is
        Write (*, '("Adjust maxlorb in modmain and recompile code")')
        Write (*,*)
        Stop
     End If

     Do ilx = 1, size(speziesdeflist(is)%sp%basis%loarray)
        ilo = nlo+ilx
        lorbl(ilo, is) = speziesdeflist(is)%sp%basis%loarray(ilx)%lo%l
        lorbord(ilo, is) = size(speziesdeflist(is)%sp%basis%loarray(ilx)%lo%wfarray)
        If (lorbl(ilo, is) .Lt. 0) Then
           Write (*,*)
           Write (*, '("Error(readinput): lorbl < 0 : ", I8)') lorbl(ilo, is)
           Write (*, '(" for species ", I4)') is
           Write (*, '(" and local-orbital ", I4)') ilo
           Write (*,*)
           Stop
        End If
        If (lorbl(ilo, is) .Gt. input%groundstate%lmaxmat) Then
           Write (*,*)
           Write (*, '("Error(readinput): lorbl > lmaxmat : ", 2I8)') &
          &  lorbl(ilo, is), input%groundstate%lmaxmat
           Write (*, '(" for species ", I4)') is
           Write (*, '(" and local-orbital ", I4)') ilo
           Write (*,*)
           Stop
        End If
        If (lorbord(ilo, is) .Lt. 2) Then
           Write (*,*)
           Write (*, '("Warning(readinput): lorbord < 2 : ", I8)') lorbord (ilo, is)
           Write (*, '(" for species ", I4)') is
           Write (*, '(" and local-orbital ", I4)') ilo
           Write (*,*)
        End If
        If (lorbord(ilo, is) .Gt. maxlorbord) Then
           Write (*,*)
           Write (*, '("Error(readinput): lorbord too large : ", I8)') lorbord(ilo, is)
           Write (*, '(" for species ", I4)') is
           Write (*, '(" and local-orbital ", I4)') ilo
           Write (*, '("Adjust maxlorbord in modmain and recompile code")')
           Write (*,*)
           Stop
        End If
        Do io = 1, lorbord(ilo, is)
           lorbe0(io, ilo, is) = speziesdeflist(is)%sp%basis%loarray(ilx)%lo%wfarray(io)%wf%trialEnergy
           lorbdm(io, ilo, is) = speziesdeflist(is)%sp%basis%loarray(ilx)%lo%wfarray(io)%wf%matchingOrder
           lorbve(io, ilo, is) = speziesdeflist(is)%sp%basis%loarray(ilx)%lo%wfarray(io)%wf%searchE
           If (lorbdm(io, ilo, is) .Lt. 0) Then
              Write (*,*)
              Write (*, '("Error(readinput): lorbdm < 0 : ", I8)') lorbdm(io, ilo, is)
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
