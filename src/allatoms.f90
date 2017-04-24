!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: allatoms
! !INTERFACE:
!
!
Subroutine allatoms(verbosity)
! !USES:
      Use modinput
      Use modmain
      Use FoX_wxml
      Use modmpi, only : rank, finitMPI
! !DESCRIPTION:
!   Solves the Kohn-Sham-Dirac equations for each atom type in the solid and
!   finds the self-consistent radial wavefunctions, eigenvalues, charge
!   densities and potentials. The atomic densities can then be used to
!   initialise the crystal densities, and the atomic self-consistent potentials
!   can be appended to the muffin-tin potentials to solve for the core states.
!   Note that, irrespective of the value of {\tt xctype}, exchange-correlation
!   functional type 3 is used. See also {\tt atoms}, {\tt rhoinit},
!   {\tt gencore} and {\tt modxcifc}.
!
! !REVISION HISTORY:
!   Created September 2002 (JKD)
!   Modified for GGA, June 2007 (JKD)
!   Modified for DFT-1/2, 2015 (Ronaldo)
!EOP
!BOC
      Implicit None
      integer :: verbosity
! always use LDA to setup atomic densities
      Integer, Parameter :: xctype_ = 3
      Integer :: xctypearray(3)
      Integer, Parameter :: xcgrad_ = 0
      Integer :: fnum_ = 333 ! file where Vs (see DFT-1/2 details) is written
      Integer :: is, i, ir, n, nshell
      Integer, Allocatable :: shell(:)
      Logical :: dirac_eq
      character(100)::buffer
      Character (256) :: fname
      Real (8) :: ampl, cut, cutfunction, aux
! allocatable arrays
      Real (8), Allocatable :: rwf (:, :, :)
      Real (8), Allocatable :: newspocc(:)
      Real (8), Allocatable :: ionization(:)
      Real (8), Allocatable :: rhoslave(:)
      Type (xmlf_t), Save :: xf
      
      dirac_eq=(input%groundstate%CoreRelativity.eq."dirac")
! allocate global species charge density and potential arrays
      If (allocated(sprho)) deallocate (sprho)
      Allocate (sprho(spnrmax, nspecies))
      If (allocated(spvr)) deallocate (spvr)
      Allocate (spvr(spnrmax, nspecies))
!write(*,*) 'howdy'
!     We need to allocate some arrays for the DFT-1/2 part
      if (associated(input%groundstate%dfthalf)) then
        if (allocated(vhalfir)) deallocate (vhalfir)
        allocate(vhalfir(ngrtot))  
        if (allocated(vhalfmt)) deallocate (vhalfmt)    
        allocate(vhalfmt(lmmaxvr, nrmtmax, natmtot))
        if (allocated(vhalfsph)) deallocate (vhalfsph)
        allocate(vhalfsph(spnrmax,nspecies))        
      endif
!write(*,*) 'approaching the loop'

      ! All libxc routines expect xctype(3) instead of a single integer
      xctypearray(1:3) = xctype_
#ifdef USEOMPallatoms
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(rwf,fnum_,fname,i,ir,n,ampl,cut,cutfunction,aux,nshell,shell,ionization,newspocc,rhoslave)
!$OMP DO
#endif
      Do is = 1, nspecies
         Allocate (rwf(spnrmax, 2, spnstmax))
         Call atom (input%groundstate%ptnucl, spzn(is), spnst(is), &
        & spn(:, is), spl(:, is), spk(:, is), spocc(:, is), xctypearray, &
        & xcgrad_, spnr(is), spr(:, is), &
        & speval(:, is), sprho(:, is), spvr(:, is), rwf,nrmt(is),dirac_eq)
         Deallocate (rwf)
         if (associated(input%groundstate%dfthalf)) then
!          Here, we check if the user has defined the DFT-1/2 correction for this species
!          We allways make the DFT correction (when it is triggered), even for the species for
!          which the corrections are not being required. The trick is to define the DFT-1/2 correction 
!          (when not specified for the species) aiming to imply a null correction (this can be handled with CUT=0)
           if (associated(input%structure%speciesarray(is)%species%dfthalfparam)) then 
             n = input%structure%speciesarray(is)%species%dfthalfparam%exponent
             ampl = input%structure%speciesarray(is)%species%dfthalfparam%ampl
             cut = input%structure%speciesarray(is)%species%dfthalfparam%cut
             nshell = size (input%structure%speciesarray(is)%species%dfthalfparam%shellarray)
             If (allocated(ionization)) deallocate (ionization) 
             Allocate (ionization(nshell))
             If (allocated(shell)) deallocate (shell) 
             Allocate (shell(nshell))
             Do i = 1, nshell
               shell(i) = input%structure%speciesarray(is)%species%dfthalfparam%shellarray(i)%shell%number
               ionization(i) = input%structure%speciesarray(is)%species%dfthalfparam%shellarray(i)%shell%ionization
             Enddo
           else
             n = 8
             ampl = 1.0
             cut = 0.0
             nshell = 1
             Allocate (ionization(nshell))
             Allocate (shell(nshell))
             shell(1) = 0
             ionization(1) = 0.5
           endif
           fnum_ = 333 + is
!          Open the file where the information about Vs will be outputed
!          E.g., the name is something like VS_S0001.OUT for the first species
           if ((input%groundstate%dfthalf%printVSfile) .and.(rank .eq. 0) ) then
             Write (fname, '("VS_S", I2.2, ".OUT")') is
             Open (fnum_, File=trim(fname), Action='WRITE', Form='FORMATTED')
             Write(fnum_,'("Species: ",I4.1,", CUT = ",F10.6,", Amplitude: ",F10.6,", Exponent: ",I2.2)') is, cut, ampl, n
             Write(fnum_,'(I4.1," shell(s) must be ionized. Its/Their ionization(s) are listed below")') nshell
             Do i = 1, nshell
               Write(fnum_,'(I4.1, F10.6)') shell(i), ionization(i)
             Enddo
             Write(fnum_,'(A20,A20,A20,A20,A20)') "r","Vatom","Vion","CUT-function","VS"
           end if
           If (allocated(rwf)) deallocate (rwf)
           Allocate (rwf(spnrmax, 2, spnstmax))
           If (allocated(rhoslave)) deallocate (rhoslave)
           Allocate ( rhoslave(spnr(is)) )
!          Now, for each species, we define a new occupation 
!          (concerning the Ionization scheme provided by the user)
           If (allocated(newspocc)) deallocate (newspocc)
           Allocate ( newspocc(spnst(is)) )
           newspocc(1:spnst(is)) = spocc(1:spnst(is), is)
           Do i = 1, nshell
             if (shell(i).le.0) then
               shell(i) = spnst(is)+shell(i)
             endif
           Enddo
           Do i = 1, nshell
             if ((shell(i).lt.0).or.(shell(i).gt.spnst(is))) then
               write(*,*) 'Error concerning the shell%number of species ', is
               write(*,*) 'It should not be negative or higher than the number of total shells'
               stop
             endif
           Enddo
!          Here, the occupation (considering the ionization) is really implemented
           Do i = 1, nshell
             newspocc(shell(i)) = newspocc(shell(i)) - ionization(i)
             if ((newspocc(shell(i))).lt.0) then
               write(*,*) 'Error concerning the ionization defined for species ', is
               write(*,*) 'A negatively occupied shell results from this ionization condition'
               stop
             endif
           Enddo
!          Now, we calculate the ionized atom
           Call atom (input%groundstate%ptnucl, spzn(is), spnst(is), &
             & spn(:, is), spl(:, is), spk(:, is), newspocc(:), xctypearray, &
             & xcgrad_, spnr(is), spr(:, is), &
             & speval(:, is), rhoslave, vhalfsph(:, is), rwf,nrmt(is), &
             & dirac_eq)
           vhalfsph(:,is) = vhalfsph(:,is)-spvr(:,is)
           Do ir = 1, spnr(is)
!          Now we can introduce the cut-function
!          This 'aux' variable just stores the value of VKS which came from the ionized atom
             aux = vhalfsph(ir,is) + spvr(ir,is)
             if (spr(ir, is).lt.cut) then
               cutfunction = ampl*( 1-((spr(ir, is)/cut)**n) )**3
               vhalfsph(ir,is) = vhalfsph(ir,is)*cutfunction
             else
               cutfunction = 0.0
               vhalfsph(ir,is) = 0.0
             endif
             if ((input%groundstate%dfthalf%printVSfile) .and.(rank .eq. 0)) then
               write (fnum_,"(ES20.6E2,ES20.6E2,ES20.6E2,F20.10,F20.10)") spr(ir,is), spvr(ir,is), aux, cutfunction, vhalfsph(ir,is)
             end if
           Enddo
           If (allocated(newspocc)) Deallocate (newspocc)
           If (allocated(rhoslave)) Deallocate (rhoslave)
           If (allocated(rwf)) Deallocate (rwf)
           If (allocated(shell)) Deallocate (shell)
           If (allocated(ionization)) Deallocate (ionization)
           if (input%groundstate%dfthalf%printVSfile) then
             Close(fnum_)
           end if
         Endif !if (associated(input%groundstate%dfthalf))
      End Do
#ifdef USEOMPallatoms
!$OMP END DO
!$OMP END PARALLEL
#endif

      if (associated(input%groundstate%dfthalf)) then 
        if (input%groundstate%dfthalf%printVSfile) then
#ifdef MPI
        call finitMPI()
#endif
        stop
        end if
      end if

!     Now, let's expand our V_S potential in the MT and in the interstitial part
      if (associated(input%groundstate%dfthalf)) then
        call vhalfinit
      endif
      
      if ((verbosity.gt.0).and.(rank.eq.0)) then
        Call xml_OpenFile ("atoms.xml", xf, replace=.True., pretty_print=.True.)
        Call xml_NewElement(xf,"atomlist")
        Call xml_NewElement (xf,"Hamiltonian")
        Call xml_AddAttribute (xf,"RelativityModel",trim(input%groundstate%CoreRelativity))
        Call xml_AddAttribute (xf,"xctype",xctype_)
        Call xml_EndElement (xf,"Hamiltonian")
        Do is = 1, nspecies
          Call xml_NewElement (xf,"atom")
          Call xml_AddAttribute (xf,"chemicalSymbol", trim(input%structure%speciesarray(is)%species%chemicalSymbol)) 
          Call xml_AddAttribute (xf,"species", trim(input%structure%speciesarray(is)%species%speciesfile))
          Call xml_NewElement (xf,"NumericalSetup")
          Call xml_AddAttribute (xf,"TotalNumberOfGridPoints",spnr(is))
          Call xml_AddAttribute (xf,"NumberOfMTGridPoints",nrmt(is))
          Call xml_AddAttribute (xf,"GridType",trim(input%groundstate%radialgridtype))
          Call xml_AddAttribute (xf,"rmin",spr(1, is))
          Call xml_AddAttribute (xf,"rmt",spr(nrmt(is), is))
          Call xml_AddAttribute (xf,"rmax",spr(spnr(is), is))
          Call xml_EndElement (xf,"NumericalSetup")
          Call xml_NewElement (xf,"spectrum")
          do i=1,spnst(is)
            Call xml_NewElement (xf,"state")
            Call xml_AddAttribute (xf,"n",spn(i, is))
            Call xml_AddAttribute (xf,"l",spl(i, is))
            Call xml_AddAttribute (xf,"kappa",spk(i, is))
            write(buffer,'(G22.12)') speval(i,is)
            Call xml_AddAttribute (xf,"energy",trim(adjustl(buffer)))
            Call xml_EndElement (xf,"state")
          enddo
          Call xml_EndElement (xf,"spectrum")
          Call xml_EndElement (xf,"atom")
        End Do
        Call xml_EndElement (xf,"atomlist")
        Call xml_close (xf)
      endif
      Return
End Subroutine
!EOC
