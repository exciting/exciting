
! Copyright (C) 2005-2010 C. Meisenbichler and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

Module scl_xml_out_Module
      Use FoX_dom
      Use mod_energy
      Use mod_LDA_LU
      Use mod_convergence
      Use mod_eigenvalue_occupancy
      Use mod_charge_and_moment
      Use mod_atoms
      Use mod_timing
      Use mod_force
      Use mod_spin
      Use mod_misc, Only: githash, filext
      Use modmpi, Only: rank
!
      Use modinput
      Implicit None
      Type (Node), Pointer :: sclDoc, root, np, npatt, energies, niter, &
     & charges, atom, xst, timing, nscl, ngroundstate
      Type (DOMConfiguration), Pointer :: configo
      Real (8) :: scltime0 = 0.d0
      real(8) :: deltae, dforcemax
      Character (512) :: buffer
!
Contains
!
      Subroutine scl_xml_out_create
         Character (10) :: dat, tim
         If (rank .Eq. 0) Then
       ! Create a new document and get a pointer to the root element, this gives you the minimum empty dom
            sclDoc => createDocument (getImplementation(), "", "info", &
           & null())
            configo => getDomConfig (sclDoc)
            Call setParameter (getDomConfig(sclDoc), "format-pretty-pri&
           &nt", .True.)
            root => getDocumentElement (sclDoc)
            xst => createProcessingInstruction (sclDoc, "xml-stylesheet&
           &", 'href="'//trim(input%xsltpath)//'/info.xsl" type="text/x&
           &sl"')
            dummy => insertBefore (sclDoc, xst, root)
!
            ngroundstate => createElementNS (sclDoc, "", "groundstate")
            dummy => appendChild (root, ngroundstate)
            nscl => createElementNS (sclDoc, "", "scl")
            dummy => appendChild (ngroundstate, nscl)
            Call date_and_time (date=dat, time=tim)
!
            Write (buffer, '(A4, "-", A2, "-", A2)') dat (1:4), dat &
           & (5:6), dat (7:8)
            Call setAttribute (root, "date", trim(adjustl(buffer)))
            Write (buffer, '( A2, ":", A2, ":", A2)') tim (1:2), tim &
           & (3:4), tim (5:6)
            Call setAttribute (root, "time", trim(adjustl(buffer)))
            buffer = githash
            Call setAttribute (root, "versionhash", &
           & trim(adjustl(buffer)))
            Call setAttribute (root, "title", &
           & trim(adjustl(input%title)))
!
         End If
!
      End Subroutine scl_xml_out_create

      Subroutine structure_xmlout
         Use mod_lattice
         Use mod_constants
         Implicit None
         Integer :: is, ia, ias, i
         Type (Node), Pointer :: structure, crystal, nbasevect, &
        & nreziprvect, species, atom, forces, text, force
         If (rank .Eq. 0) Then
            structure => createElementNS (sclDoc, "", "structure")
            dummy => appendChild (nscl, structure)
            crystal => createElementNS (sclDoc, "", "crystal")
            dummy => appendChild (structure, crystal)
            Write (buffer, '(G18.10)') omega
            Call setAttribute (crystal, "unitCellVolume", &
           & trim(adjustl(buffer)))
            Write (buffer, '(G18.10)') (twopi**3) / omega
            Call setAttribute (crystal, "BrillouinZoneVolume", &
           & trim(adjustl(buffer)))
            Do i = 1, 3
               nbasevect => createElementNS (sclDoc, "", "basevect")
               dummy => appendChild (crystal, nbasevect)
               Write (buffer, '(3G18.10)') &
              & input%structure%crystal%basevect(:, i)
               text => createTextNode (sclDoc, trim(adjustl(buffer)))
               dummy => appendChild (nbasevect, text)
            End Do
            Do i = 1, 3
               nreziprvect => createElementNS (sclDoc, "", "reciprvect")
               Write (buffer, '(3G18.10)') bvec (:, i)
               text => createTextNode (sclDoc, trim(adjustl(buffer)))
               dummy => appendChild (nreziprvect, text)
               dummy => appendChild (crystal, nreziprvect)
            End Do
            Write (buffer,'(I5)')  input%groundstate%nktot
            Call setAttribute (crystal, "nktot", &
            & trim(adjustl(buffer)))            
            Write (buffer, '(3I5)')  input%groundstate%ngridk(:)
            Call setAttribute (crystal, "ngridk", &
            & trim(adjustl(buffer)))            
            Do is = 1, nspecies
               species => createElementNS (sclDoc, "", "species")
               dummy => appendChild (structure, species)
               Write (buffer,*) trim (input%structure%speciesarray(is)%species%chemicalSymbol)
               Call setAttribute (species, "chemicalSymbol", &
              & trim(adjustl(buffer)))
               Do ia = 1, natoms (is)
                  atom => createElementNS (sclDoc, "", "atom")
                  dummy => appendChild (species, atom)
                  Call setcoord (atom, input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord)
                  If (input%groundstate%tforce) Then
                     forces => createElementNS (sclDoc, "", "forces")
                     dummy => appendChild (atom, forces)
                     ias = idxas (ia, is)
                     force => createElementNS (sclDoc, "", "Hellmann-Fe&
                    &ynman")
                     dummy => appendChild (forces, force)
                     Call setcoord (force, forcehf(:, ias))
                     force => createElementNS (sclDoc, "", "core-correc&
                    &tion")
                     dummy => appendChild (forces, force)
                     Call setcoord (force, forcecr(:, ias))
                     force => createElementNS (sclDoc, "", "IBS")
                     dummy => appendChild (forces, force)
                     Call setcoord (force, forceibs(:, ias))
                     force => createElementNS (sclDoc, "", "totalforce")
                     dummy => appendChild (forces, force)
                     Call setcoord (force, forcetot(:, ias))
                     Write (buffer, '(G22.12)') Sqrt (forcetot(1, &
                    & ias)**2+forcetot(2, ias)**2+forcetot(3, ias)**2)
                     Call setAttribute (forces, "Magnitude", &
                    & trim(adjustl(buffer)))
                  End If
               End Do
               If (input%groundstate%tforce) Then
                  Write (buffer, '(G22.12)') forcemax
                  Call setAttribute (structure, "forceMax", &
                 & trim(adjustl(buffer)))
               End If
            End Do
         End If
      End Subroutine structure_xmlout

      Subroutine scl_iter_xmlout ()
         Implicit None
         Integer :: is, ia, ias
         Real (8) :: scltime, t1
         If (rank .Eq. 0) Then
            niter => createElementNS (sclDoc, "", "iter")
            dummy => appendChild (nscl, niter)
            Write (buffer,*) iscl
            Call setAttribute (niter, "iteration", &
           & trim(adjustl(buffer)))
            Write (buffer,*) currentconvergence
            Call setAttribute (niter, "rms", trim(adjustl(buffer)))
            Write (buffer, '(G22.12)') Log10 (currentconvergence)
            Call setAttribute (niter, "rmslog10", &
           & trim(adjustl(buffer)))
            Write (buffer,*) deltae
            Call setAttribute (niter, "deltae", trim(adjustl(buffer)))
            Write (buffer, '(G22.12)') Log10 (deltae)
            Call setAttribute (niter, "deltaelog10", &
           & trim(adjustl(buffer)))
            if (input%groundstate%tforce) then
               Write (buffer,*) dforcemax
               Call setAttribute (niter, "dforcemax", trim(adjustl(buffer)))
               if (dforcemax .le. 0.d0) then
                 t1=0.d0
               else
                 t1=Log10(dforcemax)
               end if
               Write (buffer, '(G22.12)') t1
               Call setAttribute (niter, "dforcemaxlog10", &
              & trim(adjustl(buffer)))
            end if
            Write (buffer,*) chgdst
            Call setAttribute (niter, "chgdst", trim(adjustl(buffer)))
            Write (buffer, '(G22.12)') Log10 (chgdst)
            Call setAttribute (niter, "chgdstlog10", &
           & trim(adjustl(buffer)))
            Write (buffer, '(G22.12)') fermidos
            Call setAttribute (niter, "fermidos", &
           & trim(adjustl(buffer)))
            Write (buffer, '(G22.12)') efermi
            energies => createElementNS (sclDoc, "", "energies")
            dummy => appendChild (niter, energies)
            Write (buffer, '(G22.12)') engytot
            Call setAttribute (energies, "totalEnergy", &
           & trim(adjustl(buffer)))
            Write (buffer, '(G22.12)') efermi
            Call setAttribute (energies, "fermiEnergy", &
           & trim(adjustl(buffer)))
            Write (buffer, '(G22.12)') evalsum
            Call setAttribute (energies, "sum-of-eigenvalues", &
           & trim(adjustl(buffer)))
            Write (buffer, '(G22.12)') engykn
            Call setAttribute (energies, "electronic-kinetic", &
           & trim(adjustl(buffer)))
            Write (buffer, '(G22.12)') engykncr
            Call setAttribute (energies, "core-electron-kinetic", &
           & trim(adjustl(buffer)))
            Write (buffer, '(G22.12)') engycl
            Call setAttribute (energies, "Coulomb", &
           & trim(adjustl(buffer)))
            Write (buffer, '(G22.12)') engyvcl
            Call setAttribute (energies, "Coulomb-potential", &
           & trim(adjustl(buffer)))
            Write (buffer, '(G22.12)') engynn
            Call setAttribute (energies, "nuclear-nuclear", &
           & trim(adjustl(buffer)))
            Write (buffer, '(G22.12)') engyen
            Call setAttribute (energies, "electron-nuclear", &
           & trim(adjustl(buffer)))
            Write (buffer, '(G22.12)') engyhar
            Call setAttribute (energies, "Hartree", &
           & trim(adjustl(buffer)))
            Write (buffer, '(G22.12)') engymad
            Call setAttribute (energies, "Madelung", &
           & trim(adjustl(buffer)))
            If (input%groundstate%chgexs .Ne. 0.d0) Then
               Write (buffer, '(G22.12)') engycbc
               Call setAttribute (energies, "comp.-background-charge", &
              & trim(adjustl(buffer)))
            End If
            Write (buffer, '(G22.12)') engyvxc
            Call setAttribute (energies, "xc-potential", &
           & trim(adjustl(buffer)))
            If (associated(input%groundstate%spin)) Then
               Write (buffer, '(G22.12)') engybxc
               Call setAttribute (energies, "xc-effective-B-field", &
              & trim(adjustl(buffer)))
               Write (buffer, '(G22.12)') engybext
               Call setAttribute (energies, "external-B-field", &
              & trim(adjustl(buffer)))
            End If
            Write (buffer, '(G22.12)') engyx
            Call setAttribute (energies, "exchange", &
           & trim(adjustl(buffer)))
            Write (buffer, '(G22.12)') engyc
            Call setAttribute (energies, "correlation", &
           & trim(adjustl(buffer)))
            If (ldapu .Ne. 0) Then
               Write (buffer, '(G22.12)') engylu
               Call setAttribute (energies, "LDApU", &
              & trim(adjustl(buffer)))
            End If
            charges => createElementNS (sclDoc, "", "charges")
            dummy => appendChild (niter, charges)
            Write (buffer, '(G18.10)') chgcalc
            Call setAttribute (charges, "totalcharge", &
           & trim(adjustl(buffer)))
            Write (buffer, '(G18.10)') chgcr
            Call setAttribute (charges, "core", trim(adjustl(buffer)))
            Write (buffer, '(G18.10)') chgcrlk
            Call setAttribute (charges, "core_leakage", &
           & trim(adjustl(buffer)))
            Write (buffer, '(G18.10)') chgval
            Call setAttribute (charges, "valence", &
           & trim(adjustl(buffer)))
            Write (buffer, '(G18.10)') chgir
            Call setAttribute (charges, "interstitial", &
           & trim(adjustl(buffer)))
            Write (buffer, '(G18.10)') chgmttot
            Call setAttribute (charges, "muffin-tin-total", &
           & trim(adjustl(buffer)))
            Do is = 1, nspecies
               Do ia = 1, natoms (is)
                  atom => createElementNS (sclDoc, "", "atom")
                  dummy => appendChild (charges, atom)
                  ias = idxas (ia, is)
                  Write (buffer,*) spsymb (is)
                  Call setAttribute (atom, "species", &
                 & trim(adjustl(buffer)))
                  Write (buffer, '(G18.10)') chgmt (ias)
                  Call setAttribute (atom, "muffin-tin", &
                 & trim(adjustl(buffer)))
               End Do
            End Do
            If (input%groundstate%chgexs .Ne. 0.d0) Then
               Write (buffer, '(G18.10)') input%groundstate%chgexs
               Call setAttribute (charges, "excess", &
              & trim(adjustl(buffer)))
            End If
!
            timing => createElementNS (sclDoc, "", "timing")
            dummy => appendChild (niter, timing)
            Call timesec (scltime)
            Write (buffer, '(G22.12)') scltime - scltime0
            scltime0 = scltime
            If (iscl .Ge. 2) Call setAttribute (timing, "itertime", &
           & trim(adjustl(buffer)))
            Write (buffer, '(G22.12)') timeinit + timemat + timefv + &
           & timesv + timerho + timepot + timefor
            Call setAttribute (timing, "timetot", trim(adjustl(buffer)))
            Write (buffer, '(G22.12)') timeinit
            Call setAttribute (timing, "timeinit", trim(adjustl(buffer)))
            Write (buffer, '(G22.12)') timemat
            Call setAttribute (timing, "timemat", trim(adjustl(buffer)))
            Write (buffer, '(G22.12)') timefv
            Call setAttribute (timing, "timefv", trim(adjustl(buffer)))
            Write (buffer, '(G22.12)') timesv
            Call setAttribute (timing, "timesv", trim(adjustl(buffer)))
            Write (buffer, '(G22.12)') timerho
            Call setAttribute (timing, "timerho", trim(adjustl(buffer)))
            Write (buffer, '(G22.12)') timepot
            Call setAttribute (timing, "timepot", trim(adjustl(buffer)))
            Write (buffer, '(G22.12)') timefor
            Call setAttribute (timing, "timefor", trim(adjustl(buffer)))
         End If
      End Subroutine scl_iter_xmlout

      Subroutine setcoord (elementnode, coord)
         Type (Node), Pointer, Intent (In) :: elementnode
         Real (8), Intent (In) :: coord (3)
         Write (buffer, '(G22.12)') coord (1)
         Call setAttribute (elementnode, "x", trim(adjustl(buffer)))
         Write (buffer, '(G22.12)') coord (2)
         Call setAttribute (elementnode, "y", trim(adjustl(buffer)))
         Write (buffer, '(G22.12)') coord (3)
         Call setAttribute (elementnode, "z", trim(adjustl(buffer)))
      End Subroutine setcoord

      Subroutine setcoorddim (elementnode, coord, dim)
         Implicit None
         Type (Node), Pointer, Intent (In) :: elementnode
         Real (8), Intent (In) :: coord (3)
         Integer, Intent (In) :: dim
         If (dim .Gt. 0) Then
            Write (buffer, '(G22.12)') coord (1)
            Call setAttribute (elementnode, "x", trim(adjustl(buffer)))
            If (dim .Gt. 1) Then
               Write (buffer, '(G22.12)') coord (2)
               Call setAttribute (elementnode, "y", &
              & trim(adjustl(buffer)))
               If (dim .Gt. 2) Then
                  Write (buffer, '(G22.12)') coord (3)
                  Call setAttribute (elementnode, "z", &
                 & trim(adjustl(buffer)))
               End If
            End If
         End If
      End Subroutine setcoorddim

      Subroutine scl_xml_write_moments ()
         Type (Node), Pointer :: moments, moment
         Integer :: is, ia, ias
         If (rank .Eq. 0) Then
            moments => createElementNS (sclDoc, "", "moments")
            dummy => appendChild (niter, moments)
!
            moment => createElementNS (sclDoc, "", "momtot")
            dummy => appendChild (moments, moment)
            Call setcoorddim (moment, momtot(1:ndmag), ndmag)
!
            moment => createElementNS (sclDoc, "", "interstitial")
            dummy => appendChild (moments, moment)
            Call setcoorddim (moment, momir(1:ndmag), ndmag)
!
            moment => createElementNS (sclDoc, "", "mommttot")
            dummy => appendChild (moments, moment)
            Call setcoorddim (moment, mommttot(1:ndmag), ndmag)
!
            Do is = 1, nspecies
               Do ia = 1, natoms (is)
                  ias = idxas (ia, is)
                  atom => createElementNS (sclDoc, "", "atom")
                  dummy => appendChild (moments, atom)
                  moment => createElementNS (sclDoc, "", "mommt")
                  dummy => appendChild (atom, moment)
                  Call setcoorddim (moment, mommt(1:ndmag, ias), ndmag)
                  Write (buffer,*) spsymb (is)
                  Call setAttribute (atom, "species", &
                 & trim(adjustl(buffer)))
               End Do
            End Do
         End If
      End Subroutine scl_xml_write_moments

      Subroutine scl_xml_out_close ()!
         If (rank .Eq. 0) Then
            Call destroy (sclDoc)
         End If
      End Subroutine scl_xml_out_close
!
      Subroutine scl_xml_setGndstateStatus (status)
         Character (Len=*) :: status
         If (rank .Eq. 0) Then
            Call setAttribute (ngroundstate, "status", status)
         End If
      End Subroutine scl_xml_setGndstateStatus
!
      Subroutine scl_xml_out_write ()
         If (rank .Eq. 0) Then
            Call normalizeDocument (sclDoc)
            Call serialize (sclDoc, "info" // &
              filext(1:index(filext, ".OUT", .true.)-1) // ".xml")
          End If
      End Subroutine scl_xml_out_write
!
End Module scl_xml_out_Module
