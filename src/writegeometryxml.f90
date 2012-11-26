
! Copyright (C) 2005-2010 C. Meisenbichler and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: writegeometryxml
! !INTERFACE:
!
Subroutine writegeometryxml (topt)
! !USES:
      Use modmain
      Use modinput
      Use modspdeflist
      Use modspdb
      Use FoX_wxml
! !INPUT/OUTPUT PARAMETERS:
!   topt : if .true. then the filename will be {\ttgeometry_opt.xml}, otherwise
!          {\tt geometry.xml} (in,logical)
! !DESCRIPTION:
!   Outputs the lattice vectors and atomic positions to file, in a format
!   which may be then used directly in {\tt exciting.in}.
!
! !REVISION HISTORY:
!   Created January 2004 (JKD)
!EOP
!BOC
      Implicit None
      Logical, Intent (In) :: topt
! local variables
      Integer :: is, ia, i
      Character (128) :: buffer
      Type (xmlf_t), Save :: xf
      Real (8) :: v (3)
      If (topt) Then
         Call xml_OpenFile ("geometry_opt.xml", xf, &
        & replace=.True., pretty_print=.True.)
      End If
      Call xml_NewElement (xf, "input")
      Call xml_NewElement (xf, "structure")
      Call xml_AddAttribute (xf, "speciespath", &
     & trim(adjustl(input%structure%speciespath)))
! GB 6.11.2012
            If (input%structure%cartesian)Then
               Write(buffer,*) "true"
               Call xml_AddAttribute (xf,"cartesian", &
                    & trim(adjustl(buffer)))
            End If
      Call xml_NewElement (xf, "crystal")
      Do i = 1, 3
         Call xml_NewElement (xf, "basevect")
         Write (buffer, '(3G18.10)') &
        & input%structure%crystal%basevect(:, i)
         Call xml_AddCharacters (xf, trim(buffer))
         Call xml_endElement (xf, "basevect")
      End Do
      Call xml_endElement (xf, "crystal")
      Do is = 1, nspecies
         Call xml_NewElement (xf, "species")
         Call xml_AddAttribute (xf, "speciesfile", trim(adjustl(input%structure%speciesarray(is)%species%speciesfile)))
         Write (buffer,*) Int (-1.0*speziesdeflist(is)%sp%z)
         Do ia = 1, natoms (is)
            Call xml_NewElement (xf, "atom")
! GB 04.10.2012 adding the CARTESIAN switch, so that geometry_opt.xml will be consistent with Cartesian input 
            If (input%structure%cartesian) Then  !GB
! write Cartesian coordinates for the molecular case   
               Call r3mv (input%structure%crystal%basevect, input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:), &
              & v)
            Else
! otherwise write lattice coordinates  
               v (:) = input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:)
            End If
! END CHANGES
!            Write (buffer, '(3G18.10)')  input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord
            Write (buffer, '(3G18.10)') (v (:)) 
           Call xml_AddAttribute (xf, "coord", trim(adjustl(buffer)))
            If (input%structure%speciesarray(is)%species%atomarray(ia)%atom%lock) then
            Write (buffer, *)  "true"
            Call xml_AddAttribute (xf, "lock", &                                                                                                                                                                                                                      
                 & trim(adjustl(buffer)))
            End If
            If (associated(input%groundstate%spin)) Then
            Write (buffer, *)  input%structure%speciesarray(is)%species%atomarray(ia)%atom%bfcmt(:)
            Call xml_AddAttribute (xf, "bfcmt", &
                 & trim(adjustl(buffer)))
            End If

            Call xml_endElement (xf, "atom")
         End Do
         Call xml_endElement (xf, "species")
      End Do
      Call xml_close (xf)
      Return
End Subroutine
