
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
      Use modsp
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
      If (topt) Then
         Call xml_OpenFile ("geometry_opt"//trim(filext)//".xml", xf, &
        & replace=.True., pretty_print=.True.)
      Else
         Call xml_OpenFile ("geometry"//trim(filext)//".xml", xf, &
        & replace=.True., pretty_print=.True.)
      End If
      Call xml_AddXMLPI(xf,"xml-stylesheet", 'href="'//trim(input%xsltpath)//&
      &'/inputfileconverter/inputtohtml.xsl" type="text/xsl"')
      Call xml_NewElement (xf, "input")
      Call xml_NewElement (xf, "structure")
      If (input%structure%primcell) buffer = "true"
      If ( .Not. input%structure%primcell) buffer = "false"
      Call xml_AddAttribute (xf, "primcell", trim(adjustl(buffer)))
      Call xml_AddAttribute (xf, "speciespath", &
     & trim(adjustl(input%structure%speciespath)))
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
         Call xml_AddAttribute (xf, "atomicNumber", &
        & trim(adjustl(buffer)))
         Call xml_AddAttribute (xf, "chemicalSymbol", &
        & trim(adjustl(speziesdeflist(is)%sp%chemicalSymbol)))
         Do ia = 1, natoms (is)
            Call xml_NewElement (xf, "atom")
            Write (buffer, '(3G18.10)') input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord
            Call xml_AddAttribute (xf, "coord", trim(adjustl(buffer)))
            Call xml_endElement (xf, "atom")
         End Do
         Call xml_endElement (xf, "species")
      End Do
      Call xml_close (xf)
      Return
End Subroutine
