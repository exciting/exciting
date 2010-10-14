!
!
!
! Copyright (C) 2005-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
Module m_writeloss
      Implicit None
Contains
!
!
      Subroutine writeloss (iq, w, loss, fn)
         Use FoX_wxml
         Use modmain, Only: version
         Use mod_lattice
         Use mod_constants
         Use mod_charge_and_moment
         Use modxs
         Use m_getunit
         Use m_writevars
         Implicit None
    ! arguments
         Integer, Intent (In) :: iq
         Real (8), Intent (In) :: w (:)
         Real (8), Intent (In) :: loss (:)
         Character (*), Intent (In) :: fn
    ! local variables
         Character (*), Parameter :: thisnam = 'writeloss'
         Type (xmlf_t), Save :: xf
         Character (256) :: buffer
         Integer :: n1 (1), n, iw, igmt
         If (any(shape(w) .Ne. shape(loss))) Then
            Write (unitout, '(a)') 'Error(' // thisnam // '): input arr&
           &ays have diffenrent shape'
            Call terminate
         End If
         n1 = shape (w)
         n = n1 (1)
         Call getunit (unit1)
         Open (unit1, File=trim(fn), Action='write')
    ! include dynamical structure factor
    ! Dynamical structure factor; expression taken from Weissker, PRL 2006
    ! Units of dynamical structure factor are Hartree^-1
         igmt = ivgigq (ivgmt(1, iq), ivgmt(2, iq), ivgmt(3, iq), iq)
         Write (unit1, '(3g18.10)') (w(iw)*escale, loss(iw), &
        & loss(iw)*(gqc(igmt, iq)**2/(4.d0*pi**2*chgval/omega)), iw=1, &
        & n)
    ! write relevant parameters to file
         Call writevars (unit1, iq, iq)
         Close (unit1)
!
        ! write to XML file
         Call xml_OpenFile (trim(fn)//'.xml', xf, replace=.True., &
        & pretty_print=.True.)
         Call xml_NewElement (xf, "loss")
         Call xml_DeclareNamespace (xf, "http://www.w3.org/2001/XMLSche&
        &ma-instance", "xsi")
         Call xml_AddAttribute (xf, "xsi:noNamespaceSchemaLocation", ".&
        &./../xml/species.xsd")
         Call xml_NewElement (xf, "mapdef")
         Call xml_NewElement (xf, "variable1")
         Call xml_AddAttribute (xf, "name", "energy")
         Call xml_endElement (xf, "variable1")
         Call xml_NewElement (xf, "function1")
         Call xml_AddAttribute (xf, "name", "loss function")
         Call xml_endElement (xf, "function1")
         Call xml_NewElement (xf, "function2")
         Call xml_AddAttribute (xf, "name", "dynamical structure factor&
        &")
         Call xml_endElement (xf, "function2")
         Call xml_endElement (xf, "mapdef")
         Do iw = 1, n
            Call xml_NewElement (xf, "map")
            Write (buffer, '(4g18.10)') w (iw) * escale
            Call xml_AddAttribute (xf, "variable1", trim(adjustl(buffer)))
            Write (buffer, '(4g18.10)') loss (iw)
            Call xml_AddAttribute (xf, "function1", trim(adjustl(buffer)))
            Write (buffer, '(4g18.10)') loss (iw) * (gqc(igmt, &
           & iq)**2/(4.d0*pi**2*chgval/omega))
            Call xml_AddAttribute (xf, "function2", trim(adjustl(buffer)))
            Call xml_endElement (xf, "map")
         End Do
         Write (buffer, '(I1.1, ".", I1.1, ".", I3.3)') version
         Call xml_AddComment (xf, " Exciting code version : "//&
        & trim(adjustl(buffer)))
         Call xml_endElement (xf, "loss")
         Call xml_Close (xf)
!
      End Subroutine writeloss
!
End Module m_writeloss
