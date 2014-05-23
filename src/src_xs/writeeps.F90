!
!
!
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
Module m_writeeps
      Implicit None
Contains
!
!
      Subroutine writeeps (iq, iop1, iop2, w, eps, fn)
         Use modmain
         Use modxs
         Use FoX_wxml
         Use m_getunit
         Use m_writevars
         Implicit None
    ! arguments
         Integer, Intent (In) :: iq, iop1, iop2
         Real (8), Intent (In) :: w (:)
         Complex (8), Intent (In) :: eps (:)
         Character (*), Intent (In) :: fn
    ! local variables
         Type (xmlf_t), Save :: xf
         Character (256) :: buffer
         Character (*), Parameter :: thisnam = 'writeeps'
         Integer :: n1 (1), n, iw
         Real (8), Allocatable :: imeps (:), kkeps (:)
         If (any(shape(w) .Ne. shape(eps))) Then
            Write (unitout, '(a)') 'Error(' // thisnam // '): input arr&
           &ays have diffenrent shape'
            Call terminate
         End If
         n1 = shape (w)
         n = n1 (1)
         Allocate (imeps(n), kkeps(n))
    ! Kramers-Kronig transform imaginary part
         imeps (:) = aimag (eps(:))
         Call kramkron (iop1, iop2, 1.d-8, n, w, imeps, kkeps)
         Call getunit (unit1)
         Open (unit1, File=trim(fn), Action='write')
         Write (unit1, '(4g18.10)') (w(iw)*escale, eps(iw), kkeps(iw), &
        & iw=1, n)
    ! write relevant parameters to file
!        Call writevars (unit1, iq, iq)
         Close (unit1)
!
    ! write to XML file
         Call xml_OpenFile (trim(fn)//'.xml', xf, replace=.True., &
        & pretty_print=.True.)
         Call xml_NewElement (xf, "dielectric")
         Call xml_DeclareNamespace (xf, "http://www.w3.org/2001/XMLSche&
        &ma-instance", "xsi")
         Call xml_AddAttribute (xf, "xsi:noNamespaceSchemaLocation", ".&
        &./../xml/species.xsd")
         Call xml_NewElement (xf, "mapdef")
         Call xml_NewElement (xf, "variable1")
         Call xml_AddAttribute (xf, "name", "energy")
         Call xml_endElement (xf, "variable1")
         Call xml_NewElement (xf, "function1")
         Call xml_AddAttribute (xf, "name", "macroscopic dielectric fun&
        &ction, real part")
         Call xml_endElement (xf, "function1")
         Call xml_NewElement (xf, "function2")
         Call xml_AddAttribute (xf, "name", "macroscopic dielectric fun&
        &ction, imaginary part")
         Call xml_endElement (xf, "function2")
         Call xml_NewElement (xf, "function3")
         Call xml_AddAttribute (xf, "name", "macroscopic dielectric fun&
        &ction, real part (Kramers Kronig Transform of imaginary part)"&
        & )
         Call xml_endElement (xf, "function3")
         Call xml_endElement (xf, "mapdef")
         Do iw = 1, n
            Call xml_NewElement (xf, "map")
            Write (buffer, '(4g18.10)') w (iw) * escale
            Call xml_AddAttribute (xf, "variable1", trim(adjustl(buffer)))
            Write (buffer, '(4g18.10)') dble (eps(iw))
            Call xml_AddAttribute (xf, "function1", trim(adjustl(buffer)))
            Write (buffer, '(4g18.10)') aimag (eps(iw))
            Call xml_AddAttribute (xf, "function2", trim(adjustl(buffer)))
            Write (buffer, '(4g18.10)') kkeps (iw)
            Call xml_AddAttribute (xf, "function3", trim(adjustl(buffer)))
            Call xml_endElement (xf, "map")
         End Do
         Write (buffer, '(I1.1, ".", I1.1, ".", I3.3)') version
         Call xml_AddComment (xf, " Exciting code version : "//&
        & trim(adjustl(buffer)))
         Call xml_endElement (xf, "dielectric")
         Call xml_Close (xf)
!
         Deallocate (imeps, kkeps)
      End Subroutine writeeps
!
End Module m_writeeps
