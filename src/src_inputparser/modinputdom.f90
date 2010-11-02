
! Copyright (C) 2005-2010 C. Meisenbichler and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

Module inputdom
      Use FoX_dom
!
      Implicit None
      Type (Node), Pointer :: doc, inputnp, nullnode, emptynode, dummy
      Type (DOMConfiguration), Pointer :: config
      Logical :: parseerror
! Request full canonicalization
! ie convert CDATA sections to text sections, remove all entity references etc.
!
!
!
Contains
!
      Subroutine loadinputDOM (deffilename)
      implicit none
      character(*),intent(in)::deffilename
      integer errorcode
      integer iargc
      character(512) filename
         config => newDOMConfig ()
         parseerror = .False.
         Call setParameter (config, "validate-if-schema", .True.)
 !get commandline argument and use as filename
      if (iargc().eq.1)then
          CALL GETARG(1 , filename)
          write(*,*) "### Using specified input file: " // trim(filename)
      else
          filename=deffilename
      endif
       doc => parseFile (trim(filename), config,iostat=errorcode)
         if(errorcode.ne.0) then
        	 write(*,*) "### Could not open ", trim(filename), " file."
      		 write(*,*) "### Check if file exists and if it is valid XML."
       	     stop
         endif
         inputnp => getDocumentElement (doc)
         nullnode => getattributenode (inputnp, "shouldneverexist")
         parseerror = .False.
         dummy => createDocument (getImplementation(), "", "info", &
        & null())
         emptynode => createElementNS (dummy, "", "empty")
      End Subroutine
!
      Subroutine handleunknownnodes (np)
         Type (Node), Pointer :: np, unknownnode
         Type (NodeList), Pointer :: nl
         Type (NamedNodeMap), Pointer :: nnm
         Integer :: Len, i
         Len = 0
         nnm => getAttributes (np)
         Len = getLength (nnm)
         If ( .Not. getNodeName(np) .Eq. "input") Then
            If (len .Gt. 0) Then
               Do i = 0, len - 1
                  unknownnode => item (nnm, i)
                  Write (*,*) ">>>>> unrecognized attribute:"
                  Write (*,*) "      ", getName (unknownnode), '="', &
                 & getValue (unknownnode), '" in element ', getNodeName &
                 & (np)
                  parseerror = .True.
               End Do
            End If
         End If
         nl => getChildNodes (np)
         Len = getLength (nl)
         If (len .Gt. 0) Then
!
            Do i = 0, len - 1
               unknownnode => item (nl, i)
               Select Case (getnodetype(unknownnode))
               Case (ELEMENT_NODE)
                  Write (*,*) ">>>>> unrecognized element:"
                  Write (*,*) "     <", getNodeName (unknownnode), "> i&
                 &n element <", getNodeName (np), ">"
                  parseerror = .True.
               Case (TEXT_NODE)
               End Select
            End Do
         End If
      End Subroutine
!
      Subroutine ifparseerrorstop ()
         If (parseerror) Then
            Write (*,*) "Stopping because of parse error"
            Stop 1
         End If
      End Subroutine
!
!
      Subroutine destroyDOM ()
         Call destroy (doc)
         Call destroy (config)
      End Subroutine
!
End Module
