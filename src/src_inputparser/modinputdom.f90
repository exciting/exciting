! Copyright (C) 2005-2010 C. Meisenbichler and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

Module inputdom
      Use FoX_dom
! Module for handling the issues connected to parsing
! the input file and construcing the input object.
      Implicit None
      Type (Node), Pointer :: doc, inputnp, nullnode, emptynode, dummy
      Type (DOMConfiguration), Pointer :: config
      Logical :: parseerror
! 
!
!
Contains

!BOP
! !ROUTINE: loadinputDOM
! !INTERFACE:
      Subroutine loadinputDOM (deffilename)
! !INPUT/OUTPUT PARAMETERS:
! deffilename: default filename if no file name argument is in  argv()
      implicit none

      character(*),intent(in)::deffilename
! !DESCRIPTION:
! Loads the contents of the inputfile into the DOM and initializes some data
! that is used by getstructinput(inputnp)

!EOP
!BOC

      integer errorcode
      integer iargc
      character(512) filename

         config => newDOMConfig ()
         parseerror = .False.
         Call setParameter (config, "validate-if-schema", .True.)
 ! get commandline argument and use as filename
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
    ! initialise some fields of this module
    ! this saves the pointer to the root of the DOM into the inputnp
         inputnp => getDocumentElement (doc)
    ! the nullnode is an version of the node type that says it doesnt exist.
    ! used in other functions to be able to return an empy node
         nullnode => getattributenode (inputnp, "shouldneverexist")
         parseerror = .False.
         dummy => createDocument (getImplementation(), "", "info", &
        & null())
    ! This is a node without content. It is also kept around to be asigned.
         emptynode => createElementNS (dummy, "", "empty")
      End Subroutine
!BOP
! !ROUTINE: handleunknownnodes
! !INTERFACE:
      Subroutine handleunknownnodes (np)
! !INPUT/OUTPUT PARAMETERS:
! np: node pointer type containing node with all
!     known nodes removed, thus containing only unknown that need to be reported.
! !DESCRIPTION:
! This writes the error message when the getstruct... function sees an unknown (illegal) entry
 
!EOP
!BOC
         Type (Node), Pointer :: np, unknownnode
         Type (NodeList), Pointer :: nl
         Type (NamedNodeMap), Pointer :: nnm
         Integer :: Len, i
         Len = 0
         nnm => getAttributes (np)
         Len = getLength (nnm)
         If ( .Not. getNodeName(np) .Eq. "input") Then
         ! if there are any attributes in the node
            If (len .Gt. 0) Then
            !compose error message
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
         !if there are any child element in the node
         If (len .Gt. 0) Then
			! compose error message
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
