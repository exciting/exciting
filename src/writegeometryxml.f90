
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
!   which may be then used directly in {\tt input.xml}.
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
      Logical :: lock(3)

      if (topt) then
          call xml_OpenFile ("geometry_opt.xml", xf, replace=.True., pretty_print=.True.)
      else
          call xml_OpenFile ("geometry"//filext(1:index(filext, ".OUT", .true.)-1)//".xml", &
         &     xf, replace=.True., pretty_print=.True.)
      end if

      call xml_NewElement (xf, "input")

      ! ----- Title Element -----
      call xml_NewElement(xf,"title")
      call xml_AddCharacters(xf,trim(adjustl(input%title)))
      if (topt) call xml_AddCharacters(xf," (optimized)")
      call xml_EndElement(xf,"title")
      
      call xml_NewElement (xf, "structure")
      call xml_AddAttribute (xf, "speciespath", trim(adjustl(input%structure%speciespath)))
      if (input%structure%cartesian) call xml_AddAttribute (xf,"cartesian", trim(adjustl("true")))
      call xml_NewElement (xf, "crystal")

      do i = 1, 3
          call xml_NewElement (xf, "basevect")
          write (buffer, '(3G18.10)') input%structure%crystal%basevect(:, i)
          call xml_AddCharacters (xf, trim(buffer))
          call xml_endElement (xf, "basevect")
      end do

      call xml_endElement (xf, "crystal")

      do is = 1, nspecies
          call xml_NewElement (xf, "species")
          call xml_AddAttribute (xf, "speciesfile", &
         &  trim(adjustl(input%structure%speciesarray(is)%species%speciesfile)))
          write(buffer,'(F7.4)') rmt(is)
          call xml_AddAttribute (xf, "rmt", trim(adjustl(buffer)))

          do ia = 1, natoms (is)
              call xml_NewElement (xf, "atom")

              if (input%structure%cartesian) then   
                  call r3mv (input%structure%crystal%basevect, &
                 &  input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:), v)
              else
                  v(:) = input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:)
              end if
              write (buffer, '(3F18.10)') (v (:)) 
              call xml_AddAttribute (xf, "coord", trim(adjustl(buffer)))            

              lock(:) = .False.
              if ( associated(input%relax) ) then
                  lock(:) = input%structure%speciesarray(is)%species%atomarray(ia)%atom%lockxyz(:)
              end if

              if (lock(1).or.lock(2).or.lock(3)) then
                  write(buffer,*) printLogical(lock(1)), &
                 &                printLogical(lock(2)), &
                 &                printLogical(lock(3))
                  call xml_AddAttribute (xf, "lockxyz", trim(adjustl(buffer)))
              end if

              if (associated(input%groundstate%spin)) then
                  write (buffer, *) input%structure%speciesarray(is)%species%atomarray(ia)%atom%bfcmt(:)
                  call xml_AddAttribute (xf, "bfcmt", trim(adjustl(buffer)))
              end if

              call xml_endElement (xf, "atom")
          end do

          call xml_endElement (xf, "species")

      end do

      call xml_close (xf)
      return

contains      

      function printLogical(flag)
          logical, intent(IN) :: flag
          character(6) :: printLogical
          printLogical="false"
          if (flag) write(printLogical,'("true")')
      end function

End Subroutine
