!
!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: plot1d
! !INTERFACE:
!
!
Subroutine plot1d (labels, nf, lmax, ld, rfmt, rfir, plotdef)
  ! !USES:
  Use modinput
  use modmpi
  Use FoX_wxml
  use mod_muffin_tin
  use mod_atoms
  use mod_Gvector
  use modplotlabels
  use mod_plotting
  ! !INPUT/OUTPUT PARAMETERS:
  !   lables : plot labels (character*)
  !   nf    : number of functions (in,integer)
  !   lmax  : maximum angular momentum (in,integer)
  !   ld    : leading dimension (in,integer)
  !   rfmt  : real muffin-tin function (in,real(ld,nrmtmax,natmtot,nf))
  !   rfir  : real intersitial function (in,real(ngrtot,nf))
  !   plotdef: type(plot1d) defining plotregion
  ! !DESCRIPTION:
  !   Produces a 1D plot of the real functions contained in arrays {\tt rfmt} and
  !   {\tt rfir} along the lines connecting the vertices in the global array
  !   {\tt vvlp1d}. See routine {\tt rfarray}.
  !
  ! !REVISION HISTORY:
  !   Created June 2003 (JKD)
  !EOP
  !BOC
  Implicit None
  ! arguments
  type(plotlabels), Intent (In) :: labels
  Integer, Intent (In) :: nf
  Integer, Intent (In) :: lmax
  Integer, Intent (In) :: ld
  Real (8), Intent (In) :: rfmt (ld, nrmtmax, natmtot, nf)
  Real (8), Intent (In) :: rfir (ngrtot, nf)
  Type (plot1d_type) :: plotdef

  ! local variables
  Integer :: i, ip, iv
  Real (8) :: fmin, fmax, t1
  Character (128) :: buffer, buffer1
  Type (xmlf_t), Save :: xf
  ! allocatable arrays
  Real (8), Allocatable :: fp (:, :)
  If (rank .Eq. 0) Then
     If ((nf .Lt. 1) .Or. (nf .Gt. 4)) Then
        Write (*,*)
        Write (*, '("Error(plot1d): invalid number of functions : ", I&
             &8)') nf
        Write (*,*)
        Stop
     End If
     write(buffer,*)  labels%filename , ".xml"
     Call xml_OpenFile ( adjustl(trim(buffer)), xf, replace=.True. ,  pretty_print=.True.)
     Call xml_NewElement (xf, "plot1d")
     Call xml_NewElement (xf, "title")
     Call xml_AddCharacters (xf, trim(input%title))
     Call xml_endElement (xf, "title")
     Call xml_NewElement (xf, "grid")
     Call xml_NewElement (xf, "axis")
     Call xml_AddAttribute (xf, "label", get_label(labels,1))
     Call xml_AddAttribute (xf, "latexunit", get_latexunit(labels,1))
     Call xml_AddAttribute (xf, "graceunit", get_graceunit(labels,1))
     Call xml_endElement (xf, "axis")
     Call xml_NewElement (xf, "axis")
     Call xml_AddAttribute (xf, "label", get_label(labels,2))
     Call xml_AddAttribute (xf, "latexunit", get_latexunit(labels,2))
     Call xml_AddAttribute (xf, "graceunit", get_graceunit(labels,2))
     Call xml_endElement (xf, "axis")
     ! connect the plotting vertices
     nvp1d = size (plotdef%path%pointarray)
     npp1d = plotdef%path%steps
     If (allocated(dvp1d)) deallocate (dvp1d)
     Allocate (dvp1d(nvp1d))
     If (allocated(vplp1d)) deallocate (vplp1d)
     Allocate (vplp1d(3, npp1d))
     If (allocated(dpp1d)) deallocate (dpp1d)
     Allocate (dpp1d(npp1d))
     Call connect(input%structure%crystal%basevect, plotdef, &
     &            size(plotdef%path%pointarray), plotdef%path%steps, &
     &            vplp1d, dvp1d, dpp1d)
     ! evaluate function at each point 
     Allocate (fp(npp1d, nf))
     fp = 0.d0
     Do i = 1, nf
        Call rfarray(lmax, ld, rfmt(:, :, :, i), rfir(:, i), npp1d, vplp1d, fp(:, i))
     End Do
     ! write the point distances and function to file
     do i = 1, nf
        call xml_NewElement(xf, "function")
        call xml_AddAttribute(xf, "name", "")
        Do ip = 1, npp1d
          Call xml_NewElement (xf, "point")
          Write (buffer, '(G18.10)') dpp1d (ip)
          Call xml_AddAttribute (xf, "distance", trim(adjustl(buffer)))
          Write (buffer, '(G18.10)') fp (ip, i)
          Call xml_AddAttribute (xf, "value", trim(adjustl(buffer)))
          Call xml_endElement (xf, "point")
        End Do ! ip
        Call xml_endElement (xf, "function")
     End Do ! nf
     ! write the vertex location lines
     fmin = minval(fp)
     fmax = maxval(fp)    
     Do iv = 1, size (plotdef%path%pointarray)
        Call xml_NewElement (xf, "vertex")
        Write (buffer, '(5G18.10)') dvp1d (iv)
        Call xml_AddAttribute (xf, "distance", trim(adjustl(buffer)))
        Write (buffer, '(5G18.10)') fmax
        Call xml_AddAttribute (xf, "upperboundary", &
             & trim(adjustl(buffer)))
        Write (buffer, '(5G18.10)') fmin
        Call xml_AddAttribute (xf, "lowerboundary", &
             & trim(adjustl(buffer)))
        Call xml_AddAttribute (xf, "label", &
             & trim(adjustl(plotdef%path%pointarray(iv)%point%label)))
        Write (buffer, '(5G18.10)') &
             & plotdef%path%pointarray(iv)%point%coord
        Call xml_AddAttribute (xf, "coord", trim(adjustl(buffer)))
        Call xml_endElement (xf, "vertex")
     End Do
     Call xml_close (xf)
     
     Deallocate (fp)
     Deallocate (dvp1d)
     Deallocate (vplp1d)
     Deallocate (dpp1d)
  endif
  Return
End Subroutine plot1d
!EOC
