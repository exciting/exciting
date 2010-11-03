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
Subroutine plot1d (fname, nf, lmax, ld, rfmt, rfir, plotdef)
! !USES:
      Use modinput
      Use modmain
      use modmpi
      Use FoX_wxml
! !INPUT/OUTPUT PARAMETERS:
!   fnum1 : plot file name (character*)
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
      Character (Len=*), Intent (In) :: fname
      Integer, Intent (In) :: nf
      Integer, Intent (In) :: lmax
      Integer, Intent (In) :: ld
      Real (8), Intent (In) :: rfmt (ld, nrmtmax, natmtot, nf)
      Real (8), Intent (In) :: rfir (ngrtot, nf)
      Type (plot1d_type) :: plotdef

! local variables
      Integer :: i, ip, iv, fnum1 = 50, fnum2 = 51
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
      buffer = fname // "1D.OUT"
      Open (fnum1, File=trim(buffer), Action='WRITE', Form='FORMATTED')
      buffer = fname // "LINES.OUT"
      Open (fnum2, File=trim(buffer), Action='WRITE', Form='FORMATTED')
      Call xml_OpenFile (fname//"1d.xml", xf, replace=.True., &
     & pretty_print=.True.)
      Call xml_NewElement (xf, "plot1d")
      Call xml_NewElement (xf, "title")
      Call xml_AddCharacters (xf, trim(input%title))
      Call xml_endElement (xf, "title")
! connect the plotting vertices
      nvp1d = size (plotdef%path%pointarray)
      npp1d = plotdef%path%steps
      Allocate (fp(npp1d, nf))
      If (allocated(dvp1d)) deallocate (dvp1d)
      Allocate (dvp1d(nvp1d))
      If (allocated(vplp1d)) deallocate (vplp1d)
      Allocate (vplp1d(3, npp1d))
      If (allocated(dpp1d)) deallocate (dpp1d)
      Allocate (dpp1d(npp1d))
      Call connect (input%structure%crystal%basevect, plotdef, &
     & size(plotdef%path%pointarray), plotdef%path%steps, vplp1d, &
     & dvp1d, dpp1d)
      Do i = 1, nf
! evaluate function at each point
         Call rfarray (lmax, ld, rfmt(:, :, :, i), rfir(:, i), npp1d, &
        & vplp1d, fp(:, i))
      End Do
      fmin = fp (1, 1)
      fmax = fp (1, 1)
      Do ip = 1, npp1d
         Do i = 1, nf
            fmin = Min (fmin, fp(ip, i))
            fmax = Max (fmax, fp(ip, i))
         End Do
! write the point distances and function to file
         Write (fnum1, '(5G18.10)') dpp1d (ip), (fp(ip, i), i=1, nf)
         Call xml_NewElement (xf, "point")
         Write (buffer, '(5G18.10)') dpp1d (ip)
         Call xml_AddAttribute (xf, "distance", trim(adjustl(buffer)))
         Do i = 1, nf
            Write (buffer, '(5G18.10)') fp (ip, i)
            Write (buffer1,*) i
            Call xml_AddAttribute (xf, "function"//&
           & trim(adjustl(buffer1)), trim(adjustl(buffer)))
         End Do
         Call xml_endElement (xf, "point")
      End Do
! write the vertex location lines
      t1 = 0.5d0 * (fmax-fmin)
      Do iv = 1, nvp1d
         Write (fnum2, '(2G18.10)') dvp1d (iv), fmax + t1
         Write (fnum2, '(2G18.10)') dvp1d (iv), fmin - t1
         Write (fnum2, '("     ")')
      End Do
      Do iv = 1, size (plotdef%path%pointarray)
         Call xml_NewElement (xf, "vertex")
!
         Write (buffer, '(5G18.10)') dvp1d (iv)
         Call xml_AddAttribute (xf, "distance", trim(adjustl(buffer)))
         Write (buffer, '(5G18.10)') fmax + t1
         Call xml_AddAttribute (xf, "upperboundary", &
        & trim(adjustl(buffer)))
         Write (buffer, '(5G18.10)') fmin - t1
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
      Close (fnum1)
      Close (fnum2)
      Deallocate (fp)
      Deallocate (dvp1d)
      Deallocate (vplp1d)
      Deallocate (dpp1d)
      endif
      Return
End Subroutine
!EOC
