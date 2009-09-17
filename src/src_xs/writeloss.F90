


! Copyright (C) 2005-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_writeloss
  implicit none
contains


subroutine writeloss(iq, w, loss, fn)
    use FoX_wxml
    use modmain, only: version
    use mod_lattice
    use mod_constants
    use mod_charge_and_moment
    use modxs
    use m_getunit
    use m_writevars
    implicit none
    ! arguments
    integer, intent(in) :: iq
    real(8), intent(in) :: w(:)
    real(8), intent(in) :: loss(:)
    character(*), intent(in) :: fn
    ! local variables
    character(*), parameter :: thisnam='writeloss'
    type(xmlf_t), save::xf
    character(256)::buffer,buffer2
    integer :: n1(1), n, iw, igmt
    if (any(shape(w).ne.shape(loss))) then
       write(unitout, '(a)') 'Error('//thisnam//'): input arrays have &
	    &diffenrent shape'
       call terminate
    end if
    n1=shape(w); n=n1(1)
    call getunit(unit1)
    open(unit1, file=trim(fn), action='write')
    ! include dynamical structure factor
    ! Dynamical structure factor; expression taken from Weissker, PRL 2006
    ! Units of dynamical structure factor are Hartree^-1
    igmt=ivgigq(ivgmt(1, iq), ivgmt(2, iq), ivgmt(3, iq), iq)
    write(unit1, '(3g18.10)') (w(iw) * escale, loss(iw), loss(iw)* &
	 (gqc(igmt, iq) ** 2/(4.d0 * pi ** 2 * chgval/omega)), iw = 1, n)
    ! write relevant parameters to file
    call writevars(unit1, iq, iq)
    close(unit1)

        ! write to XML file
    call xml_OpenFile (trim(fn)//'.xml', xf, replace=.true.,pretty_print=.true.)
    call xml_NewElement (xf, "loss")
    call xml_DeclareNamespace(xf, "http://www.w3.org/2001/XMLSchema-instance", "xsi")
    call xml_AddAttribute(xf, "xsi:noNamespaceSchemaLocation", "../../xml/species.xsd" )
    call xml_NewElement (xf, "mapdef")
    call xml_NewElement (xf, "variable1")
    call xml_AddAttribute(xf, "name", "energy")
    call xml_endElement (xf, "variable1")
    call xml_NewElement (xf, "function1")
    call xml_AddAttribute(xf, "name", "loss function")
    call xml_endElement (xf, "function1")
    call xml_NewElement (xf, "function2")
    call xml_AddAttribute(xf, "name", "dynamical structure factor")
    call xml_endElement (xf, "function2")
    call xml_endElement (xf, "mapdef")
    do iw=1,n
      call xml_NewElement (xf, "map")
      write(buffer,'(4g18.10)') w(iw)*escale
      call xml_AddAttribute(xf, "variable1", trim(buffer))
      write(buffer,'(4g18.10)') loss(iw)
      call xml_AddAttribute(xf, "function1", trim(buffer))
      write(buffer,'(4g18.10)') loss(iw)* &
     (gqc(igmt, iq) ** 2/(4.d0 * pi ** 2 * chgval/omega))
      call xml_AddAttribute(xf, "function2", trim(buffer))
      call xml_endElement (xf, "map")
    end do
    write(buffer, '(I1.1, ".", I1.1, ".", I3.3)') version
    call xml_AddComment(xf," Exciting code version : "//trim(buffer))
    call xml_endElement (xf, "loss")
    call xml_Close(xf)

  end subroutine writeloss

end module m_writeloss
