


! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_writeeps
  implicit none
contains


subroutine writeeps(iq, iop1, iop2, w, eps, fn)
    use modmain
    use modxs
    use FoX_wxml
    use m_getunit
    use m_writevars
    implicit none
    ! arguments
    integer, intent(in) :: iq, iop1, iop2
    real(8), intent(in) :: w(:)
    complex(8), intent(in) :: eps(:)
    character(*), intent(in) :: fn
    ! local variables
    type(xmlf_t), save::xf
    character(256)::buffer,buffer2
    character(*), parameter :: thisnam='writeeps'
    integer :: n1(1), n, iw
    real(8), allocatable :: imeps(:), kkeps(:)
    if (any(shape(w).ne.shape(eps))) then
       write(unitout, '(a)') 'Error('//thisnam//'): input arrays have &
	    &diffenrent shape'
       call terminate
    end if
    n1=shape(w); n=n1(1)
    allocate(imeps(n), kkeps(n))
    ! Kramers-Kronig transform imaginary part
    imeps(:)=aimag(eps(:))
    call kramkron(iop1, iop2, 1.d-8, n, w, imeps, kkeps)
    call getunit(unit1)
    open(unit1, file=trim(fn), action='write')
    write(unit1, '(4g18.10)') (w(iw)*escale, eps(iw), kkeps(iw), iw=1, n)
    ! write relevant parameters to file
    call writevars(unit1, iq, iq)
    close(unit1)

    ! write to XML file
    call xml_OpenFile (trim(fn)//'.xml', xf, replace=.true.,pretty_print=.true.)
    call xml_NewElement (xf, "dielectric")
    call xml_DeclareNamespace(xf, "http://www.w3.org/2001/XMLSchema-instance", "xsi")
    call xml_AddAttribute(xf, "xsi:noNamespaceSchemaLocation", "../../xml/species.xsd" )
    call xml_NewElement (xf, "mapdef")
    call xml_NewElement (xf, "variable1")
    call xml_AddAttribute(xf, "name", "energy")
    call xml_endElement (xf, "variable1")
    call xml_NewElement (xf, "function1")
    call xml_AddAttribute(xf, "name", "macroscopic dielectric function, real part")
    call xml_endElement (xf, "function1")
    call xml_NewElement (xf, "function2")
    call xml_AddAttribute(xf, "name", "macroscopic dielectric function, imaginary part")
    call xml_endElement (xf, "function2")
    call xml_NewElement (xf, "function3")
    call xml_AddAttribute(xf, "name", "macroscopic dielectric function, real part &
    &(Kramers Kronig Transform of imaginary part)")
    call xml_endElement (xf, "function3")
    call xml_endElement (xf, "mapdef")
    do iw=1,n
      call xml_NewElement (xf, "map")
      write(buffer,'(4g18.10)') w(iw)*escale
      call xml_AddAttribute(xf, "variable1", trim(buffer))
      write(buffer,'(4g18.10)') dble(eps(iw))
      call xml_AddAttribute(xf, "function1", trim(buffer))
      write(buffer,'(4g18.10)') aimag(eps(iw))
      call xml_AddAttribute(xf, "function2", trim(buffer))
      write(buffer,'(4g18.10)') kkeps(iw)
      call xml_AddAttribute(xf, "function3", trim(buffer))
      call xml_endElement (xf, "map")
    end do
    write(buffer, '(I1.1, ".", I1.1, ".", I3.3)') version
    call xml_AddComment(xf," Exciting code version : "//trim(buffer))
    call xml_endElement (xf, "dielectric")
    call xml_Close(xf)

    deallocate(imeps, kkeps)
  end subroutine writeeps

end module m_writeeps
