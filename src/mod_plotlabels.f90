module modplotlabels
    implicit none
    type axisdesc
	character,pointer,dimension(:) ::label=>null(),unit=>null()

end type

type plotlabels
	    character, pointer, dimension(:) ::title=>null(),filename=>null()
	type(axisdesc),allocatable :: axis(:)
end type

contains
  pure function vs_str_alloc(s) result(vs)
    character(len=*), intent(in) :: s
    character, dimension(:), pointer :: vs
    allocate(vs(len(s)))
    vs = s
  end function vs_str_alloc

!BOP
! !ROUTINE: create_plotlablels
! !INTERFACE:
function create_plotlablels(title,filename,dimension) result (plotlabels)
! !DESCRIPTION:
!   This function creates a plotlabels type and initializes the axis descriptions
!
! !INPUT/OUTPUT PARAMETERS:
!   title     :  string for plot title (in, character(512))
!   dimension :  number of dimensions. e.g 1d plot has one dimension but 2 axes (in, integer)
!   create_plotlablels : pointer to plotlabels type  (out, type(plotlabels))
! !REVISION HISTORY:
!   Created April 2007 (JKD)
!EOP
	type(plotlabels),pointer:: plotlabels
	character(len=*), intent(in)::title,filename
	integer, intent(in)::dimension

	allocate(plotlabels)
    plotlabels%title=>   vs_str_alloc(title)
    plotlabels%filename=> vs_str_alloc(filename)
	allocate(plotlabels%axis(dimension+1))

end function

subroutine destroy_plotlablels(self)
	type(plotlabels),pointer::self

	deallocate(self%axis,self%title,self%filename)
	deallocate(self)
end subroutine



end module
