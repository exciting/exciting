module modplotlabels
  implicit none

  type axisdesc
    character,pointer,dimension(:) ::label=>null(),latexunit=>null(),graceunit=>null()
  end type


   type plotlabels
     type(axisdesc), dimension(:), pointer :: axis
     character, pointer, dimension(:) ::title=>null(),filename=>null()
  end type plotlabels

contains
  pure function vs_str_alloc(s) result(vs)
    character(len=*), intent(in) :: s
    character, dimension(:), pointer :: vs
    allocate(vs(len(s)))

    vs = vs_str(s)
  end function vs_str_alloc

  pure function vs_str(s) result(vs)
    character(len=*), intent(in) :: s
    character, dimension(len(s)) :: vs
    vs = transfer(s, vs)
  end function vs_str
  !BOP
  ! !ROUTINE: create_plotlablels
  ! !INTERFACE:
  function create_plotlablels(title,filename,dimension) result (thisplotlabels)
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
    type(plotlabels),pointer:: thisplotlabels
    character(len=*), intent(in)::title,filename
    integer, intent(in)::dimension
	integer::i
    allocate(thisplotlabels)
    thisplotlabels%title=>   vs_str_alloc(title)
    thisplotlabels%filename=> vs_str_alloc(filename)
    allocate(thisplotlabels%axis(dimension+1))

  end function create_plotlablels
!BOP
  ! !ROUTINE: create_plotlablels
  ! !INTERFACE:
  subroutine set_plotlabel_axis(plotlabels_,axis,label,latexunit,graceunit)
  ! !DESCRIPTION:
    !   This subroutine sets the axis labels for dinmension axis
    !
    ! !INPUT/OUTPUT PARAMETERS:
    !   plotlabels_  :  plotlabels type of which an axis schould be labeled (inout, type(plotlabels))
    !   axis :  number of the axis to set (in,integer)
    !   label: (in,character*)
    !   unit : (in,character*)
    ! !REVISION HISTORY:
    !   Created April 2007 (JKD)
    !EOP
    type(plotlabels),pointer::plotlabels_
    integer, intent(in)::axis
    character(len=*),intent(in)::label,latexunit,graceunit

    plotlabels_%axis(axis)%label=>vs_str_alloc(label)
    plotlabels_%axis(axis)%latexunit=>vs_str_alloc(latexunit)
    plotlabels_%axis(axis)%graceunit=>vs_str_alloc(graceunit)
   end subroutine set_plotlabel_axis

  subroutine destroy_plotlablels(self)
    type(plotlabels),pointer::self
    integer::i
    deallocate(self%title,self%filename)
    do  i=1,size(self%axis)
    deallocate(self%axis(i)%label, self%axis(i)%latexunit,self%axis(i)%graceunit)

    end do
    deallocate(self%axis)
    deallocate(self)
  end subroutine destroy_plotlablels

function get_unit(self,axis) result(unit)
 type(plotlabels)::self
 integer::axis
 character(size(self%axis(axis)%latexunit)+1)::unit
 write(unit,*) self%axis(axis)%latexunit
 end function

 function get_label(self,axis) result(label)
 type(plotlabels)::self
 integer::axis
 character(size(self%axis(axis)%label)+1)::label
 write(label,*) self%axis(axis)%label

 end function

end module modplotlabels
