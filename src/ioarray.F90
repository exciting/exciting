
! Copyright (C) 2007 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module ioarray
  implicit none
  integer, parameter :: kndi=4
  integer, parameter :: kndr=8
  integer, parameter :: kndc=8
  integer, parameter :: fmtlen=1024
  integer, parameter :: und=1
  logical, parameter :: tparend=.true.
  character(*), parameter :: fmtdl='l2,1x'
  character(*), parameter :: fmtdi='i12'
  character(*), parameter :: fmtdr='es23.15E3'
  character(*), parameter :: fmtdc='es23.15E3'
  private :: kndi,kndr,kndc,fmtlen,und
  private :: fmtdl,fmtdi,fmtdr,fmtdc
  private :: ioaction
contains
  subroutine ioarr( &
       arr1dl,arr1di,arr1dr,arr1dc,  &
       arr2dl,arr2di,arr2dr,arr2dc,  &
       arr3dl,arr3di,arr3dr,arr3dc,  &
       arr4dl,arr4di,arr4dr,arr4dc,  &
       arr5dl,arr5di,arr5dr,arr5dc,  &
       arr6dl,arr6di,arr6dr,arr6dc,  &
       arr7dl,arr7di,arr7dr,arr7dc,  &
       ioa,un,fmt,tparen)
    implicit none
    ! arguments
    character(*), intent(in) :: ioa
    ! * 1D arrays *
    logical,       optional :: arr1dl(:)
    integer(kndi), optional :: arr1di(:)
    real(kndr),    optional :: arr1dr(:)
    complex(kndc), optional :: arr1dc(:)
    ! * 2D arrays *
    logical,       optional :: arr2dl(:,:)
    integer(kndi), optional :: arr2di(:,:)
    real(kndr),    optional :: arr2dr(:,:)
    complex(kndc), optional :: arr2dc(:,:)
    ! * 3D arrays *
    logical,       optional :: arr3dl(:,:,:)
    integer(kndi), optional :: arr3di(:,:,:)
    real(kndr),    optional :: arr3dr(:,:,:)
    complex(kndc), optional :: arr3dc(:,:,:)
    ! * 4D arrays *
    logical,       optional :: arr4dl(:,:,:,:)
    integer(kndi), optional :: arr4di(:,:,:,:)
    real(kndr),    optional :: arr4dr(:,:,:,:)
    complex(kndc), optional :: arr4dc(:,:,:,:)
    ! * 5D arrays *
    logical,       optional :: arr5dl(:,:,:,:,:)
    integer(kndi), optional :: arr5di(:,:,:,:,:)
    real(kndr),    optional :: arr5dr(:,:,:,:,:)
    complex(kndc), optional :: arr5dc(:,:,:,:,:)
    ! * 6D arrays *
    logical,       optional :: arr6dl(:,:,:,:,:,:)
    integer(kndi), optional :: arr6di(:,:,:,:,:,:)
    real(kndr),    optional :: arr6dr(:,:,:,:,:,:)
    complex(kndc), optional :: arr6dc(:,:,:,:,:,:)
    ! * 7D arrays *
    logical,       optional :: arr7dl(:,:,:,:,:,:,:)
    integer(kndi), optional :: arr7di(:,:,:,:,:,:,:)
    real(kndr),    optional :: arr7dr(:,:,:,:,:,:,:)
    complex(kndc), optional :: arr7dc(:,:,:,:,:,:,:)
    ! * other arguments *
    integer, optional, intent(in) :: un
    character(*), optional, intent(in) :: fmt
    logical, optional, intent(in) :: tparen
    ! local variables
    integer :: i1,i2,i3,i4,i5,i6,i7,npr,unt,ndim
    integer :: j1,j2,j3,j4,j5,j6,j7
    integer :: sh(7),lb(7),ub(7)
    character(1) :: iot
    character(fmtlen) :: frmt,str
    logical :: tpr(28),twrite,tun,tfmt,tparent
    logical :: lt
    integer :: it
    real(kndr) :: rt
    complex(kndc) :: zt
    logical, external :: ioaction
    character(1), external :: iotype
    ! check if only one array is present
    tpr=(/ &
         present(arr1dl),present(arr1di),present(arr1dr),present(arr1dc), &
         present(arr2dl),present(arr2di),present(arr2dr),present(arr2dc), &
         present(arr3dl),present(arr3di),present(arr3dr),present(arr3dc), &
         present(arr4dl),present(arr4di),present(arr4dr),present(arr4dc), &
         present(arr5dl),present(arr5di),present(arr5dr),present(arr5dc), &
         present(arr6dl),present(arr6di),present(arr6dr),present(arr6dc), &
         present(arr7dl),present(arr7di),present(arr7dr),present(arr7dc) &
         /)
    npr=count(tpr)
    if (npr.ne.1) then
       write(*,*)
       write(*,*) 'Error(ioarray::ioarr): more than one array specified'
       write(*,*)
       stop
    end if
    ! determine dimension
    if (any(tpr(1:4))) ndim=1
    if (any(tpr(5:8))) ndim=2
    if (any(tpr(9:12))) ndim=3
    if (any(tpr(13:16))) ndim=4
    if (any(tpr(17:20))) ndim=5
    if (any(tpr(21:24))) ndim=6
    if (any(tpr(25:28))) ndim=7
    ! determine type
    if (any(tpr(1::4))) iot='l'
    if (any(tpr(2::4))) iot='i'
    if (any(tpr(3::4))) iot='r'
    if (any(tpr(4::4))) iot='c'
    ! determine IO action (read or write)
    twrite=ioaction(ioa)
    ! optional parentheses for complex numbers
    tparent=tparend
    if (present(tparen)) tparent=tparen
    if (iot.ne.'c') tparent=.false.
    ! determine file unit
    unt=und
    if (present(un)) unt=un
    ! determine format
    select case (iot)
    case('l')
       frmt=trim(fmtdl)
    case('i')
       frmt=trim(fmtdi)
    case('r')
       frmt=trim(fmtdr)
    case('c')
       frmt=trim(fmtdc)
    end select
    if (present(fmt)) frmt=trim(adjustl(fmt))
    !*********************************************************************
    write(*,*) 'npr,iot,ndim',npr,iot,ndim
    !*********************************************************************
    ! select by dimension
    select case(ndim)
    case(1)
       !---------------------!
       !     1 dimension     !
       !---------------------!
       select case(iot)
       case('l')
          sh(1:1)=shape(arr1dl)
          lb(1:1)=lbound(arr1dl); ub(1:1)=ubound(arr1dl)
          str='('//fmtdi//','//trim(frmt)//')'
          do i1=lb(1),ub(1)
             if (twrite) then
                lt=arr1dl(i1)
                write(unit=unt,fmt=trim(str)) i1,lt
             else
                read(unit=unt,fmt=*) j1,lt
                arr1dl(j1)=lt
             end if
          end do
       case('i')
          sh(1:1)=shape(arr1di)
          lb(1:1)=lbound(arr1di); ub(1:1)=ubound(arr1di)
          str='('//fmtdi//','//trim(frmt)//')'
          do i1=lb(1),ub(1)
             if (twrite) then
                it=arr1di(i1)
                write(unit=unt,fmt=trim(str)) i1,it
             else
                read(unit=unt,fmt=*) j1,it
                arr1di(j1)=it
             end if
          end do
       case('r')
          sh(1:1)=shape(arr1dr)
          lb(1:1)=lbound(arr1dr); ub(1:1)=ubound(arr1dr)
          str='('//fmtdi//'," ",'//trim(frmt)//')'
          do i1=lb(1),ub(1)
             if (twrite) then
                rt=arr1dr(i1)
                write(unit=unt,fmt=trim(str)) i1,rt
             else
                read(unit=unt,fmt=*) j1,rt
                arr1dr(j1)=rt
             end if
          end do
       case('c')
          sh(1:1)=shape(arr1dc)
          lb(1:1)=lbound(arr1dc); ub(1:1)=ubound(arr1dc) 
          str='('//fmtdi//'," ",'//trim(frmt)//'," ",'//trim(frmt)//')'
          if (tparent) str='('//fmtdi//'," ("'//trim(frmt)//',","'// &
               trim(frmt)//',")")'
          do i1=lb(1),ub(1)
             if (twrite) then
                zt=arr1dc(i1)
                write(unit=unt,fmt=trim(str)) i1,zt
             else
                read(unit=unt,fmt=*) j1,zt
                arr1dc(j1)=zt
             end if
          end do
       end select
    case(2)
       !----------------------!
       !     2 dimensions     !
       !----------------------!
       select case(iot)
       case('l')
          sh(1:2)=shape(arr2dl)
          lb(1:2)=lbound(arr2dl); ub(1:2)=ubound(arr2dl)
          str='(2'//fmtdi//','//trim(frmt)//')'
          do i1=lb(1),ub(1)
             do i2=lb(2),ub(2)
                if (twrite) then
                   lt=arr2dl(i1,i2)
                   write(unit=unt,fmt=trim(str)) &
                        i1,i2,lt
                else
                   read(unit=unt,fmt=*) j1,j2,lt
                   arr2dl(j1,j2)=lt
                end if
             end do
          end do
       case('i')
          sh(1:2)=shape(arr2di)
          lb(1:2)=lbound(arr2di); ub(1:2)=ubound(arr2di) 
          str='(2'//fmtdi//','//trim(frmt)//')'
          do i1=lb(1),ub(1)
             do i2=lb(2),ub(2)
                if (twrite) then
                   it=arr2di(i1,i2)
                   write(unit=unt,fmt=trim(str)) &
                        i1,i2,it
                else
                   read(unit=unt,fmt=*) j1,j2,it
                   arr2di(j1,j2)=it
                end if
             end do
          end do
       case('r')
          sh(1:2)=shape(arr2dr)
          lb(1:2)=lbound(arr2dr); ub(1:2)=ubound(arr2dr)
          str='(2'//fmtdi//'," ",'//trim(frmt)//')'
          do i1=lb(1),ub(1)
             do i2=lb(2),ub(2)
                if (twrite) then
                   rt=arr2dr(i1,i2)
                   write(unit=unt,fmt=trim(str)) &
                        i1,i2,rt
                else
                   read(unit=unt,fmt=*) j1,j2,rt
                   arr2dr(j1,j2)=rt
                end if
             end do
          end do
       case('c')
          sh(1:2)=shape(arr2dc)
          lb(1:2)=lbound(arr2dc); ub(1:2)=ubound(arr2dc)
          str='(2'//fmtdi//'," ",'//trim(frmt)//'," ",'//trim(frmt)//')'
          if (tparent) str='(2'//fmtdi//'," ("'//trim(frmt)//',","'// &
               trim(frmt)//',")")'
          do i1=lb(1),ub(1)
             do i2=lb(2),ub(2)
                if (twrite) then
                   zt=arr2dc(i1,i2)
                   if (tparent) then
                      write(unit=unt,fmt=trim(str)) i1,i2,dble(zt), &
                           aimag(zt)
                   else
                      write(unit=unt,fmt=trim(str)) i1,i2,zt
                   end if
                else
                   read(unit=unt,fmt=*) j1,j2,zt
                   arr2dc(j1,j2)=zt
                end if
             end do
          end do
       end select
    case(3)
       !----------------------!
       !     3 dimensions     !
       !----------------------!
       select case(iot)
       case('l')
          sh(1:3)=shape(arr3dl)
          lb(1:3)=lbound(arr3dl); ub(1:3)=ubound(arr3dl)
          str='(3'//fmtdi//','//trim(frmt)//')'
          do i1=lb(1),ub(1)
             do i2=lb(2),ub(2)
                do i3=lb(3),ub(3)
                   if (twrite) then
                      lt=arr3dl(i1,i2,i3)
                      write(unit=unt,fmt=trim(str)) &
                           i1,i2,i3,lt
                   else
                      read(unit=unt,fmt=*) j1,j2,j3,lt
                      arr3dl(j1,j2,j3)=lt
                   end if
                end do
             end do
          end do
       case('i')
          sh(1:3)=shape(arr3di)
          lb(1:3)=lbound(arr3di); ub(1:3)=ubound(arr3di)
          str='(3'//fmtdi//','//trim(frmt)//')'
          do i1=lb(1),ub(1)
             do i2=lb(2),ub(2)
                do i3=lb(3),ub(3)
                   if (twrite) then
                      it=arr3di(i1,i2,i3)
                      write(unit=unt,fmt=trim(str)) &
                           i1,i2,i3,it
                   else
                      read(unit=unt,fmt=*) j1,j2,j3,it
                      arr3di(j1,j2,j3)=it
                   end if
                end do
             end do
          end do
       case('r')
          sh(1:3)=shape(arr3dr)
          lb(1:3)=lbound(arr3dr); ub(1:3)=ubound(arr3dr)
          str='(3'//fmtdi//'," ",'//trim(frmt)//')'
          do i1=lb(1),ub(1)
             do i2=lb(2),ub(2)
                do i3=lb(3),ub(3)
                   if (twrite) then
                      rt=arr3dr(i1,i2,i3)
                      write(unit=unt,fmt=trim(str)) &
                           i1,i2,i3,rt
                   else
                      read(unit=unt,fmt=*) j1,j2,j3,rt
                      arr3dr(j1,j2,j3)=rt
                   end if
                end do
             end do
          end do
       case('c')
          sh(1:3)=shape(arr3dc)
          lb(1:3)=lbound(arr3dc); ub(1:3)=ubound(arr3dc)
          str='(3'//fmtdi//'," ",'//trim(frmt)//'," ",'//trim(frmt)//')'
          if (tparent) str='(3'//fmtdi//'," ("'//trim(frmt)//',","'// &
               trim(frmt)//',")")'
          do i1=lb(1),ub(1)
             do i2=lb(2),ub(2)
                do i3=lb(3),ub(3)
                   if (twrite) then
                      zt=arr3dc(i1,i2,i3)
                      write(unit=unt,fmt=trim(str)) &
                           i1,i2,i3,zt
                   else
                      read(unit=unt,fmt=*) j1,j2,j3,zt
                      arr3dc(j1,j2,j3)=zt
                   end if
                end do
             end do
          end do
       end select
    case(4)
       !----------------------!
       !     4 dimensions     !
       !----------------------!
       select case(iot)
       case('l')
          sh(1:4)=shape(arr4dl)
          lb(1:4)=lbound(arr4dl); ub(1:4)=ubound(arr4dl)
          str='(4'//fmtdi//','//trim(frmt)//')'
          do i1=lb(1),ub(1)
             do i2=lb(2),ub(2)
                do i3=lb(3),ub(3)
                   do i4=lb(4),ub(4)
                      if (twrite) then
                         lt=arr4dl(i1,i2,i3,i4)
                         write(unit=unt,fmt=trim(str)) &
                              i1,i2,i3,i4,lt
                      else
                         read(unit=unt,fmt=*) j1,j2,j3,j4,lt
                         arr4dl(j1,j2,j3,4)=lt
                      end if
                   end do
                end do
             end do
          end do

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
       case('i')
          sh(1:3)=shape(arr3di)
          lb(1:3)=lbound(arr3di); ub(1:3)=ubound(arr3di)
          str='(3'//fmtdi//','//trim(frmt)//')'
          do i1=lb(1),ub(1)
             do i2=lb(2),ub(2)
                do i3=lb(3),ub(3)
                   if (twrite) then
                      it=arr3di(i1,i2,i3)
                      write(unit=unt,fmt=trim(str)) &
                        i1,i2,i3,it
                   else
                      read(unit=unt,fmt=*) j1,j2,j3,it
                      arr3di(j1,j2,j3)=it
                   end if
                end do
             end do
          end do
       case('r')
          sh(1:3)=shape(arr3dr)
          lb(1:3)=lbound(arr3dr); ub(1:3)=ubound(arr3dr)
          str='(3'//fmtdi//'," ",'//trim(frmt)//')'
          do i1=lb(1),ub(1)
             do i2=lb(2),ub(2)
                do i3=lb(3),ub(3)
                   if (twrite) then
                      rt=arr3dr(i1,i2,i3)
                      write(unit=unt,fmt=trim(str)) &
                        i1,i2,i3,rt
                   else
                      read(unit=unt,fmt=*) j1,j2,j3,rt
                      arr3dr(j1,j2,j3)=rt
                   end if
                end do
             end do
          end do
       case('c')
          sh(1:3)=shape(arr3dc)
          lb(1:3)=lbound(arr3dc); ub(1:3)=ubound(arr3dc)
          str='(3'//fmtdi//'," ",'//trim(frmt)//'," ",'//trim(frmt)//')'
          if (tparent) str='(3'//fmtdi//'," ("'//trim(frmt)//',","'// &
               trim(frmt)//',")")'
          do i1=lb(1),ub(1)
             do i2=lb(2),ub(2)
                do i3=lb(3),ub(3)
                   if (twrite) then
                      zt=arr3dc(i1,i2,i3)
                      write(unit=unt,fmt=trim(str)) &
                        i1,i2,i3,zt
                   else
                      read(unit=unt,fmt=*) j1,j2,j3,zt
                      arr3dc(j1,j2,j3)=zt
                   end if
                end do
             end do
          end do
       end select









       goto 10
    case(5)
       goto 10
    case(6)
       goto 10
    case(7)
       goto 10
    end select
    goto 20
10  continue
    write(*,*)
    write(*,*) 'Error(ioarray::ioarr): currently arrays limited to 3 dimensions'
    write(*,*)
    stop
20 continue
  end subroutine ioarr

  ! determine IO (read/write) action
  logical function ioaction(ioa)
    character(*), intent(in) :: ioa
    select case(trim(adjustl(ioa)))
    case('r','R','read','Read','READ')
       ioaction=.false.
    case('w','W','write','Write','WRITE')
       ioaction=.true.
    case default
       write(*,*) 'Error(ioarray::ioaction): unknown IO action: '// &
            trim(adjustl(ioa))
       stop
    end select
  end function ioaction

!!$  character(1) function iotype(iot)
!!$    character(*), intent(in) :: iot
!!$    select case(trim(adjustl(iot)))
!!$    case ('l','L','logical','Logical','LOGICAL')
!!$       iotype='l'
!!$    case ('i','I','integer','Integer','INTEGER')
!!$       iotype='i'
!!$    case ('r','R','real','Real','REAL')
!!$       iotype='r'
!!$    case ('c','C','complex','Complex','COMPLEX')
!!$       iotype='c'
!!$    case default
!!$       write(*,*) 'Error(ioarray::iotype): unknown IO type: '// &
!!$            trim(adjustl(iot))
!!$       stop
!!$    end select
!!$  end function iotype

end module ioarray
