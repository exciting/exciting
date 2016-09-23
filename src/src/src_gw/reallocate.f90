! !MODULE: reallocate
       module reallocate

!
! !INTERFACE:
         interface doreallocate
            module procedure doreallocate_r8_d1
            module procedure doreallocate_r8_d2
            module procedure doreallocate_r8_d3
            module procedure doreallocate_r8_d4
            module procedure doreallocate_r8_d5
            module procedure doreallocate_r8_d6
            module procedure doreallocate_i4_d1
            module procedure doreallocate_i4_d2
            module procedure doreallocate_i4_d3
            module procedure doreallocate_z8_d1
            module procedure doreallocate_z8_d2
            module procedure doreallocate_z8_d3
            module procedure hugo     !   ;)
          end interface
          
! !DESCRIPTION:
!
! Performs the reallocation in memory of a matrix that has been
! "overallocated"
!
!EOP
        contains
!BOP
!
! !IROUTINE: doreallocate_r8_d1
!
! !INTERFACE:
          subroutine doreallocate_r8_d1(tf, newdimension)
!
! !INPUT PARAMETERS:
            implicit none

            real(8), pointer :: tf(:)  ! vector to reallocate
            
            integer(4) :: newdimension ! new dimension of the vector
!
! !DESCRIPTION:
!
! Reallocates a real(8) pointer vector of range 1
!            
!            
! !LOCAL VARIABLES:            
            real(8), pointer :: hilfsfeld(:)

            integer(4) :: min1

!EOP
!BOC
            allocate(hilfsfeld(newdimension))
            hilfsfeld=0.0D0
!            
!           It is enough to copy only once            
!
            min1=min(newdimension,size(tf,1))
            hilfsfeld(1:min1)=tf(1:min1)
            deallocate(tf)
!            
!           the pointer is redirected to the new field, not reallocated
!
            tf=>hilfsfeld
          end subroutine
!EOC

!BOP
!
! !IROUTINE: doreallocate_r8_d2
!
! !INTERFACE:
          subroutine doreallocate_r8_d2(tf, newdimension1, newdimension2)

!
! !INPUT PARAMETERS:          
            implicit none
            
            real(8), pointer :: tf(:,:) ! the matrix to reallocate
            
            integer(4) :: newdimension1, newdimension2 ! new dimensions
            
!
! !DESCRIPTION:
!
! Reallocates a real(8) pointer matrix of range 2
!            
! !LOCAL VARIABLES:
            real(8), pointer :: hilfsfeld(:,:)
            
            integer(4) :: min1, min2
!EOP
!BOC
            allocate(hilfsfeld(newdimension1,newdimension2))
            hilfsfeld=0.0D0
            min1=min(newdimension1,size(tf,1))
            min2=min(newdimension2,size(tf,2))
            hilfsfeld(1:min1,1:min2)=tf(1:min1,1:min2)
            deallocate(tf)
            tf=>hilfsfeld
          end subroutine
!EOC
!BOP
!
! !IROUTINE: doreallocate_r8_d3
!
! !INTERFACE:
          subroutine doreallocate_r8_d3(tf, newdim1, newdim2,          &
     &                                  newdim3)

!
! !INPUT PARAMETERS:     
            implicit none

            real(8), pointer :: tf(:,:,:) 
            
            integer(4) :: newdim1, newdim2, newdim3

!
! !LOCAL VARIABLES:            
            
            real(8), pointer :: hilfsfeld(:,:,:)
            
            integer(4) :: min1, min2, min3
!
!EOP
!BOC            
            allocate(hilfsfeld(newdim1,newdim2,newdim3))
            hilfsfeld=0.0D0
            min1=min(newdim1,size(tf,1))
            min2=min(newdim2,size(tf,2))
            min3=min(newdim3,size(tf,3))
            hilfsfeld(1:min1,1:min2,1:min3)=tf(1:min1,1:min2,1:min3)
            deallocate(tf)
            tf=>hilfsfeld
          end subroutine
!EOC
!BOP
!
! !IROUTINE: doreallocate_r8_d4
!
! !INTERFACE:
          subroutine doreallocate_r8_d4(tf, newdim1, newdim2, newdim3,&
     &    newdim4)

! !INPUT PARAMETERS:
            implicit none

            integer(4) :: newdim1,newdim2,newdim3,newdim4

            real(8), pointer :: tf(:,:,:,:)

! !LOCAL VARIABLES:
            real(8), pointer :: hilfsfeld(:,:,:,:)

            integer(4) :: min1, min2, min3, min4

!EOP
!BOC            
            allocate(hilfsfeld(newdim1,newdim2,newdim3,newdim4))
            hilfsfeld=0.0D0
            min1=min(newdim1,size(tf,1))
            min2=min(newdim2,size(tf,2))
            min3=min(newdim3,size(tf,3))
            min4=min(newdim4,size(tf,4))
            hilfsfeld(1:min1,1:min2,1:min3,1:min4)=&
     &        tf(1:min1,1:min2,1:min3,1:min4)
            deallocate(tf)
            tf=>hilfsfeld
          end subroutine
!EOC
!BOP
!
! !IROUTINE: doreallocate_r8_d5
!
! !INTERFACE:
          subroutine doreallocate_r8_d5(tf, newdim1, newdim2, newdim3,&
     &    newdim4,newdim5)

! !INPUT PARAMETERS:     
            implicit none
            
            real(8), pointer :: tf(:,:,:,:,:)
            
            integer(4) :: newdim1,newdim2,newdim3,newdim4, newdim5
            
!
! !LOCAL VARIABLES:
            integer(4) :: min1, min2, min3, min4, min5

            real(8), pointer :: hilfsfeld(:,:,:,:,:)
!EOP
!BOC
            allocate(hilfsfeld(newdim1,newdim2,newdim3,newdim4,newdim5))
            hilfsfeld=0.0D0
            min1=min(newdim1,size(tf,1))
            min2=min(newdim2,size(tf,2))
            min3=min(newdim3,size(tf,3))
            min4=min(newdim4,size(tf,4))
            min5=min(newdim5,size(tf,5))
            hilfsfeld(1:min1,1:min2,1:min3,1:min4,1:min5)=&
     &        tf(1:min1,1:min2,1:min3,1:min4,1:min5)
            deallocate(tf)
            tf=>hilfsfeld
          end subroutine
!EOC
!BOP
!
! !IROUTINE: doreallocate_r8_d6
!
! !INTERFACE:
          subroutine doreallocate_r8_d6(tf, newdim1, newdim2, newdim3,&
     &    newdim4,newdim5,newdim6)
!
! !INPUT PARAMETERS:     
            implicit none

            real(8), pointer :: tf(:,:,:,:,:,:)

            integer(4) :: newdim1, newdim2, newdim3
            integer(4) :: newdim4, newdim5, newdim6
!
! !LOCAL VARIABLES:            
            real(8), pointer :: hilfsfeld(:,:,:,:,:,:)

            integer(4) :: min1, min2, min3, min4, min5, min6
!
!EOP
!BOC            
            allocate(hilfsfeld(newdim1,newdim2,newdim3,newdim4,        &
     &               newdim5,newdim6))
            hilfsfeld=0.0D0
            min1=min(newdim1,size(tf,1))
            min2=min(newdim2,size(tf,2))
            min3=min(newdim3,size(tf,3))
            min4=min(newdim4,size(tf,4))
            min5=min(newdim5,size(tf,5))
            min6=min(newdim6,size(tf,6))
            hilfsfeld(1:min1,1:min2,1:min3,1:min4,1:min5,1:min6)=      &
     &        tf(1:min1,1:min2,1:min3,1:min4,1:min5,1:min6)
            deallocate(tf)
            tf=>hilfsfeld
          end subroutine
!EOC
!BOP
!
! !IROUTINE: doreallocate_i4_d1
!
! !INTERFACE:
          subroutine doreallocate_i4_d1(tf, newdimension)
!
! !INPUT PARAMETERS:          
            implicit none
            
            integer(4), pointer :: tf(:)
            
            integer(4) :: newdimension
!
! !LOCAL VARIABLES:
            integer(4), pointer :: hilfsfeld(:)
            
            integer(4) :: min1
!
!EOP
!BOC            
            allocate(hilfsfeld(newdimension))
            hilfsfeld=0
            min1=min(newdimension,size(tf,1))
            hilfsfeld(1:min1)=tf(1:min1)
            deallocate(tf)
            tf=>hilfsfeld
          end subroutine
!EOC
!BOP
!
! !IROUTINE: doreallocate_I4_d2
!
! !INTERFACE:
          subroutine doreallocate_i4_d2(tf, newdimension1, newdimension2)
!
! !INPUT PARAMETERS:
          
            implicit none

            integer(4), pointer :: tf(:,:)

            integer(4) :: newdimension1, newdimension2
!
! !LOCAL VARIABLES:            

            integer(4), pointer :: hilfsfeld(:,:)

            integer(4) :: min1, min2

!EOP
!BOC            
            allocate(hilfsfeld(newdimension1,newdimension2))
            hilfsfeld=0
            min1=min(newdimension1,size(tf,1))
            min2=min(newdimension2,size(tf,2))
            hilfsfeld(1:min1,1:min2)=tf(1:min1,1:min2)
            deallocate(tf)
            tf=>hilfsfeld
          end subroutine
!EOC
!BOP
!
! !IROUTINE: doreallocate_i4_d3
!
! !INTERFACE:
          subroutine doreallocate_i4_d3(tf, newdim1, newdim2, newdim3)
!
! !INPUT PARAMETERS:          
            implicit none
            
            integer(4), pointer :: tf(:,:,:)
            
            integer(4) :: newdim1, newdim2, newdim3

! !LOCAL VARIABLES:            
            
            integer(4), pointer :: hilfsfeld(:,:,:)
            
            integer(4) :: min1, min2, min3

!EOP
!BOC
            allocate(hilfsfeld(newdim1,newdim2,newdim3))
            hilfsfeld=0
            min1=min(newdim1,size(tf,1))
            min2=min(newdim2,size(tf,2))
            min3=min(newdim3,size(tf,3))
            hilfsfeld(1:min1,1:min2,1:min3)=tf(1:min1,1:min2,1:min3)
            deallocate(tf)
            tf=>hilfsfeld
          end subroutine


!EOC

!BOP
!
! !IROUTINE: doreallocate_z8_d1
!
! !INTERFACE:
          subroutine doreallocate_z8_d1(tf, newdimension)
!
! !INPUT PARAMETERS:
            implicit none

            complex(8), pointer :: tf(:)  ! vector to reallocate
            
            integer(4) :: newdimension ! new dimension of the vector
!
! !DESCRIPTION:
!
! Reallocates a complex(8) pointer vector of range 1
!            
!            
! !LOCAL VARIABLES:            
            complex(8), pointer :: hilfsfeld(:)

            integer(4) :: min1

!EOP
!BOC
            allocate(hilfsfeld(newdimension))
            hilfsfeld=0.0D0
!            
!           It is enough to copy only once            
!
            min1=min(newdimension,size(tf,1))
            hilfsfeld(1:min1)=tf(1:min1)
            deallocate(tf)
!            
!           the pointer is redirected to the new field, not reallocated
!
            tf=>hilfsfeld
          end subroutine
!EOC

!BOP
!
! !IROUTINE: doreallocate_z8_d2
!
! !INTERFACE:
          subroutine doreallocate_z8_d2(tf, newdim1,newdim2)

! !DESCRIPTION:
!
! Reallocates a complex(8) pointer vector of range 2
!            
!
! !INPUT PARAMETERS:
            implicit none

            integer(4), intent(in) :: newdim1 ! new dimension of the vector
            integer(4), intent(in) :: newdim2 ! new dimension of the vector
            
! !INPUT/OUTPUT PARAMETERS:

            complex(8),  pointer :: tf(:,:)  ! vector to reallocate
            
!
!            
! !LOCAL VARIABLES:            
            complex(8), pointer :: hilfsfeld(:,:)

            integer(4) :: min1,min2

!EOP
!BOC
            allocate(hilfsfeld(newdim1,newdim2))
            hilfsfeld=0.0D0
!            
!           It is enough to copy only once            
!
            min1=min(newdim1,size(tf,1))
            min2=min(newdim2,size(tf,2))
            hilfsfeld(1:min1,1:min2)=tf(1:min1,1:min2)
            deallocate(tf)
!            
!           the pointer is redirected to the new field, not reallocated
!
            tf=>hilfsfeld
          end subroutine
!EOC

!BOP
!
! !IROUTINE: doreallocate_z8_d3
!
! !INTERFACE:
          subroutine doreallocate_z8_d3(tf,newdim1,newdim2,newdim3)

! !DESCRIPTION:
!
! Reallocates a complex(8) pointer vector of range 3
!            
!
! !INPUT PARAMETERS:
            implicit none

            integer(4), intent(in) :: newdim1 ! new dimension of the vector
            integer(4), intent(in) :: newdim2 ! new dimension of the vector
            integer(4), intent(in) :: newdim3 ! new dimension of the vector
            
! !INPUT/OUTPUT PARAMETERS:

            complex(8), pointer :: tf(:,:,:)  ! vector to reallocate
            
!
!            
! !LOCAL VARIABLES:            
            complex(8), pointer :: hilfsfeld(:,:,:)

            integer(4) :: min1,min2,min3

!EOP
!BOC
            allocate(hilfsfeld(newdim1,newdim2,newdim3))
            hilfsfeld=0.0D0
!            
!           It is enough to copy only once            
!
            min1=min(newdim1,size(tf,1))
            min2=min(newdim2,size(tf,2))
            min3=min(newdim3,size(tf,3))
            hilfsfeld(1:min1,1:min2,1:min3)=tf(1:min1,1:min2,1:min3)
            deallocate(tf)
!            
!           the pointer is redirected to the new field, not reallocated
!
            tf=>hilfsfeld
          end subroutine
!EOC

!     Es gibt auch Methoden, um das Programm unleserlich zu 
!     machen :-) das sollten wir besser vermeiden ;-)
          subroutine hugo
            write(*,*) " Hier ist Hugo"
          end subroutine


        end module reallocate
