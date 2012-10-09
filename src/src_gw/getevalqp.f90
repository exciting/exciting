!
!BOP
! !ROUTINE: getevalqp
! !INTERFACE:
!
subroutine getevalqp (vpl, evalsvp)
! !USES:
      use modmain
      use modinput
! !DESCRIPTION:
!   The file where the quasiparticle energies are stored is {\tt EVALQP.OUT}.
!   It is a direct-access binary file, the record length of which can be
!   determined with the help of the array sizes and data type information.
!   One record of this file has the following structure
!
!EOP
!BOC
      implicit none
  ! arguments
      real(8), intent (In) :: vpl(3)
      real(8), intent (Out):: evalsvp(nstsv)
  ! local variables
      logical :: exist
      integer :: isym, ik
      integer :: recl, nstsv_
      real(8) :: vkl_ (3), t1
      character (256) :: file
  ! to access arrays for only a subset of bands
      Real (8), Allocatable :: evalsv_(:)

  ! find the k-point number
      Call findkpt(vpl,isym,ik)

  ! find the record length
      Inquire(IoLength=Recl) vkl_, nstsv_, evalsvp

      file='EVALQP.OUT'
      inquire(File=file, Exist=Exist)
      if (exist) then
        open (70, File=file, Action='READ', Form='UNFORMATTED', &
       &  Access='DIRECT', Recl=Recl)
      else
        write(*,*)'ERROR(getevalqp) File EVALQP.OUT does not exist!'
        stop
      end if
      
      read(70, Rec=1) vkl_, nstsv_
      close(70)
      if (nstsv .gt. nstsv_) Then
         write (*,*)
         write (*, '("Error(getevalqp): invalid nstsv")')
         write (*, '(" current      : ",I8)') nstsv
         write (*, '(" EVALSVQP.OUT : ",I8)') nstsv_
         write (*,*)
         stop
      end if

      allocate(evalsv_(nstsv_))
      inquire(IoLength=Recl) vkl_, nstsv_, evalsv_
      open(70, File=file, Action='READ', Form='UNFORMATTED', &
     &  access='DIRECT', Recl=Recl)
      read(70, Rec=ik) vkl_, nstsv_, evalsv_
  ! retreive subset
      evalsvp (:) = evalsv_ (:nstsv)
      deallocate(evalsv_)
      close (70)
!
      t1 = Abs(vkl(1, ik)-vkl_(1)) + &
     &     Abs(vkl(2, ik)-vkl_(2)) + &
     &     abs(vkl(3, ik)-vkl_(3))
      if (t1 .gt. input%structure%epslat) then
         write (*,*)
         write (*, '("Error(getevalsv): differing vectors for k-point "&
        &,i8)') ik
         write (*, '(" current    : ",3G18.10)') vkl (:, ik)
         write (*, '(" EVALQP.OUT : ",3G18.10)') vkl_
         write (*,*)
         stop
      end if
!
      return
end subroutine
!EOC
