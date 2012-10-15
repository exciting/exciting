!
!BOP
! !ROUTINE: getevalqp
! !INTERFACE:
!
subroutine getevalqp(nkp2)
! !USES:
      use modinput
      use modmain
      use modgw,    only: ibgw, nbgw, nkp1, kvecs1, eks1, eqp1

! !DESCRIPTION:
!   The file where the quasiparticle energies are stored is {\tt EVALQP.OUT}.
!   It is a direct-access binary file, the record length of which can be
!   determined with the help of the array sizes and data type information.
!
!EOP
!BOC
      implicit none
      
      integer, optional, intent(in) :: nkp2
      ! if not specified - only read the QP energies
      ! otherwise        - read and interpolate QPE

!     local variables
      logical :: exist
      integer :: ik, ib
      integer :: recl
      real(8) :: v(3), t1
      character (256) :: file
      complex(8), allocatable :: de1(:,:), de2(:,:)

!-----------------------------------------------------------------------------
!     Just read the file
!-----------------------------------------------------------------------------      

      file='EVALQP.OUT'
      inquire(File=file,Exist=exist)
      if (.not.exist) then
        write(*,*)'ERROR(getevalqp): File EVALQP.OUT does not exist!'
        stop
      end if
      
      inquire(IoLength=Recl) nkp1, v, ibgw, nbgw
      open (70, File=file, Action='READ', Form='UNFORMATTED', &
     &  Access='DIRECT', Recl=Recl)
      read(70, Rec=1) nkp1, v, ibgw, nbgw
      close(70)
      
      allocate(kvecs1(1:3,nkp1))
      allocate(eqp1(ibgw:nbgw,nkp1))
      allocate(eks1(ibgw:nbgw,nkp1))
      inquire(IoLength=Recl) nkp1, kvecs1(:,1), ibgw, nbgw, eqp1(:,1), eks1(:,1)
      open (70, File=file, Action='READ', Form='UNFORMATTED', &
     &  Access='DIRECT', Recl=Recl)
      do ik = 1, nkp1
        read(70, Rec=ik) nkp1, kvecs1(:,ik), ibgw, nbgw, eqp1(:,ik), eks1(:,ik)
      end do ! ik
      close(70)
      
!-----------------------------------------------------------------------------
!     Interpolate the energies
!-----------------------------------------------------------------------------      
      
      if (.not.present(nkp2)) nkp2=nkpt
      
      if ((ibgw.ne.1).or.(nbgw.lt.nstsv)) then
        write(*,*)
        write(*,*)'WARNING(getevalqp):'
        write(*,*)'  Quasiparticle energies has been calculated for the interval of bands:'
        write(*,*)'  [ibgw,nbgw]=[', ibgw, nbgw,']' 
        write(*,*)'  Regenerate your EVALQP.OUT file if it is needed'
        write(*,*)
      end if
!
!     if k-point grid is the same for both GW and BSE
!
      if (nkp1.eq.nkp2) then
        t1=0.0d0
        do ik = 1, nkp2
          t1 = t1 + abs(vkl(1,ik)-kvecs1(1,ik)) + &
         &          abs(vkl(2,ik)-kvecs1(2,ik)) + &
         &          abs(vkl(3,ik)-kvecs1(3,ik))
        end do
        if (t1.lt.input%structure%epslat) then
          evalsv(ibgw:nbgw,1:nkp1)=eqp1(ibgw:nbgw,1:nkp1)
          return
        end if
      end if

!--------------------------------------------------
!       Fourier interpolation
!--------------------------------------------------

      allocate(de1(ibgw:nbgw,nkp1),de2(ibgw:nbgw,nkp2))

      do ik = 1, nkp1
         de1(ibgw:nbgw,ik)=cmplx(eqp1(ibgw:nbgw,ik)-eks1(ibgw:nbgw,ik),0.d0,8)
      enddo

      de2(:,:)=0.0d0
      call fourintp(nbgw-ibgw+1,nkp1,kvecs1,de1,nkp2,vkl,de2)

      do ib = ibgw, min(nbgw,nstsv)
         do ik = 1, nkp2
            evalsv(ib,ik)=evalsv(ib,ik)+real(de2(ib,ik))
         enddo 
      enddo
      
      deallocate(de1,de2)
      deallocate(kvecs1,eqp1,eks1)

      return
end subroutine
!EOC
