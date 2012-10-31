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
      integer(4) :: ik, ik1, ik2, isym
      integer(4) :: ib, nb, nk
      integer(8) :: recl
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
      
      inquire(IoLength=Recl) nkp1, kvecs1(1:3,1), ibgw, nbgw, &
        eqp1(ibgw:nbgw,1), eks1(ibgw:nbgw,1)
      open (70, File=file, Action='READ', Form='UNFORMATTED', &
     &  Access='DIRECT', Recl=Recl)
      
      do ik = 1, nkp1
        read(70, Rec=ik) nk, kvecs1(:,ik), ib, nb, &
          eqp1(ibgw:nbgw,ik), eks1(ibgw:nbgw,ik)
      end do ! ik
      close(70)
      
      if (.not.present(nkp2)) return
      
!-----------------------------------------------------------------------------
!     Interpolate the energies
!-----------------------------------------------------------------------------      
      
      if ((ibgw.ne.1).or.(nbgw.lt.nstsv)) then
        write(*,*)
        write(*,*)'WARNING(getevalqp):'
        write(*,*)'  Quasiparticle energies has been calculated for the interval of bands:'
        write(*,*)'  [ibgw,nbgw]=[', ibgw, nbgw,']' 
        write(*,*)'  If it creates problems - regenerate your EVALQP.OUT.'
        write(*,*)
      end if

      allocate(de1(nkp1,ibgw:nbgw),de2(nkpt,ibgw:nbgw))
      do ik = 1, nkp1
         de1(ik,:)=cmplx(eqp1(ibgw:nbgw,ik)-eks1(ibgw:nbgw,ik),0.d0,8)
      enddo

      de2(:,:)=zzero
      call fourintp(de1,nkp1,kvecs1,de2,nkp2,vkl,nbgw-ibgw+1)

      do ib = ibgw, min(nbgw,nstsv)
         do ik = 1, nkpt
            evalsv(ib,ik)=evalsv(ib,ik)+real(de2(ik,ib))
         enddo 
      enddo
      
      deallocate(de1,de2)
      deallocate(kvecs1,eqp1,eks1)

      return
end subroutine
!EOC
