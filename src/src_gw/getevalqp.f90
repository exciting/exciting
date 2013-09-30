!
!BOP
! !ROUTINE: getevalqp
! !INTERFACE:
!
subroutine getevalqp(nkp2)
! !USES:
      use modinput
      use modmain
      use modgw,    only: fgw, ibgw, nbgw, nkp1, kvecs1, eks1, eqp1

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
      
      inquire(IoLength=Recl) nkp1, ibgw, nbgw
      open (70, File=file, Action='READ', Form='UNFORMATTED', &
     &  Access='DIRECT', Recl=Recl)
      read(70, Rec=1) nkp1, ibgw, nbgw
      close(70)
      
      allocate(kvecs1(1:3,nkp1))
      allocate(eqp1(ibgw:nbgw,nkp1))
      allocate(eks1(ibgw:nbgw,nkp1))
      
      inquire(IoLength=Recl) nkp1, ibgw, nbgw, &
     &  kvecs1(1:3,1), eqp1(ibgw:nbgw,1), eks1(ibgw:nbgw,1)
      open (70, File=file, Action='READ', Form='UNFORMATTED', &
     &  Access='DIRECT', Recl=Recl)
      
      do ik = 1, nkp1
        read(70, Rec=ik) nk, ib, nb, &
       &  kvecs1(:,ik), eqp1(ibgw:nbgw,ik), eks1(ibgw:nbgw,ik)
        
        !write(fgw,*) '# ik    kvecs1    ibgw,    nbgw'
        !write(fgw,*) ik, kvecs1(:,ik), ib, nb
        !write(fgw,*) '# ib    eqp1    eks1'
        !do ib = ibgw, nbgw
        !  write(fgw,*) ib, eqp1(ib,ik), eks1(ib,ik)
        !end do 
          
      end do ! ik
      close(70)
      
      if (.not.present(nkp2)) goto 10
      
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
      
! if both k-sets are identical - no interpolation is needed
      if (nkp1 == nkp2) then
        t1 = 0.d0
        do ik = 1, nkp1
            t1 = t1+abs(kvecs1(1,ik)-vkl(1,ik))+ &
           &        abs(kvecs1(2,ik)-vkl(2,ik))+ &
           &        abs(kvecs1(3,ik)-vkl(3,ik))
        end do
        if (t1 < 1.d-6) then
! sets are identical, just copy QP energies
            do ib = ibgw, min(nbgw,nstsv)
                do ik = 1, nkpt
                    evalsv(ib,ik)=eqp1(ib,ik)
                enddo 
            enddo
! skip interpolation            
            goto 10
        end if
      end if

! Special case of only one k-point      
      if (nkp1==1) then
        write(*,*) 'ERROR(getevalqp):' 
        write(*,*) '  Interpolation is not possible!'
        write(*,*) '  EVALQP.OUT file contains data only for a single k-point.'
        stop
      end if
      
! otherwise perform interpolation
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

10    continue
      deallocate(kvecs1,eqp1,eks1)
      return

end subroutine
!EOC
