!
!BOP
! !ROUTINE: getevalqp
! !INTERFACE:
!
subroutine getevalqp(nkp2,kvecs2,eqp2)
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
      
      integer, intent(in) :: nkp2
      real(8), intent(in) :: kvecs2(3,nkp2)
      real(8), intent(inout):: eqp2(nstsv,nkp2)

!     local variables
      logical :: exist
      integer(4) :: ik, ik1, ik2, isym, lspl, iv(3)
      integer(4) :: ib, nb, nk
      integer(4) :: recl
      integer :: nkpt0, nstsv0 
      real(8) :: s(3,3), v1(3), v2(3), t1
      character (256) :: file

      integer, allocatable :: map(:)
      complex(8), allocatable :: de1(:,:), de2(:,:)
  real(8),    Allocatable :: vkl0(:,:), ehf(:,:)
  complex(8), allocatable :: e0(:,:),e1(:,:)

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
        read(70, Rec=ik) nk, ib, nb, kvecs1(:,ik), & 
        &  eqp1(ibgw:nbgw,ik), eks1(ibgw:nbgw,ik)
        
        ! debugging info
        !write(fgw,*) '# ik    kvecs1    ibgw,    nbgw'
        !write(fgw,*) ik, kvecs1(:,ik), ib, nb
        !write(fgw,*) '# ib    eqp1    eks1'
        !do ib = ibgw, nbgw
        !  write(fgw,*) ib, eqp1(ib,ik), eks1(ib,ik)
        !end do 
                
      end do ! ik
      close(70)

! read DFT eigenvalues from EVALHF.OUT in case of HF-Hybrids
if (associated(input%groundstate%Hybrid)) then
    if (input%groundstate%Hybrid%exchangetypenumber== 1) then 
    file='EVALHF.OUT'
    inquire(File=file,Exist=exist)
    if (.not.exist) then
      write(*,*)'ERROR(getevalqp): File EVALHF.OUT does not exist!'
      stop
    end if
    inquire(IoLength=Recl) nkpt0, nstsv0
    open (70, File=file, Action='READ', Form='UNFORMATTED', &
    &  Access='DIRECT', Recl=Recl)
    read(70, Rec=1) nkpt0, nstsv0
    close(70)
    nstsv=min(nstsv,nstsv0)
    allocate(vkl0(3,nkpt0))
    allocate(ehf(nstsv0,nkpt0))
    allocate(e0(nkpt0,nstsv))
    allocate(e1(nkpt,nstsv))
    inquire(IoLength=Recl) nkpt0, nstsv0, vkl0(:,1), ehf(:,1)
    open (70, File=file, Action='READ', Form='UNFORMATTED', &
    &  Access='DIRECT', Recl=Recl)
    do ik = 1, nkpt0
      read(70, Rec=ik) nkpt0, nstsv0, vkl0(:,ik), ehf(1:,ik)
    end do ! ik
    close(70)
    ! Fourier Interpolation for DFT eigenvalues
    do ik = 1, nkpt0
      e0(ik,1:nstsv)=cmplx(ehf(1:nstsv,ik),0.d0,8)
    enddo
    e1(:,:)=zzero
    call fourintp(e0,nkpt0,vkl0,e1,nkp2,vkl,nstsv)
    do ik2 = 1, nkp2
        do ib =  1,nstsv 
          eqp2(ib,ik2) = e1(ik2,ib)
      end do
    end do
    deallocate(vkl0,ehf,e0,e1)
  end if
end if

!-------------------------------------------------------------------------------
!     Data-set consistency check
!-------------------------------------------------------------------------------      

      if ((ibgw.ne.1).or.(nbgw.lt.nstsv)) then
        write(*,*)
        write(*,*)'WARNING(getevalqp):'
        write(*,*)'  Quasiparticle energies has been calculated for the interval of bands:'
        write(*,*)'  [ibgw,nbgw]=[', ibgw, nbgw,']'
        write(*,*)'  Check whether it is that you need ...'
        write(*,*)
      end if

! Special case of only one k-point (molecules)
      if (nkp1==1) then
        if (nkp2==1) then
          do ib = ibgw, min(nbgw,nstsv)
            eqp2(ib,1) = eqp1(ib,1)
          end do
          goto 10
        else
          write(*,*) 'ERROR(getevalqp):' 
          write(*,*) '  Interpolation is not possible!'
          write(*,*) '  EVALQP.OUT file contains data only for a single k-point.'
          stop
        end if
      end if 

! Check if k-sets are (symmetry) identical
      allocate(map(nkp2))
      map(:) = 0
      do ik2 = 1, nkp2
        do isym = 1, nsymcrys
          lspl = lsplsymc(isym)
          s(:, :) = dble(symlat(:, :, lspl))
          call r3mtv(s, kvecs2(:,ik2), v2)
          call r3frac(input%structure%epslat, v1, iv)
          do ik1 = 1, nkp1
            v1(:) = kvecs1(:,ik1)
            call r3frac (input%structure%epslat, v1, iv)
            t1 = Abs(v2(1)-v1(1)) + Abs(v2(2)-v1(2)) + Abs(v2(3)-v1(3))
            if (t1 > input%structure%epslat) then
              ! not equivalent, proceed with interpolation
              deallocate(map)
              goto 20
            else
              map(ik2) = ik1 
            end if
          end do
        end do ! isym
      end do ! ik2

! reached this points = sets are equivalent
      do ik2 = 1, nkp2
        do ib = ibgw, min(nbgw,nstsv)
          eqp2(ib,ik2) = eqp1(ib,map(ik2))
        end do
      end do
      deallocate(map)   
      goto 10

!-----------------------------------------------------------------------------
!     Interpolate the energies
!-----------------------------------------------------------------------------      
20    continue
      
      allocate(de1(nkp1,ibgw:nbgw),de2(nkpt,ibgw:nbgw))
      do ik = 1, nkp1
         de1(ik,:)=cmplx(eqp1(ibgw:nbgw,ik)-eks1(ibgw:nbgw,ik),0.d0,8)
      enddo

      de2(:,:)=zzero
      call fourintp(de1,nkp1,kvecs1,de2,nkp2,vkl,nbgw-ibgw+1)

      do ib = ibgw, min(nbgw,nstsv)
         do ik = 1, nkpt
            eqp2(ib,ik)=eqp2(ib,ik)+real(de2(ik,ib))
         enddo 
      enddo
      deallocate(de1,de2)

10    continue
      deallocate(kvecs1,eqp1,eks1)
      return

end subroutine
!EOC
