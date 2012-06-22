!BOP
!
! !ROUTINE: plotevec
!
! !INTERFACE:
      subroutine plotevec

! !DESCRIPTION:
!
! This subroutine calculates the real space representation of the
! (L)APW+lo basis functions in the line joining the two given atoms 
! for ploting
!
! !USES:

      use modinput
      use modmain
      use modgw

! !INPUT PARAMETERS:

      implicit none

! !LOCAL VARIABLES:
      integer(4) :: ib
      integer(4) :: i
      integer(4) :: ia
      integer(4) :: ia1
      integer(4) :: ia2
      integer(4) :: irc
      integer(4) :: is
      integer(4) :: is1
      integer(4) :: is2
      integer(4) :: j
      integer(4) :: jrc
      integer(4) :: krc
      integer(4) :: np
      integer(4) :: ik, ikp
      
      character(len=64) :: filename

      real(8) :: irlen
      real(8) :: rd(3)
      real(8) :: rdlen
      real(8) :: rr(99)
      real(8), allocatable :: rs(:)
      
      complex(8), allocatable :: evec(:)

!
! !EXTERNAL ROUTINES: 
!
      external calcevecplot
      external ylm
! 
! !INTRINSIC ROUTINES: 
!
      intrinsic achar
      intrinsic dcos
      intrinsic dsin
      intrinsic mod
! 
! !REVISION HISTORY:
!
! Created: 9th. July 2004 by RGA
! Last modified: 15th. July 2004 by RGA
! Revisited 29.04.2011 by DIN
!
!EOP
!BOC
!
!     Set the name of the output file
!
      
      call boxmsg(6,'-','PLOTEVEC')
      
      ik=input%gw%iik
      ikp=indkp(input%gw%iik) 
      write(*,*) 'Parameters:'
      write(*,*) 'k-point number (iik): ', ik
      write(*,*) 'irreducible k-point number (ikp): ', ikp
      write(*,*) 'lower bound for band number (ibmin): ', input%gw%ibmin
      write(*,*) 'upper bound for band number (ibmax): ', input%gw%ibmax
      write(*,*) 'atom 1 (at1): ', input%gw%at1
      write(*,*) 'atom 2 (at2): ', input%gw%at2
      write(*,*)

      if((input%gw%at1.lt.1).or.(input%gw%at1.gt.natmtot)) stop 'atom1 is wrong'
      if((input%gw%at2.lt.1).or.(input%gw%at2.gt.natmtot)) stop 'atom2 is wrong'
!
!     Set the indexes of the two atoms
!
      ia1=0; is1=0
      ia2=0; is2=0
      do is = 1,nspecies
        do ia=1,natoms(is)
          if(idxas(ia,is).eq.input%gw%at1)then
            ia1=ia
            is1=is
          endif  
        enddo 
      enddo
      if(input%gw%at1.eq.input%gw%at2)then
        rd(1:3)=avec(1:3,1)
        is2=is1
        ia2=ia1
      else  
        do is = 1,nspecies
          do ia=1,natoms(is)
            if(idxas(ia,is).eq.input%gw%at2)then
              ia2=ia
              is2=is
            endif  
          enddo 
        enddo
        rd(1:3)=atposc(1:3,ia2,is2)-atposc(1:3,ia1,is1)
      endif  
      rdlen=sqrt(rd(1)*rd(1)+rd(2)*rd(2)+rd(3)*rd(3))

!     Create the path grid from at1 -> at2
      np=nrcmt(is1)+nrcmt(is2)+99
      allocate(rs(np))
      do irc=1,nrcmt(is1)
        rs(irc)=rcmt(irc,is1)
      enddo
!
      irlen=rdlen-rmt(is1)-rmt(is2)
      do i=1,99
        rr(i)=rmt(is1)+dble(i)*irlen/1.0d+2
        j=i+nrcmt(is1)
        rs(j)=rr(i)
      enddo
!      
      do irc=1,nrcmt(is2)
        jrc=nrcmt(is2)-irc+1
        krc=irc+99+nrcmt(is1)
        rs(krc)=rdlen-rcmt(jrc,is2)
      enddo
!
!     Loop over eigenvectors
!      
      do ib = input%gw%ibmin, input%gw%ibmax
!
!        open output file
!       
         write(filename,5) ik, ib, input%gw%at1, input%gw%at2
         call str_strip(filename)
         open(unit=71,file=filename,status='unknown')

         allocate(evec(np))
         evec(:)=zzero
         call calcevecplot(ik,ib,ia1,is1,ia2,is2,evec)
!
!        write to file
!
         do i=1,np
           write(71,*) rs(i),real(evec(i)),aimag(evec(i))
         enddo
         deallocate(evec)
         
         close(71)
      end do ! ib

   5  format('evec-',i4,'-',i4,'-',i4,'-',i4,'.out')      
      return
      end subroutine
!EOC      
