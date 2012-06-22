!BOP
!
! !ROUTINE: plotevecprod
!
! !INTERFACE:
      subroutine plotevecprod(ik,jk,ib1,ib2,atom1,atom2)

! !DESCRIPTION:
!
! This subroutine calculates the real space representation of the product
! of two eigenvectors in the line joining the two given atoms 
! for ploting. 
!
! !USES:
      use modmain
      use modgw

! !INPUT PARAMETERS:

      implicit none

      integer(4), intent(in) :: ik  ! The k-point for which the (L)APW+lo
!                                     function is ploted
      
      integer(4), intent(in) :: jk  ! The k-point for which the (L)APW+lo
!                                     function is ploted
      
      integer(4), intent(in) :: ib1 ! THe band index of the function

      integer(4), intent(in) :: ib2 ! THe band index of the function

      integer(4), intent(in) :: atom1 ! the atom used as origin

      integer(4), intent(in) :: atom2 ! the atom used as final position
      
      
      
! !LOCAL VARIABLES:

      integer(4) :: ia
      integer(4) :: ia1
      integer(4) :: ia2
      integer(4) :: irc
      integer(4) :: is
      integer(4) :: is1
      integer(4) :: is2
      integer(4) :: jrc
      integer(4) :: krc
      
      real(8)    :: rd(3)
      
      integer(4) :: i,j,np
      
      real(8) :: rdlen,rr(99),irlen
      real(8), allocatable :: rs(:)
      
      complex(8), allocatable :: evecprod(:)
      complex(8), allocatable :: evecprod1(:)
      
      character(len=64) :: filename

!
! !EXTERNAL ROUTINES: 

! !INTRINSIC ROUTINES: 
! 
! !REVISION HISTORY:
!
! Created: 9th. July 2004 by RGA
! Last modified:  June 2006 by RGA
! Revisited 29.04.2011 by DIN
!
!EOP
!BOC

      call boxmsg(6,'-','PLOTEVECPROD')
      
      write(*,*) 'Parameters:'
      write(*,*) '1 k-point number (iik): ', ik
      write(*,*) '2 k-point number (jjk): ', jk
      write(*,*) '1 band index (ib1): ', ib1
      write(*,*) '2 band index (ib2): ', ib2
      write(*,*) 'atom 1 (at1): ', atom1
      write(*,*) 'atom 2 (at2): ', atom2
      write(*,*)

      if((atom1.lt.1).or.(atom1.gt.natmtot)) stop 'atom1 is wrong'
      if((atom2.lt.1).or.(atom2.gt.natmtot)) stop 'atom2 is wrong'
!
!     Set the name of the output file
!
   5  format('evecs-',i4,'-',i4,'-',i4,'-',i4,'-',i4,'-',i4,'.out')      
      write(filename,5) ik, jk, ib1, ib2, atom1, atom2
      call str_strip(filename)
      open(unit=71,file=filename,status='unknown')

   6  format('eprod-',i4,'-',i4,'-',i4,'-',i4,'-',i4,'-',i4,'.out')      
      write(filename,6) ik, jk, ib1, ib2, atom1, atom2
      call str_strip(filename)
      open(unit=72,file=filename,status='unknown')
!
!     Set the indexes of the two atoms
!
      do is = 1,nspecies
        do ia=1,natoms(is)
          if(idxas(ia,is).eq.atom1)then
            ia1=ia
            is1=is
          endif  
        enddo 
      enddo  
      if(atom1.eq.atom2)then
        rd(1:3)=avec(1:3,1)
        is2=is1
        ia2=ia1
      else  
        do is = 1,nspecies
          do ia=1,natoms(is)
            if(idxas(ia,is).eq.atom2)then
              ia2=ia
              is2=is
            endif  
          enddo 
        enddo  
        rd(1:3)=atposc(1:3,ia2,is2)-atposc(1:3,ia1,is1)
      endif  
      rdlen=sqrt(rd(1)*rd(1)+rd(2)*rd(2)+rd(3)*rd(3))
!
!     Prepair the path at1 -> at2 !!!!!
!     
      np=nrcmt(is1)+nrcmt(is2)+99
      allocate(rs(np)) 
      do irc=1,nrcmt(is1)
        rs(irc)=rcmt(irc,is1)
      enddo
      irlen=rdlen-rmt(is1)-rmt(is2)
      do i=1,99
        rr(i)=rmt(is1)+dble(i)*irlen/1.0d+2
      enddo  
      do i=1,99
        j=i+nrcmt(is1)
        rs(j)=rr(i)
      enddo
      do irc=1,nrcmt(is2)
        jrc=nrcmt(is2)-irc+1
        krc=irc+99+nrcmt(is1)
        rs(krc)=rdlen-rcmt(jrc,is2)
      enddo
!
      allocate(evecprod(np))
      allocate(evecprod1(np))
      call calcevecplot(ik,ib1,ia1,is1,ia2,is2,evecprod)
      call calcevecplot(jk,ib2,ia1,is1,ia2,is2,evecprod1)
!      
      do i=1,np
        write(71,'(5g18.10)') rs(i), real(evecprod(i)), aimag(evecprod(i)), &
     &                               real(evecprod1(i)),aimag(evecprod1(i))
        evecprod(i)=evecprod(i)*conjg(evecprod1(i))
        write(72,'(4g18.10)') rs(i), real(evecprod(i)), aimag(evecprod(i)), &
     &      abs(evecprod(i))   
      enddo
      deallocate(evecprod)
      deallocate(evecprod1)
      deallocate(rs)
      close(71)
      close(72)  
      
      return
  
      end subroutine plotevecprod  
!EOC
