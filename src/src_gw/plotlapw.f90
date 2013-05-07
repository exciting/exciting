!BOP
!
! !ROUTINE: plotlapw
!
! !INTERFACE:
      subroutine plotlapw

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

      implicit none

! !LOCAL VARIABLES:
      integer(4) :: igp ! The G vector of the function to be ploted
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
      integer(4) :: ispn
      integer(4) :: ikp
      
      character(len=64) :: filename

      real(8) :: irlen
      real(8) :: rd(3),rd2(3)
      real(8) :: rdlen
      real(8) :: rr
      real(8) :: ri(3)
      real(8), allocatable :: rs(:)
      real(8) :: kgvec(3), kgvecl(3)
      real(8) :: phs
      real(8) :: phsat
      
      complex(8) :: pw
      complex(8), allocatable :: yl(:)
      complex(8), allocatable :: apwalm(:,:,:,:,:)
      complex(8), allocatable :: evecmtlm(:,:)
      complex(8), allocatable :: evecmt(:)
      complex(8), allocatable :: evecfv(:)
      complex(8), allocatable :: evec(:)      
      
!
! !EXTERNAL ROUTINES: 
!
      external str_strip
      external boxmsg
      external match
      external ylm
      external wavefmt
      external zgemv
 
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
! Last modified: 16th. Sept 2005 by RGA
! Revisited: 02.05.2011 by DIN
!
!EOP
!BOC

      call boxmsg(6,'-','PLOTLAPW')
      
      ikp=indkp(input%gw%iik)
      write(*,*) 'Parameters:'
      write(*,*) 'k-point number (iik): ', input%gw%iik
      write(*,*) 'irreducible k-point number (ikp): ', ikp
      write(*,*) 'lower bound for k+g-vector number (igmin): ', input%gw%igmin
      write(*,*) 'upper bound for k+g-vector number (igmax): ', input%gw%igmax
      write(*,*) 'atom 1 (at1): ', input%gw%at1
      write(*,*) 'atom 2 (at2): ', input%gw%at2
      write(*,*)

      if((input%gw%at1.lt.1).or.(input%gw%at1.gt.natmtot)) stop 'atom1 is wrong'
      if((input%gw%at2.lt.1).or.(input%gw%at2.gt.natmtot)) stop 'atom2 is wrong'
!
!     Set the indexes of the two atoms
!
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
!
!        Create the path grid from at1 -> at2
!
      np=nrcmt(is1)+nrcmt(is2)+99
      allocate(rs(np))
      do irc=1,nrcmt(is1)
        rs(irc)=rcmt(irc,is1)
      enddo
      irlen=rdlen-rmt(is1)-rmt(is2)
      do i=1,99
        rr=rmt(is1)+dble(i)*irlen/100.0
        j=i+nrcmt(is1)
        rs(j)=rr
      enddo
      do irc=1,nrcmt(is2)
        jrc=nrcmt(is2)-irc+1
        krc=irc+99+nrcmt(is1)
        rs(krc)=rdlen-rcmt(jrc,is2)
      enddo
!
!     Loop over G vectors
!
      do igp = input%gw%igmin, input%gw%igmax
!
!       open output file
!
        write(filename,5) ikp, igp, input%gw%at1, input%gw%at2
        call str_strip(filename)
        open(unit=71,file=filename,status='unknown')
        
!       Allocate the local arrays
        allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))

        allocate(yl(lmmaxapw))
        allocate(evecmtlm(lmmaxapw,nrmtmax))
        allocate(evecmt(nrmtmax))
        allocate(evecfv(nmatmax))
        allocate(evec(np))

!       find the matching coefficients
        do ispn = 1, nspnfv
           call match(ngk(1,ikp),gkc(:,ispn,ikp),tpgkc(:,:,ispn,ikp), &
       &       sfacgk(:,:,ispn,ikp),apwalm(:,:,:,:,ispn))
        end do
       
!       to calculate just a single basis function
        evecfv(:)=0.0d0
        evecfv(igp)=1.0d0

!       calculate the values of the spherical harmonics for atom 1
        call ylm(rd,input%groundstate%lmaxapw,yl)

!       calculate the radial wavefunctions of atom 1
        call wavefmt(input%groundstate%lradstep,input%groundstate%lmaxapw, &
       &   is1,ia1,ngk(1,ikp),apwalm,evecfv,lmmaxapw,evecmtlm)
       
!       convert from spherical harmonics to spherical coordinates
        call zgemv('T',lmmaxapw,nrcmt(is1),zone,evecmtlm,lmmaxapw,yl,1, &
       &   zzero,evecmt,1)
        do irc=1,nrcmt(is1)
           evec(irc)=evecmt(irc)
        enddo
!      
!       Calculate the phase of the plane waves due to the change of origin
!
        kgvec(:)=vgkc(:,igp,1,ikp)
        kgvecl(:)=vgkl(:,igp,1,ikp)
       
        phsat=2.0d0*pi*(kgvecl(1)*atposl(1,ia1,is1)+ &
       &                kgvecl(2)*atposl(2,ia1,is1)+ &
       &                kgvecl(3)*atposl(3,ia1,is1))
        
        irlen=rdlen-rmt(is1)-rmt(is2)
        do i=1,99
          rr=rmt(is1)+dble(i)*irlen/1.0d2
          ri(:)=rd(:)*rr/rdlen
          phs=phsat+(kgvec(1)*ri(1)+kgvec(2)*ri(2)+kgvec(3)*ri(3))
          pw=cmplx(dcos(phs),dsin(phs),8)/sqrt(omega)
          j=i+nrcmt(is1)
          evec(j)=pw
        enddo  
!     
!       calculate the radial wavefunctions of atom 2 (if needed)
!
        if ((is1.eq.is2).and.(ia1.eq.ia2)) then
            do irc=1,nrcmt(is1)
              jrc=nrcmt(is1)-irc+1
              krc=irc+nrcmt(is1)+99
              evec(krc)=evecmt(jrc)
            enddo
        else
            rd2(:)=-1.0d0*rd(:)
            call ylm(rd2,input%groundstate%lmaxapw,yl)
            call wavefmt(input%groundstate%lradstep, &
           &             input%groundstate%lmaxapw,  &
           &             is2,ia2,ngk(1,ikp),apwalm,evecfv,lmmaxapw,evecmtlm)
            call zgemv('t',lmmaxapw,nrcmt(is2),zone,evecmtlm,lmmaxapw,yl,1, &
           &   zzero,evecmt,1)
            do irc=1,nrcmt(is2)
               jrc=nrcmt(is2)-irc+1
               krc=irc+nrcmt(is1)+99
               evec(krc)=evecmt(jrc)
            enddo
        end if
!
!        write to file
!
        do i=1,np
           write(71,*) rs(i),real(evec(i)),aimag(evec(i))
        enddo
        close(71)
        
        deallocate(apwalm)
        deallocate(yl)
        deallocate(evecmtlm)
        deallocate(evecmt)
        deallocate(evecfv)
        deallocate(evec)        

      enddo ! igp    
      
      deallocate(rs)
      return
 
   5  format('lapw-',i4,'-',i4,'-',i4,'-',i4,'.out')
      end subroutine plotlapw  
!EOC
