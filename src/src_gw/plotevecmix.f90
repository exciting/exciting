!BOP
!
! !ROUTINE: plotevecmix
!
! !INTERFACE:
      subroutine plotevecmix(iqp,ikp,jkp,ib1,ib2,atom1,atom2)

! !DESCRIPTION:
!
! This subroutine calculates the real space representation of the product
! of two eigenvectors using the expansion in mixed basis functions:
!
!\begin{equation}
!\Psi_{n\vec{k}}(\vec{r})\Psi^*_{m\vec{k}-\vec{q}}(\vec{r})=\sum\limits_i
!M^i_{nm}(\vec{k},\vec{q})\tilde{\chi}_i^{\vec{q}}(\vec{r})
!\end{equation}
!In the line joining the two given atoms 
! for ploting. 
!
! !USES:

      use modmain
      use modgw

! !INPUT PARAMETERS:

      implicit none
      integer(4), intent(in) :: iqp ! q-point index
      integer(4), intent(in) :: ikp ! The k-point for which the (L)APW+lo
!                                     function is ploted
      integer(4), intent(in) :: jkp ! The k-point for which the (L)APW+lo
!                                     function is ploted
      integer(4), intent(in) :: ib1 ! THe band index of the function
      integer(4), intent(in) :: ib2 ! THe band index of the function
      integer(4), intent(in) :: atom1 ! the atom used as origin
      integer(4), intent(in) :: atom2 ! the atom used as final position
      
! !LOCAL VARIABLES:

      integer(4) :: i
      integer(4) :: ia
      integer(4) :: ia1
      integer(4) :: ia2
      integer(4) :: ias1
      integer(4) :: ias2
      integer(4) :: igp
      integer(4) :: im
      integer(4) :: imix
      integer(4) :: imr
      integer(4) :: ir
      integer(4) :: irc      
      integer(4) :: is
      integer(4) :: is1
      integer(4) :: is2
      integer(4) :: j
      integer(4) :: jgp
      integer(4) :: jr
      integer(4) :: jrc
      integer(4) :: krc
      integer(4) :: l
      integer(4) :: lm
      integer(4) :: m
      integer(4) :: np

      real(8) :: irlen
      real(8) :: kgvec(3)
      real(8) :: kgvecl(3)
      real(8) :: qvec(3)
      real(8) :: phs
      real(8) :: phsat
      real(8) :: rd(3)
      real(8) :: rd2(3)
      real(8) :: rdlen
      real(8) :: ri(3)
      real(8) :: rr(99)
      real(8) :: t1
      real(8), allocatable :: rs(:)
      
      complex(8) :: eph
      complex(8), allocatable :: evecmix(:)
      complex(8) :: pw(99)
      complex(8), allocatable :: yl(:)

      character(len=64) :: filename

! !EXTERNAL ROUTINES: 

! !INTRINSIC ROUTINES: 

      intrinsic achar
      intrinsic cos
      intrinsic sin
      intrinsic mod

! !REVISION HISTORY:
!
! Created: 9th. July 2004 by RGA
! Last modified: 9th. June 2006 by RGA
! Revisited: May 2011 by DIN
!
!EOP
!BOC
!BOC
      call boxmsg(6,'-','PLOTEVECMIX')
      
      write(*,*) 'Parameters:'
      write(*,*) '1 k-point number (iik): ', ikp
      write(*,*) '2 k-point number (jjk): ', jkp
      write(*,*) '1 band index (ib1): ', ib1
      write(*,*) '2 band index (ib2): ', ib2
      write(*,*) 'atom 1 (at1): ', atom1
      write(*,*) 'atom 2 (at2): ', atom2
      write(*,*)

      if((atom1.lt.1).or.(atom1.gt.natmtot))stop 'atom1 is wrong'
      if((atom2.lt.1).or.(atom2.gt.natmtot))stop 'atom2 is wrong'
!
!     Set the name of the output file
!
   5  format('evmix-',i4,'-',i4,'-',i4,'-',i4,'-',i4,'-',i4,'.out')
      write(filename,5) ikp, jkp, ib1, ib2, atom1, atom2
      call str_strip(filename)
!
!     open output files
!     
      open(unit=72,file=filename,status='unknown')
!
!     Set the indexes for two atoms
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

      write(*,'(a)')'# plotevecmix: atposl,rd,rlen'
      write(*,'(a,3f12.4)')'# atposl 1:', atposl(:,ia1,is1)
      write(*,'(a,3f12.4)')'# atposl 1:', atposl(:,ia2,is2)
      write(*,'(a,3f12.4)')'# rd:', rd
      write(*,'(a,f12.4)')'# rdlen:', rdlen

      ias1=idxas(ia1,is1)
      ias2=idxas(ia2,is2)

      np=nrmt(is1)+nrmt(is2)+99
      allocate(evecmix(np))
      allocate(rs(np))
      allocate(yl((maxbigl+1)*(maxbigl+1)))
      
      evecmix(1:np)=zzero
      do irc=1,nrmt(is1)
        rs(irc)=spr(irc,is1)
      enddo

!     calculate the values of the spherical harmonics for atom 1
      call ylm(rd,maxbigl,yl)

!     Calculate the phase of the plane waves due to the change of origin
      qvec(1:3)=vql(1:3,iqp)
      phsat=2.0d0*pi*(qvec(1)*atposl(1,ia1,is1)+ &
                      qvec(2)*atposl(2,ia1,is1)+ &
                      qvec(3)*atposl(3,ia1,is1))
      eph=cmplx(cos(phsat),sin(phsat),8)

      im=0
      do imr=1,nmix(ias1)
        l=bigl(ias1,imr)  
        do m=-l,l
          im=im+1
          imix=locmixind(ias1,im)
          lm=l*l+l+m+1
          do ir=1,nrmt(is1)
            t1=1.0d0/spr(ir,is1)
            evecmix(ir)=evecmix(ir)+minmmat(imix,ib1,ib2)*yl(lm)*      &
     &                  cmplx(umix(ias1,imr,ir)*t1,0.0d0,8)*eph
          enddo
        enddo
      enddo
!      
      irlen=rdlen-rmt(is1)-rmt(is2)
      do i=1,99
        rr(i)=rmt(is1)+dble(i)*irlen/1.0d+2
      enddo  

      pw(1:99)=zzero
      do igp=1,ngq(iqp)
        imix=locmatsiz+igp
        do jgp=1,ngq(iqp)
          kgvec(1:3)=vgqc(1:3,jgp,iqp)
          kgvecl(1:3)=vgql(1:3,jgp,iqp)
          phsat=2.0d0*pi*(kgvecl(1)*atposl(1,ia1,is1)+ &
         &                kgvecl(2)*atposl(2,ia1,is1)+ &
         &                kgvecl(3)*atposl(3,ia1,is1))
          do i=1,99
            ri(1:3)=rd(1:3)*rr(i)/rdlen
            phs=phsat+(kgvec(1)*ri(1)+kgvec(2)*ri(2)+kgvec(3)*ri(3))
            eph=cmplx(cos(phs),sin(phs),8)/sqrt(omega)
            pw(i)=pw(i)+minmmat(imix,ib1,ib2)*conjg(sgi(jgp,igp))*eph
          enddo  
        enddo  
      enddo  
      do i=1,99
        j=i+nrmt(is1)
        rs(j)=rr(i)
        evecmix(j)=pw(i)
      enddo

!     calculate the values of the spherical harmonics for atom 2
      rd2(1:3)=-1.0d0*rd(1:3)
      call ylm(rd2,maxbigl,yl)
      do irc=1,nrmt(is2)
        jrc=nrmt(is2)-irc+1
        krc=irc+99+nrmt(is1)
        rs(krc)=rdlen-spr(jrc,is2)
      enddo

      phsat=2.0d0*pi*(qvec(1)*atposl(1,ia2,is2)+ &
                      qvec(2)*atposl(2,ia2,is2)+ &
                      qvec(3)*atposl(3,ia2,is2))
      eph=cmplx(cos(phsat),sin(phsat),8)

      im=0
      do imr=1,nmix(ias2)
        l=bigl(ias2,imr)  
        do m=-l,l
          im=im+1
          imix=locmixind(ias2,im)
          lm=l*l+l+m+1
          do ir=1,nrmt(is2)
            jr=nrmt(is2)-ir+1
            krc=ir+99+nrmt(is1)
            t1=1.0d0/spr(jr,is2)
            evecmix(krc)=evecmix(krc)+minmmat(imix,ib1,ib2)*yl(lm)*      &
     &                   cmplx(umix(ias2,imr,jr)*t1,0.0d0,8)*eph
          enddo
        enddo
      enddo

      do i=1,np
        write(72,'(4g18.10)')rs(i),real(evecmix(i)),aimag(evecmix(i)),&
     &    abs(evecmix(i))   
      enddo

      deallocate(evecmix)
      deallocate(rs)
      deallocate(yl)
      
      close(71)
      close(72)
      
      return
      end subroutine plotevecmix
!EOC
