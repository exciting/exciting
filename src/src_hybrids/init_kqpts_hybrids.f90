!
!BOP
! !ROUTINE: init_kqpts_hybrids
! !INTERFACE:
!
subroutine init_kqpts_hybrids
! !USES:
      use modinput
      use modmain
      use modgw
      use modmpi, only: rank
!
! !DESCRIPTION:
!   Initialization for k- and q-point grids.
!EOP
!BOC

      implicit none
      integer(4) :: ik, iq, ig
      integer(4) :: igq, igb
      integer(4) :: i1, i2, i3
      real(8)    :: gpq(3), gqlen
      integer(4) :: ispn
      integer(4) :: iv(3)
      real(8) :: v1(3), v2(3), t1
      real(8), external :: r3taxi
      
      if (allocated(indkp)) deallocate(indkp)
      allocate(indkp(nkptnr))
      ik = 0
      do i3 = 0, input%groundstate%ngridk(3)-1
      do i2 = 0, input%groundstate%ngridk(2)-1
      do i1 = 0, input%groundstate%ngridk(1)-1
         ik = ik+1
         indkp(ik) = ikmap(i1,i2,i3) 
      end do
      end do
      end do

      if(allocated(idikp))deallocate(idikp)
      allocate(idikp(nkpt))
      do ik = 1, nkptnr
        v1(:)=vklnr(:,ik)
        v2(:)=vkl(:,indkp(ik))
        t1 = r3taxi(v1,v2)
        if (t1.lt.input%structure%epslat) then
          idikp(indkp(ik))=ik
        end if
      end do ! ik
      
!     allocate corresponding G+k-vector arrays
      if (allocated(ngknr)) deallocate(ngknr)
      allocate (ngknr(nspnfv,nkptnr))
      if (allocated(igkignr)) deallocate(igkignr)
      allocate (igkignr(ngkmax,nspnfv,nkptnr))
      if (allocated(vgklnr)) deallocate(vgklnr)
      allocate (vgklnr(3,ngkmax,nspnfv,nkptnr))
      if (allocated(vgkcnr)) deallocate(vgkcnr)
      allocate (vgkcnr(3,ngkmax,nspnfv,nkptnr))
      if (allocated(gkcnr)) deallocate(gkcnr)
      allocate (gkcnr(ngkmax,nspnfv,nkptnr))
      if (allocated(tpgkcnr)) deallocate(tpgkcnr)
      allocate (tpgkcnr(2,ngkmax,nspnfv,nkptnr))
      if (allocated(sfacgknr)) deallocate(sfacgknr)
      allocate (sfacgknr(ngkmax,natmtot,nspnfv,nkptnr))
      do ik = 1, nkptnr
        do ispn = 1, nspnfv
          ! generate the G+k-vectors
          call gengpvec(vklnr(:,ik),vkcnr(:,ik),ngknr(ispn,ik),igkignr(:,ispn,ik), &
         &  vgklnr(:,:,ispn,ik),vgkcnr(:,:,ispn,ik),gkcnr(:,ispn,ik),tpgkcnr(:,:,ispn,ik))
          ! generate structure factors for G+k-vectors
          call gensfacgp(ngknr(ispn,ik),vgkcnr(:,:,ispn,ik),ngkmax,sfacgknr(:,:,ispn,ik))
        end do
      end do
      
!-------------------------------------------------!
!     Q-point set and corresponding q+G vectors   !
!-------------------------------------------------!

      ! Determine G-vector cutoff parameters
      gqmax = input%gw%MixBasis%gmb*gkmax
      gmaxbarc = min(input%gw%BareCoul%pwm*gqmax,input%groundstate%gmaxvr)
      
      if (gmaxbarc.gt.input%groundstate%gmaxvr) then
        write(*,*)'WARNING(initgw)! One should increase the value of gmaxvr:'
        write(*,*) 'gkmax=',gkmax,'    gqmax=', gqmax
        write(*,*) 'gmaxvr',input%groundstate%gmaxvr,'    gmaxbarc=', gmaxbarc
      end if
      
      nqptnr = nkptnr
      if (allocated(vql)) deallocate(vql)
      allocate(vql(3,nqptnr))
!     q-grid must not be shifted
      do iq=1,nqptnr
           vql(:,iq) = vklnr(:,iq) - vklnr(:,1)
      end do
      if (allocated(vqc)) deallocate(vqc)
      allocate(vqc(3,nqptnr))
      do iq=1,nqptnr
            vqc(:,iq) = vkcnr(:,iq) - vkcnr(:,1)
      end do
      if (allocated(kqid)) deallocate(kqid)
      allocate(kqid(nkptnr,nkptnr))
      kqid(:,:) = 0
      do ik = 1, nkptnr
      do iq = 1, nkptnr
        iv(:) = ivknr(:,ik)-ivknr(:,iq)
        iv(:) = modulo(iv(1:3),input%groundstate%ngridk(:))
        kqid(ik,iq) = ikmapnr(iv(1),iv(2),iv(3))
      end do
      end do  
 
      if (allocated(ngq)) deallocate(ngq)
      if (allocated(ngbarc)) deallocate(ngbarc)
      allocate(ngq(nqptnr),ngbarc(nqptnr))
      do iq = 1, nqptnr
        ngq(iq) = 0
        ngbarc(iq) = 0
        do ig = 1, ngrtot
          gpq(:) = vgc(:,ig)+vqc(:,iq)
          gqlen = sqrt(gpq(1)*gpq(1)+gpq(2)*gpq(2)+gpq(3)*gpq(3))
          if (gqlen.lt.gqmax) ngq(iq) = ngq(iq)+1
          if (gqlen.lt.gmaxbarc) ngbarc(iq) = ngbarc(iq)+1
        enddo ! ig
      enddo ! iq

      ngqmax=maxval(ngq)
      ngbarcmax=maxval(ngbarc)
      if(ngqmax.eq.ngrtot) write(*,*)'WARNING !! ngqmax = ngrtot !!!'
      if(ngbarcmax.eq.ngrtot) write(*,*)'WARNING !! ngbarcmax = ngrtot !!!'

!     index from G+q to G vector for IPW and coulomb potential
      if(allocated(igqig)) deallocate(igqig)
      allocate(igqig(ngqmax,nqptnr))
      if(allocated(igqigb)) deallocate(igqigb)
      allocate(igqigb(ngbarcmax,nqptnr))
!     G+q vectors
      if(allocated(vgql)) deallocate(vgql)
      allocate(vgql(3,ngqmax,nqptnr))
      if(allocated(vgqc)) deallocate(vgqc)
      allocate(vgqc(3,ngqmax,nqptnr))

      if(allocated(igigq)) deallocate(igigq)
      allocate(igigq(ngrtot,nqptnr))
      igigq(:,:)=0

      if(allocated(igigqb)) deallocate(igigqb)
      allocate(igigqb(ngrtot,nqptnr))
      igigqb(:,:)=0

      igqig(:,:)=0
      igqigb(:,:)=0
      do iq=1,nqptnr
         igq=0
         igb=0
         do ig=1,ngrtot
            gpq(:)=vgc(:,ig)+vqc(:,iq)
            gqlen=sqrt(gpq(1)*gpq(1)+gpq(2)*gpq(2)+gpq(3)*gpq(3))
            if(gqlen.lt.gqmax) then
               igq=igq+1
               vgql(:,igq,iq)=dble(ivg(:,ig))+vql(:,iq)
               vgqc(:,igq,iq)=gpq(:)
               igqig(igq,iq)=ig
            endif
            if(gqlen.lt.gmaxbarc) then
               igb=igb+1
               igqigb(igb,iq)=ig
            endif
         enddo ! ig
         if(igq.ne.ngq(iq))write(*,*)'WARNING(initkqpts): igq.ne.ngqmax'
         if(igb.ne.ngbarc(iq))write(*,*)'WARNING(initkqpts): igb.ne.ngbarcmax'
!        map from iG to iGq vectors
         do igq=1,ngq(iq)
            igigq(igqig(igq,iq),iq)=igq
         end do
!        map from iG to iGq(barc) vectors
         do igq=1,ngbarc(iq)
            igigqb(igqigb(igq,iq),iq)=igq
         end do
      enddo ! iq

      if (.true.) then
         open(99,file='INITKQPTS.OUT',action='write')
         ! Write the list of k-points
         write(99,*) '# Irreducible k-points:'
         call writeklist(nkpt,vkl,vkc)
         write(99,*) '# All k-points:'
         call writeklist(nkptnr,vklnr,vkcnr)
         write(99,*) '# ngrtot, ngkmax, ngqmax, ngbarcmax'
         write(99,*) ngrtot, ngkmax, ngqmax, ngbarcmax
         write(99,*) 
         ! Mapping between irreducible and complete grid k-points
         write(99,*)
         write(99,*) '# INDKP: mapping from the complete k-point set to the corresponding irreducible:'
         write(99,*) 'indkp: ', indkp(:)
         write(99,*)
         write(99,*) '# IDIKP: index of the irreducible k-point in the complete set:'
         write(99,*) 'idikp: ', idikp(:)
         write(99,*)
         close(99)
      end if
      
!-------------------------------------------------!
!     k/q points summary
!-------------------------------------------------!     
      if (rank == 0) then
        call linmsg(fgw,'-','K-points')
        write(fgw,*) 'Total number (NKPTNR) = ', nkptnr
        write(fgw,*) 'Irreducible (NKPT) = ', nkpt
        write(fgw,*) 'Mapping from the complete k-point set to the corresponding irreducible (INDKP):'
        write(fgw,*) indkp(:)
        write(fgw,*) 'Index of the irreducible k-point in the complete set (IDIKP):'
        write(fgw,*) idikp(:)
        write(fgw,*)
      end if

      return
end subroutine
