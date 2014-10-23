
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: initkqpts
! !INTERFACE:
subroutine init_kqpts
! !USES:
      use modinput
      use modmain
      use modgw
      use modmpi, only: rank
!
! !DESCRIPTION:
!   Generates the ${\bf k}$- and ${\bf q}$-point set and then allocates 
!   and initialises global variables which depend on the ${\bf k}$-point set.
! !LOCAL VARIABLES:

      implicit none
      integer(4) :: ik, iq, ig
      integer(4) :: igq, igb
      real(8)    :: gpq(3), gqlen
      real(8)    :: len
      integer(4) :: ispn
      real(8) :: v1(3), v2(3), t1
      Real (8), External :: r3taxi  
!
! !REVISION HISTORY:
!   Created May 2006 (RGA)
!   Revisited: May (DIN)
!EOP
!BOC
!
!     k-grid shifts
      len=input%groundstate%vkloff(1)**2+ &
     &    input%groundstate%vkloff(2)**2+ &
     &    input%groundstate%vkloff(3)**2
     
      call init1

!---------------------------------------------------!
!     Generate non-reduced k- and q-points meshes   !
!---------------------------------------------------!
      nkptnr=input%groundstate%ngridk(1) * &
     &       input%groundstate%ngridk(2) * &
     &       input%groundstate%ngridk(3)
      ntetnr=6*nkptnr
      nqptnr=nkptnr

      if (allocated(wtetnr)) deallocate(wtetnr)
      allocate(wtetnr(ntetnr))
      if (allocated(tnodesnr)) deallocate(tnodesnr)
      allocate(tnodesnr(4,ntetnr))
      
      If (allocated(ivq)) deallocate (ivq)
      allocate(ivq(3,nqptnr))
      if (allocated(vql)) deallocate(vql)
      allocate(vql(3,nqptnr))
      if (allocated(vqc)) deallocate(vqc)
      allocate(vqc(3,nqptnr))
      if (allocated(linkq)) deallocate(linkq)
      allocate(linkq(ntetnr,nqptnr))
      if (allocated(kqid)) deallocate(kqid)
      allocate(kqid(nqptnr,nqptnr))

!     Generate the k- and q-points meshes
      call kqgen(bvec,input%groundstate%ngridk,ikloff,dkloff,nkptnr, &
        ivk,ivq,dvk,dvq,kqid,ntetnr,tnodesnr,wtetnr,linkq,tvol)

!-------------------------------------------------!
!     K-point set and corresponding k+G vectors   !
!-------------------------------------------------!

!     non-reduced k-points
      if (allocated(vklnr)) deallocate(vklnr)
      allocate (vklnr(3,nkptnr))
      if (allocated(vkcnr)) deallocate(vkcnr)
      allocate (vkcnr(3,nkptnr))

      do ik=1,nkptnr
        vklnr(:,ik)=dble(ivk(:,ik))/dble(dvk)
        call r3mv(bvec,vklnr(:,ik),vkcnr(:,ik))
      end do
      deallocate(ivk)

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
      
!-------------------------------------------------!
!     Q-point set and corresponding q+G vectors   !
!-------------------------------------------------!
!     Determine G-vector cutoff parameters
      gqmax=input%gw%MixBasis%gmb*gkmax
      gmaxbarc=min(input%gw%BareCoul%pwm*gqmax,input%groundstate%gmaxvr)
      
      if(gmaxbarc.gt.input%groundstate%gmaxvr)then
        write(*,*)'WARNING(initgw)! One should increase the value of gmaxvr:'
        write(*,*) 'gkmax=',gkmax,'    gqmax=', gqmax
        write(*,*) 'gmaxvr',input%groundstate%gmaxvr,'    gmaxbarc=', gmaxbarc
      end if

      do iq = 1, nqptnr
        vql(:,iq)= dble(ivq(:,iq))/dble(dvq)
        call r3mv(bvec,vql(:,iq),vqc(:,iq))
      end do

      if (allocated(ngq)) deallocate(ngq)
      if (allocated(ngbarc)) deallocate(ngbarc)
      allocate(ngq(nqptnr),ngbarc(nqptnr))
      do iq = 1, nqptnr
        ngq(iq)=0
        ngbarc(iq)=0
        do ig = 1, ngrtot
          gpq(:)=vgc(:,ig)+vqc(:,iq)
          gqlen=sqrt(gpq(1)*gpq(1)+gpq(2)*gpq(2)+gpq(3)*gpq(3))
          if(gqlen.lt.gqmax) ngq(iq)=ngq(iq)+1
          if(gqlen.lt.gmaxbarc) ngbarc(iq)=ngbarc(iq)+1
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

      if (rank==0) then
         open(99,file='INITKQPTS.OUT',action='write')
!
!        Write the list of k-points
!       
         write(99,*) '# Irreducible k-points:'
         call writeklist(nkpt,vkl,vkc)

         write(99,*) '# All k-points:'
         call writeklist(nkptnr,vklnr,vkcnr)

         write(99,*) '# ngrtot, ngkmax, ngqmax, ngbarcmax'
         write(99,*) ngrtot, ngkmax, ngqmax, ngbarcmax
         write(99,*) 
!
!        Write the list of q-points in WIEN2k format
!      
         write(99,*) "# All q-points:"
         call writeklist(nqptnr,vql,vqc)
!
!        Write the tetrahedra data to file
!
         call writeqgen
      
!        Mapping between irreducible and complete grid k-points
         write(99,*)
         write(99,*) '# INDKP: mapping from the complete k-point set to the corresponding irreducible:'
         write(99,*) 'indkp: ', indkp(:)
         write(99,*)
         write(99,*) '# IDIKP: index of the irreducible k-point in the complete set:'
         write(99,*) 'idikp: ', idikp(:)
         write(99,*)
      end if

!-------------------------------------------------!
!     Generate the small group of q               !
!-------------------------------------------------!

      call gensmallq
      if (rank==0) close(99)
      
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
        call linmsg(fgw,'-','Q-points')
        write(fgw,*) 'Total number (NQPTNR) = ', nqptnr
        write(fgw,*) 'Reduced (NQPT) = ', nqpt
        write(fgw,*) 'Small group of q-vectors:'
        do iq = 1, nqpt
            write(fgw,*) 'iq=', iq, '    nsymq=', nsymq(iq), '    nkptq=', nkptq(iq)
        end do
      end if

      return
end subroutine
!EOC

