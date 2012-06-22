
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: initkqpts
! !INTERFACE:
subroutine initkqpts
! !USES:
      use modinput
      use modmain
      use modgw
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
      real (8), allocatable :: weight(:)
      integer(4) :: ikvec(3)
      integer(4) :: irkvec(3)
      logical :: found
      
!
! !REVISION HISTORY:
!   Created May 2006 (RGA)
!   Revisited: May (DIN)
!EOP
!BOC
!

!    No k-grid shift is allowed
     len=input%groundstate%vkloff(1)**2+ &
    &    input%groundstate%vkloff(2)**2+ &
    &    input%groundstate%vkloff(3)**2
     if(len.gt.1.0d-10)then
        write(*,*)'WARNING: Not shifted k-grids are only valid for the current version of GW code.'
        write(*,*)'         Please change vkloff option to (0,0,0) and recalculate the eigenvalues and eigenvectors.'
        stop '(initkqpts) vkloff .ne. (0,0,0)'
     end if

!-------------------------------!
!     Irreducible k-mesh data   !
!-------------------------------!
      nirtet = ntet
      allocate(ivkir(3,nkpt))
      allocate(wkir(nkpt))
      allocate(wirtet(nirtet))
      allocate(tndi(4,nirtet))      
      ivkir(1:3,1:nkpt)=ivk(1:3,1:nkpt)
      wkir(1:nkpt)=iwkp(1:nkpt)
      wirtet(1:nirtet)=wtet(1:nirtet)
      tndi(1:4,1:nirtet)=tnodes(1:4,1:nirtet)

!---------------------------------------------------!
!     Generate non-reduced k- and q-points meshes   !
!---------------------------------------------------!
      nqpt=input%groundstate%ngridk(1) * &
     &     input%groundstate%ngridk(2) * &
     &     input%groundstate%ngridk(3)
      ntet=6*nqpt
      
      allocate(ivq(3,nqpt))
      if (allocated(vql)) deallocate(vql)
      allocate(vql(3,nqpt))
      if (allocated(vqc)) deallocate(vqc)
      allocate(vqc(3,nqpt))
      if (allocated(linkq)) deallocate(linkq)
      allocate(linkq(ntet,nqpt))
      if (allocated(kqid)) deallocate(kqid)
      allocate(kqid(nqpt,nqpt))

!     Generate the k- and q-points meshes
      call kqgen(bvec,input%groundstate%ngridk,ikloff,dkloff,nqpt, &
        ivk,ivq,dvk,dvq,kqid,ntet,tnodes,wtet,linkq,tvol)
 
!-------------------------------------------------------------!
!     Non-reduced K-point set and corresponding G+k vectors   !
!-------------------------------------------------------------!
      nkptnr=input%groundstate%ngridk(1) * &
     &       input%groundstate%ngridk(2) * &
     &       input%groundstate%ngridk(3)

!     non-reduced k-points
      if (allocated(vklnr)) deallocate(vklnr)
      allocate (vklnr(3,nkptnr))
      if (allocated(vkcnr)) deallocate(vkcnr)
      allocate (vkcnr(3,nkptnr))

      do ik=1,nkptnr
        vklnr(:,ik)=dble(ivk(:,ik))/dble(dvk)
        call r3mv(bvec,vklnr(:,ik),vkcnr(:,ik))
      end do

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
         ! generate the G+k-vectors
         call gengpvec(vklnr(:,ik),vkcnr(:,ik),ngknr(1,ik),igkignr(:,1,ik), &
        &  vgklnr(:,:,1,ik),vgkcnr(:,:,1,ik),gkcnr(:,1,ik),tpgkcnr(:,:,1,ik))
         ! generate structure factors for G+k-vectors
         call gensfacgp(ngknr(1,ik),vgkcnr(:,:,1,ik),ngkmax,sfacgknr(:,:,1,ik))
      end do

      if(allocated(idikp))deallocate(idikp)
      allocate(idikp(nkpt))
      
      do ik = 1, nkptnr
         ikvec(1:3)=ivk(1:3,ik)
         irkvec(1:3)=ivkir(1:3,indkp(ik))
         found=(ikvec(1).eq.irkvec(1)).and. &
        &      (ikvec(2).eq.irkvec(2)).and. &
        &      (ikvec(3).eq.irkvec(3))  
         if(found)then
            idikp(indkp(ik))=ik
         end if
      end do ! ik
      
!-------------------------------------------------!
!     Q-point set and corresponding G+q vectors   !
!-------------------------------------------------!
!
!     K- and Q-grids are essentially the same    
! 
      nqptnr=nkptnr

      do iq = 1, nqptnr
        vql(:,iq)= dble(ivq(:,iq))/dble(dvq)
        call r3mv(bvec,vql(:,iq),vqc(:,iq))
      end do

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
      if(ngqmax.eq.ngrtot) write(6,*)'WARNING !! ngqmax = ngrtot !!!'
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

!---------------
!     DEBUG
!---------------

      if (debug) then
        
        open(99,file='INITKQPTS.OUT',action='write')
        
        write(99,*) "### ngrtot, ngkmax, ngqmax, ngbarcmax ###"
        write(99,*) ngrtot, ngkmax, ngqmax, ngbarcmax
        write(99,*) 
        write(99,*) "### igqig (indgw) ###"
        do iq=1,nqptnr
          write(99,*) 
          write(99,*) "iq= ",iq
          do ig=1,ngq(iq),ngq(iq)/20
            write(99,'(2i5)') ig, igqig(ig,iq)
          enddo
        enddo
        do iq=1,nqptnr
          write(99,*) 
          write(99,*) "iq= ",iq
          do ig=1,ngbarc(iq),ngbarc(iq)/20
            write(99,'(2i5)') ig, igqigb(ig,iq)
          enddo
        enddo
!
!       Write the list of k-points in WIEN2k format
!         
        write(99,*) "# Reducible k-points "
        call writeklist(nkpt,input%groundstate%ngridk,dvk,wkpt,ivkir)
!
!       Write the list of k-points in WIEN2k format
!         
        allocate(weight(nkptnr))
        weight(:)=1.0d0/dble(nkptnr)
        write(99,*) "# All k-points "
        call writeklist(nkptnr,input%groundstate%ngridk,dvk,weight,ivk)
!
!       Write the list of q-points in WIEN2k format
!        
        write(99,*) "# All q-points "
        call writeklist(nqptnr,input%groundstate%ngridk,dvq,weight,ivq)
        deallocate(weight)
!
!       Write the tetrahedra data to file
!
        call writeqgen
        
!       Symmetry
        write(99,*) "indkp: ", indkp(:)
        write(99,*) "idikp: ", idikp(:)
        
      end if ! DEBUG

!-------------------------------------------------!
!     Generate the small/little group of q        !
!-------------------------------------------------!

      call gensmallq
      
      ! close file INITKQPTS.OUT
      if (debug) close(99) 

      return
end subroutine
!EOC

