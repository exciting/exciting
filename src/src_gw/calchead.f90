!BOP
!
! !ROUTINE: calchead
!
! !INTERFACE:
      subroutine calchead

!
! !DESCRIPTION: 
!
!This subroutine calculate the dielectric matrix at the $\Gamma$ point.      
!

! !USES:
      use modmain
      use modgw

! !LOCAL VARIABLES:

      implicit none
      
      integer(4) :: icg
      integer(4) :: ia
      integer(4) :: is
      integer(4) :: ic   ! Counter, runs over core states
      integer(4) :: ias  ! Counter, runs over equivalent atoms
      integer(4) :: ie1  ! Counter, runs over bands.
      integer(4) :: ie2  ! Counter, runs over bands.
      integer(4) :: ikp, ik  ! Counter, runs over kpoints
      integer(4) :: iom  ! Counter, runs over frequencies
      integer(4) :: dimtk
      
      real(8) :: cpivi    
      real(8) :: edif    ! energy differences
      real(8) :: edsq    ! edif^2
      real(8) :: tstart, tend
      
      complex(8) :: coef
      complex(8) :: psq, pvec(3)
      complex(8) :: sumfs
      complex(8), allocatable :: fnm(:), pnm(:)
      complex(8), allocatable :: pmat(:,:,:)
      complex(8), allocatable :: pmatc(:,:,:)
!
! !EXTERNAL ROUTINES: 
!
      complex(8), external :: zdotu
!
! !INTRINSIC ROUTINES: 
!
      intrinsic cpu_time
!
! !REVISION HISTORY:
!
! Created 11.02.05 by RGA
! Revisited July 2011 by DIN
!
!EOP
!BOC
      call cpu_time(tstart)

      allocate(pmat(3,nstsv,nstsv))
      if(iopcore.eq.0)then
        allocate(pmatc(3,ncg,nstsv))
      end if 
!
!     Local arrays
!
      dimtk=(nstfv-numin+1)
      allocate(fnm(numin:nstfv))
      allocate(pnm(numin:nstfv))
!
!     Loop over k-points
!
      cpivi=4.0d0*pi*vi
      do ikp = 1, nkpt
      
        ik=idikp(ikp)
        coef=zone*iwkp(ikp)*sfact*cpivi
        
!---------------------------------          
!       Valence-valence
!---------------------------------                    

        read(50,rec=ikp) pmat

        do ie1 = 1, nomax
          do ie2 = numin, nstfv
            pvec(1:3)=pmat(1:3,ie1,ie2)
            psq=sum(pvec(1:3)*q0eps(1:3))
            psq=psq*conjg(psq)
            !psq=sum(pvec(1:3)*conjg(pvec(1:3)))/3.0d0
            edif=evaldft(ie2,ikp)-evaldft(ie1,ikp)
            edsq=edif*edif
            if(edsq.gt.1.0d-10)then
              pnm(ie2)=psq/edsq
            else
              write(fgw,*)'WARNING(calchead): Degenerated CB and VB'
              pnm(ie2)=zzero  
            endif
          enddo ! ie2
          do iom = 1, nomeg
            fnm(numin:nstfv)=cmplx(kcw(ie1,numin:nstfv,iom,ik),0.0d0,8)
            head(iom)=head(iom)-coef*zdotu(dimtk,fnm,1,pnm,1)
          end do
        enddo ! ie1
!
!       the part is related to metallic systems
!       sumfs := sum of |P_nnk|^2 \delta(e_F - e_nk) over n and k on Fermi surface
!
        if(metallic) then 
          sumfs=zzero
          do ie1 = numin, nomax
            pvec(1:3)=pmat(1:3,ie1,ie1)
            psq=sum(pvec(1:3)*q0eps(1:3))
            psq=psq*conjg(psq)
            !psq=sum(pvec(1:3)*conjg(pvec(1:3)))/3.0d0
            sumfs=sumfs+kwfer(ie1,ik)*psq
          enddo ! ie1
          !! for imaginary frequency, a negative sign is needed 
          if(fflg.eq.3) sumfs=-sumfs
          do iom = 1, nomeg
            head(iom)=head(iom)-coef*(1.d0/freqs(iom)**2)*sumfs
          enddo 
        endif 
        
!---------------------------------
!         core-valence
!---------------------------------
        if(iopcore.eq.0)then
 
          read(51,rec=ikp) pmatc

          do icg = 1, ncg
            is=corind(icg,1)
            ia=corind(icg,2)
            ias=idxas(ia,is)
            ic=corind(icg,3)
            do ie2 = numin, nstfv
              pvec(1:3)=pmatc(1:3,icg,ie2)
              psq=sum(pvec(1:3)*q0eps(1:3))
              psq=psq*conjg(psq)
              !psq=sum(pvec(1:3)*conjg(pvec(1:3)))/3.0d0
              edif=evaldft(ie2,ikp)-evalcr(ic,ias)
              edsq=edif*edif
              if(edsq.gt.1.0d-10)then
                pnm(ie2)=psq/edsq
              else
                write(fgw,*)'WARNING(calchead): Degenerated Core and VB!'
                pnm(ie2)=zzero
              endif
            enddo ! ie2
            do iom = 1, nomeg
              fnm(numin:nstfv)=cmplx(unw(ias,ic,numin:nstfv,iom,ik),0.0d0,8)
              head(iom)=head(iom)-coef*zdotu(dimtk,fnm,1,pnm,1)
            end do
          enddo ! icg
       
        endif ! iopcore.eq.0

      enddo ! ikp

      deallocate(fnm,pnm)
      deallocate(pmat)
      if(iopcore.eq.0)then
        deallocate(pmatc)
      end if
!      
!     Write head to file
!       
      open(43,file='HEAD.OUT',action='WRITE',form='FORMATTED')
      do iom=1,nomeg
        write(43,1)iom,freqs(iom),head(iom)
      enddo
      close(43)
    1 format(i4,3f16.8)

      call cpu_time(tend)
      if(tend.lt.0.0d0)write(fgw,*)'warning, tend < 0'
      call write_cputime(fgw,tend-tstart,'CALCHEAD')

      return
      end subroutine calchead
!EOC      
