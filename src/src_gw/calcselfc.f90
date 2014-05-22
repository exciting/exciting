!BOP
!
! !ROUTINE: selfesing
!
! !INTERFACE: 
      subroutine calcselfc(ikp,iqp)
      
! !DESCRIPTION:
!
! This subroutine calculates the singular terms of the selfenergy.
!
! !USES:

      use modmain
      use modgw

! !LOCAL VARIABLES:            

      implicit none
      integer(4), intent(in) :: ikp
      integer(4), intent(in) :: iqp
      
      integer(4) :: ik, iq, jk, jkp
      integer(4) :: ia 
      integer(4) :: is
      integer(4) :: ic       ! (Counter) Runs over core states.
      integer(4) :: ias      ! (Counter) Runs over all atoms
      integer(4) :: ie1      ! (Counter) Runs over bands.
      integer(4) :: ie12
      integer(4) :: ie2      ! (Counter) Runs over bands.
      integer(4) :: iom      ! (Counter) Runs over frequencies.
      integer(4) :: m, m1, m2
      integer(4) :: dimtk
      
      real(8)    :: tstart, tend, tcore
      real(8)    :: enk
      real(8)    :: vi4pi, sqvi4pi
      real(8)    :: coefs1, coefs2
      
      complex(8) :: sum
      complex(8), allocatable :: wmat(:,:), mmat(:,:), wm(:,:)
      complex(8), allocatable :: mwm(:,:,:)  ! Sum_ij{M^i*Wc_ij*conjg(M^j)}
      complex(8), allocatable :: scqval(:,:), scqcor(:,:)
      
! !INTRINSIC ROUTINES: 

      
      intrinsic abs
      intrinsic sign
      intrinsic cmplx
      intrinsic conjg
      intrinsic cpu_time
      

! !EXTERNAL ROUTINES: 
      
      complex(8), external :: zdotu, zdotc
      external zhemm
      complex(8), external :: freqconvl
 
! !REVISION HISTORY:
!
! Created 23.06.05 by RGA.
! Revisted July 2011 by DIN
!
!EOP
!BOC
!
      call cpu_time(tstart)

      vi4pi=4.d0*pi*vi
      sqvi4pi=sqrt(vi4pi)

      coefs1=singc1*sqvi4pi
      coefs2=singc2*vi4pi

      ik=idikp(ikp)
      iq=idikpq(iqp,ikp)
      
      ! k-q point
      jk=kqid(ik,iq)
      jkp=indkp(jk)

!-------------------------------- 
!       Valence contribution
!-------------------------------- 

      dimtk=nbandsgw*nstfv
      allocate(mmat(mbsiz,dimtk))
      ie12=0
      do ie2 = 1, nstfv
        do ie1 = ibgw, nbgw
          ie12=ie12+1
          mmat(1:mbsiz,ie12)=minmmat(1:mbsiz,ie1,ie2)
        end do
      end do
!
!     Sum_ij{M^i*W^c_{ij}*conjg(M^j)}
!
      allocate(mwm(nomeg,ibgw:nbgw,nstfv))
      allocate(wmat(mbsiz,mbsiz))
      allocate(wm(mbsiz,dimtk))
      do iom = 1, nomeg
        wmat(1:mbsiz,1:mbsiz)=inveps(1:mbsiz,1:mbsiz,iom)
        call zhemm('l','u',mbsiz,dimtk, &
        &          zone,wmat,mbsiz,mmat,mbsiz,zzero,wm,mbsiz)
        ie12=0
        do ie2 = 1, nstfv
          do ie1 = ibgw, nbgw
            ie12=ie12+1
            mwm(iom,ie1,ie2) = wkpq(iqp,ikp)*zdotc(mbsiz,mmat(:,ie12),1,wm(:,ie12),1)
            if ((Gamma).and.(ie1.eq.ie2)) then
              mwm(iom,ie1,ie2) = mwm(iom,ie1,ie2) + &
           &    coefs2*head(iom) + &
           &    coefs1*(zdotu(mbsiz,mmat(:,ie12),1,epsw2(:,iom),1) + &
           &            zdotc(mbsiz,mmat(:,ie12),1,epsw1(:,iom),1) )
            end if ! q=\Gamma
          end do
        end do
      end do ! iom
      deallocate(mmat,wmat,wm)
!
!     Frequency convolution
!
      allocate(scqval(ibgw:nbgw,1:nomeg))
      scqval(:,:)=zzero

      do iom = 1, nomeg
        do ie1 = ibgw, nbgw
          m1=max(n12dgn(1,ie1,ikp),ibgw) ! low index
          m2=min(n12dgn(2,ie1,ikp),nbgw) ! upper index
          sum=zzero
          do m = m1, m2
            do ie2 = 1, nstfv
              enk = evaldft(ie2,jkp)-efermi
              sum = sum+freqconvl(iom,nomeg,freqs(iom),enk,mwm(1:nomeg,m,ie2),freqs,womeg)
            end do ! ie2
          end do ! m
          scqval(ie1,iom) = scqval(ie1,iom)+sum/dble(m2-m1+1)
        enddo ! ie1
      enddo ! iom
      deallocate(mwm)

!------------------------------------
!     Core-valence contribution
!------------------------------------

      call cpu_time(tcore)

      if(iopcore.eq.0)then      

        dimtk=nbandsgw*ncg
        allocate(mmat(mbsiz,dimtk))
        ie12=0
        do ie2 = 1, ncg
          do ie1 = ibgw, nbgw
            ie12=ie12+1
            mmat(1:mbsiz,ie12) = mincmat(1:mbsiz,ie1,ie2)
          end do
        end do
!
!       Sum_ij{M^i*W^c_{ij}*conjg(M^j)}
!
        allocate(mwm(nomeg,ibgw:nbgw,ncg))
        allocate(wmat(mbsiz,mbsiz))
        allocate(wm(mbsiz,dimtk))
        do iom = 1, nomeg
          wmat(1:mbsiz,1:mbsiz)=inveps(1:mbsiz,1:mbsiz,iom)
          call zhemm('l','u',mbsiz,dimtk, &
          &          zone,wmat,mbsiz,mmat,mbsiz,zzero,wm,mbsiz)
          ie12=0
          do ie2 = 1, ncg
            do ie1 = ibgw, nbgw
              ie12=ie12+1
              mwm(iom,ie1,ie2) = wkpq(iqp,ikp)*zdotc(mbsiz,mmat(:,ie12),1,wm(:,ie12),1)
            end do ! ie1
          end do ! ie2
        end do ! iom
        deallocate(mmat,wmat,wm)

        allocate(scqcor(ibgw:nbgw,1:nomeg))
        scqcor=zzero
        do iom = 1, nomeg
          do ie1 = ibgw, nbgw
            do m = m1, m2
              sum = zzero
              ! Sum over ie2
              do ie2 = 1, ncg
                is = corind(ie2,1)
                ia = corind(ie2,2)
                ias = idxas(ia,is)
                ic = corind(ie2,3)
                enk = evalcr(ic,ias)-efermi
                sum = sum+freqconvl(iom,nomeg,freqs(iom),enk,mwm(1:nomeg,m,ie2),freqs,womeg)
              end do ! ie2
              scqcor(ie1,iom) = scqcor(ie1,iom)+sum/dble(m2-m1+1)
            end do ! m
          enddo ! ie1
        enddo ! iom
        deallocate(mwm)

      endif ! core
      
!----------------------------------------
!     Sum up the contributions
!----------------------------------------

      write(96,*)'-------- CALCSELFC -------------, ikp = ', ikp, '    iqp = ', iqp
      if(iopcore.eq.0)then
        write(96,*)'# omega      core        valence        selfec'
        do ie1 = ibgw, nbgw
          write(96,*)'band nr. = ', ie1
          do iom = 1, nomeg
            selfec(ie1,ikp,iom)=selfec(ie1,ikp,iom)+scqcor(ie1,iom)+scqval(ie1,iom)
            write(96,10)iom,scqcor(ie1,iom),scqval(ie1,iom),selfec(ie1,ikp,iom)
          end do ! iom
        enddo ! ie1
        write(96,*)
        deallocate(scqcor)
      else 
        write(96,*)'# omega      valence        selfec'
        do ie1 = ibgw, nbgw
          write(96,*)'band nr. = ', ie1
          do iom = 1, nomeg
            selfec(ie1,ikp,iom)=selfec(ie1,ikp,iom)+scqval(ie1,iom)
            write(96,11)iom,scqval(ie1,iom),selfec(ie1,ikp,iom)
          end do ! iom
        enddo ! ie1
        write(96,*)
      endif
      deallocate(scqval)

      call cpu_time(tend)
      if(tend.lt.0.0d0)write(fgw,*)'warning, tend < 0'
      if(tstart.lt.0.0d0)write(fgw,*)'warning, tstart < 0'
      call write_cputime(fgw,tend-tstart, 'CALCSELFC')
      call write_cputime(fgw,tcore-tstart,'CALCSELFC: core')
      call write_cputime(fgw,tend-tcore,  'CALCSELFC: valence')

   10 format(i4,2f18.10,'  | ',2f18.10,' || ',2f18.10)
   11 format(i4,2f18.10,' || ',2f18.10)

      return 
      end subroutine calcselfc
!EOC        
              
                 
                    
                    
                               
                        
