!BOP
!
! !ROUTINE: calcpolmat
!
! !INTERFACE:
subroutine calcepsilon(iqp)

! !DESCRIPTION:
!
! This subroutine calculates the dielectric matrix in the RPA approximation
!
! !USES:

      use modmain
      use modgw

! !INPUT PARAMETERS:
      
      implicit none

! !LOCAL VARIABLES:
      integer(4), intent(in) :: iqp

      integer(4) :: ia
      integer(4) :: is
      integer(4) :: ias
      integer(4) :: ie1, ie2, ie12
      integer(4) :: ic   ! Counter: Runs over core states of a given atom.
      integer(4) :: icg  ! Counter: Runs over all core states.
      integer(4) :: iom  ! Counter: Runs over frequencies.
      integer(4) :: ik, jk, ikp, iq, ik0
      integer(4) :: i, isym, lspl, isym0
      integer(4) :: im, jm
      integer(4) :: dimtk
      integer(4) :: recl
      
      real(8)    :: tstart,tend
      real(8)    :: edif, pmn, pvec(3)
      real(8)    :: v1(3),v2(3),v3(3)
      
      complex(8) :: coefw
      complex(8), allocatable :: minm(:,:), micm(:,:)
      complex(8), allocatable :: temp(:,:), body(:,:)
      complex(8), allocatable :: pmat(:,:,:), pmatc(:,:,:)
      complex(8), allocatable :: pm(:), wtmp(:)
      
! !EXTERNAL ROUTINES: 

      external zgemm

! !INTRINSIC ROUTINES: 
      
      intrinsic conjg
      intrinsic cpu_time
        
! !REVISION HISTORY:
!
! Created Dec 2011 by DIN
!
!EOP
!BOC
!
      call cpu_time(tstart)
      if(tstart.lt.0.0d0)write(fgw,*)'warning, tstart < 0'
!                                                                            
!     q-dependent dielectric function
!
      epsilon(:,:,:)=zzero
      
      iq=idikpq(iqp,1)
      coefw=2.0d0*sqrt(pi*vi)

      if (iq.eq.1) then
        
        allocate(pmat(3,nstsv,nstsv))
        recl=16*(3*nstsv*nstsv)
        open(50,file='PMAT.OUT',action='READ',form='UNFORMATTED', &
       &  access='DIRECT',recl=recl)

        if(wcore)then
          allocate(pmatc(3,ncg,nstsv))
          recl=16*(3*ncg*nstsv)
          open(51,file='PMATCOR.OUT',action='READ',form='UNFORMATTED', &
         &   access='DIRECT',recl=recl)
        end if 

!=====================================================================+
!                             HEAD
!=====================================================================+

        if(allocated(head))deallocate(head)
        allocate(head(1:nomeg))
        head=zone
      
        call calchead

!       For calculating wings
        if(allocated(epsw1))deallocate(epsw1)
        allocate(epsw1(mbsiz,nomeg))
        epsw1=zzero

        if(allocated(epsw2))deallocate(epsw2)
        allocate(epsw2(mbsiz,nomeg))
        epsw2=zzero
      
        if(allocated(emac))deallocate(emac)
        allocate(emac(2,nomeg))
        emac=zzero
        
        allocate(wtmp(1:mbsiz))
        
      end if

!======================================================================+
!                             BODY
!======================================================================+
     
      allocate(body(1:mbsiz,1:mbsiz))

!---------------------------------------------------------------------!
!     BZ summation
!---------------------------------------------------------------------!     
      do ikp = 1, nkptq(iqp)

        ik = idikpq(ikp,iqp)
        jk = kqid(ik,iq)
        
        if (ik.eq.1) then
          do iom = 1, nomeg
            do im = 1, mbsiz
              epsilon(im,im,iom)=zone+epsilon(im,im,iom)
            enddo ! im
          end do ! iom
        end if
!        
!       Calculate the matrix elements M^i_{cm} and M^i_{nm}
!
        call expand_prods(ik,iq,0)
!
!       Loop over the symmetry operations
!
        do i = 1, nsymkstar(ikp,iqp)
           
          isym=isymkstar(i,ikp,iqp)

!----------------------------------------------------------------------!
!                       Valence contribution                             !
!----------------------------------------------------------------------!
          dimtk=maxoccband*(nstfv-minunoband+1)
          allocate(minm(1:matsiz,1:dimtk))
          allocate(temp(1:matsiz,1:dimtk))

          ie12=0
          do ie1 = 1, maxoccband
            do ie2 = minunoband, nstfv
              ie12=ie12+1
              minm(1:matsiz,ie12)=minmmat(1:matsiz,ie1,ie2)
            end do
          end do
                     
          ! Get M^i_{nm} by symmetry if needed
          if (isym>1) then
!
!           Matrix representation of the group symmetry operation in MB
!
            call genmbrotmat(iq,isym)
!
!           Rotate M^i_{nm}
!
            call zgemm('c','n',matsiz,dimtk,matsiz, &
           &     zone,rotmat(1:matsiz,1:matsiz),matsiz, &
           &     minm,matsiz,zzero,temp,matsiz)
          else
            temp = minm
          end if
!
!         Transform to the eigenvectors of the coulomb matrix.
! 
          deallocate(minm)
          allocate(minm(1:mbsiz,1:dimtk))

          call zgemm('c','n',mbsiz,dimtk,matsiz, &
         &  zone,barcvm,matsiz,temp,matsiz,zzero,minm,mbsiz)
         
          deallocate(temp)
          
!======================================================================+
!                             WINGS
!======================================================================+
          if (iq.eq.1) then
            
            ik0=indkp(ik)
            if (input%gw%reduceq) then
              lspl=lsplsymc(isym)
            else
              call findkpt(vklnr(:,ik),isym0,ik0)
              lspl=lsplsymc(isym0)
            end if

            read(50,rec=ik0) pmat

            allocate(pm(dimtk))

            ie12=0
            do ie1 = 1, maxoccband
              do ie2 = minunoband, nstfv
                ie12=ie12+1
!               rotate the matrix element from the reduced to non-reduced k-point
!              (note that the inverse operation is used)
                v1(:)=dble(pmat(1:3,ie1,ie2))
                call r3mv(symlatc(:,:,lspl),v1,v2)
                v1(:)=aimag(pmat(1:3,ie1,ie2))
                call r3mv(symlatc(:,:,lspl),v1,v3)
                pvec(:)=cmplx(v2(:),v3(:),8)
                !pmn=sum(pvec(1:3)*q0_eps(1:3))
                pmn=(pvec(1)+pvec(2)+pvec(3))/sqrt(3.0d0)
                edif=evaldft(ie1,ik0)-evaldft(ie2,ik0)
                if(abs(edif).gt.1.0d-10)then
                  pm(ie12)=coefw*pmn/edif
                else
                  write(fgw,*)'WARNING(calcwing): Degenerated CB and VB'
                  pm(ie12)=zzero
                endif
              enddo ! ie2
            end do ! ie1
            
          end if ! iq.eq.1

!---------------------------------------------------------------------!
!         Frequency loop
!---------------------------------------------------------------------!
          allocate(temp(1:mbsiz,1:dimtk))

          do iom = 1, nomeg
            
            ie12=0
            do ie1 = 1, maxoccband
              do ie2 = minunoband, nstfv
                ie12=ie12+1
                temp(1:mbsiz,ie12)=sfact*cmplx(kcw(ie1,ie2,iom,ik),0.d0,8)* &
               &                   minm(1:mbsiz,ie12)
              enddo ! ie2
            enddo ! ie1

            call zgemm('n','c',mbsiz,mbsiz,dimtk, &
           & -zone,temp,mbsiz,minm,mbsiz,zzero,body,mbsiz)
           
            do jm = 1, mbsiz
              do im = 1, mbsiz
                epsilon(im,jm,iom)=epsilon(im,jm,iom)+body(im,jm)
              enddo ! jm
            enddo ! im
            
            ! Wings
            if (iq.eq.1) then 
              
              call zgemv('n',mbsiz,dimtk,-zone,temp,mbsiz,pm,1,zzero,wtmp,1)
              epsw1(:,iom)=epsw1(:,iom)+wtmp(:)
              if (fflg.eq.2) then !! real freq 
                ie12=0
                do ie1 = 1, maxoccband
                  do ie2 = minunoband, nstfv
                    ie12=ie12+1
                    temp(1:mbsiz,ie12)=sfact*conjg(cmplx(kcw(ie1,ie2,iom,ik),0.d0,8))* &
                   &                   minm(1:mbsiz,ie12)
                  enddo ! ie2
                enddo ! ie1
                call zgemv('n',mbsiz,dimtk,-zone,temp,mbsiz,pm,1,zzero,wtmp,1)
              endif ! fflg.eq.2
              epsw2(:,iom)=epsw2(:,iom)+conjg(wtmp(:))

            endif ! iq.eq.1
        
          end do ! iom
          
          deallocate(minm,temp)
          if(iq.eq.1)deallocate(pm)
          
!---------------------------------------------------------------------!
!                       Core contributions                            !
!---------------------------------------------------------------------!
          if (wcore) then
            
            dimtk=ncg*(nstfv-minunoband+1)
            allocate(micm(1:locmatsiz,1:dimtk))
            allocate(temp(1:locmatsiz,1:dimtk))
            
            ie12=0
            do icg = 1, ncg
              do ie2 = minunoband, nstfv
                ie12=ie12+1
                micm(1:locmatsiz,ie12)=micmmat(1:locmatsiz,icg,ie2)
              end do
            end do
            
            ! Get M^i_{nm} by symmetry if needed
            if (isym>1) then
              !
              ! Rotate M^i_{cm}
              !
              call zgemm('c','n',locmatsiz,dimtk,locmatsiz, &
             &   zone,rotmat(1:locmatsiz,1:locmatsiz), &
             &   locmatsiz,micm,locmatsiz,zzero,temp,locmatsiz)
            else
              temp = micm
            end if
!
!           Transform to the eigenvectors of the coulomb matrix
! 
            deallocate(micm)
            allocate(micm(1:mbsiz,1:dimtk))

            call zgemm('c','n',mbsiz,dimtk,locmatsiz, &
           &  zone,barcvm,locmatsiz,temp,locmatsiz,zzero,micm,mbsiz)
           
            deallocate(temp)

!======================================================================+
!                             WINGS
!======================================================================+
            if (iq.eq.1) then

              read(51,rec=ik0) pmatc
              
              allocate(pm(dimtk))

              ie12=0
              do icg = 1, ncg
                is=corind(icg,1)
                ia=corind(icg,2)
                ias=idxas(ia,is)
                ic=corind(icg,3)
                do ie2 = minunoband, nstfv
                  ie12=ie12+1
!                 rotate the matrix element from the reduced to non-reduced k-point
!                 (note that the inverse operation is used)
                  v1(:)=dble(pmatc(1:3,icg,ie2))
                  call r3mv(symlatc(:,:,lspl),v1,v2)
                  v1(:)=aimag(pmatc(1:3,icg,ie2))
                  call r3mv(symlatc(:,:,lspl),v1,v3)
                  pvec(:)=cmplx(v2(:),v3(:),8)
                  !pmn=sum(pvec(1:3)*q0_eps(1:3))
                  pmn=(pvec(1)+pvec(2)+pvec(3))/sqrt(3.0d0)
                  edif=evalcr(ic,ias)-evaldft(ie2,ik0)
                  pm(ie12)=zzero
                  if(abs(edif).gt.1.0d-10)then
                    pm(ie12)=pmn/edif
                  else
                    write(fgw,*)'WARNING(calcepsilon): Strange degeneracy of Core and VB'
                  endif
                enddo ! ie2
              enddo ! icg

            end if ! iq.eq.1
            
            allocate(temp(1:mbsiz,1:dimtk))
            
            do iom = 1, nomeg
              
              ie12=0
              do icg = 1, ncg
                is=corind(icg,1)
                ia=corind(icg,2)
                ias=idxas(ia,is)
                ic=corind(icg,3)
                do ie2 = minunoband, nstfv
                  ie12=ie12+1
                  temp(1:mbsiz,ie12)=sfact*cmplx(unw(ias,ic,ie2,iom,jk),0.d0,8)* &
                 &                   micm(1:mbsiz,ie12)
                enddo ! jst
              enddo ! ist
  
              call zgemm('n','c',mbsiz,mbsiz,dimtk, &
             & -zone,temp,mbsiz,micm,mbsiz,zzero,body,mbsiz)
             
              do jm = 1, mbsiz
                do im = 1, mbsiz
                  epsilon(im,jm,iom)=epsilon(im,jm,iom)+body(im,jm)
                enddo ! jm
              enddo ! im
              
              ! Wings
              if (iq.eq.1) then                                               
                call zgemv('n',mbsiz,dimtk,-zone,temp,mbsiz,pm,1,zzero,wtmp,1) 
                epsw1(:,iom)=epsw1(:,iom)+wtmp(:)                             
                epsw2(:,iom)=epsw2(:,iom)+conjg(wtmp(:))                      
              endif                                                           

            end do ! iom
            
            deallocate(micm,temp)
            if(iq.eq.1)deallocate(pm)
      
          endif ! wcore
          
        end do ! i (symmetry)
        
      end do ! ik
      
      deallocate(body)
      deallocate(minmmat)
      if(wcore)deallocate(micmmat)
      if(allocated(rotmat))deallocate(rotmat)
      if(iq.eq.1)then
        deallocate(pmat)
        close(50)
        if(wcore)then
          deallocate(pmatc)
          close(51)
        end if
        deallocate(wtmp)
      end if ! iq.eq.1
 
      call cpu_time(tend)
      if(tend.lt.0.0d0)write(fgw,*)'warning, tend < 0'
      call write_cputime(fgw,tend-tstart,'CALCEPSILON')
      
      return
end subroutine
!EOC
