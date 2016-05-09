!BOP
!
!!ROUTINE: setuprod
!
!!INTERFACE:
!
subroutine setuprod(ia,is)
!       
!!DESCRIPTION:
!
! This subroutine calculates the radial product functions
! and their overlap matrix
!
!!USES:
    use modinput
    use modmain
    use modgw
    use reallocate
    use mod_mpi_gw
    
!!INPUT VARIABLES:
    implicit none
    integer(4), intent(in) :: ia
    integer(4), intent(in) :: is

!!LOCAL VARIABLES:
    integer(4) :: ias
    integer(4) :: io1, io2
    integer(4) :: ilo1, ilo2
    integer(4) :: ipr1, ipr2
    integer(4) :: ir, l1, l2, ist
    integer(4) :: nupcore
    real(8) :: fr(nrmtmax)
    real(8) :: gr(nrmtmax) 
    real(8) :: cf(3,nrmtmax)

!!REVISION HISTORY:
!
! Created 17. May 2006 by RGA
! Revisited 5.05.2011 by DIN
! Modified Dec 2013 by DIN
!
!EOP
!BOC
    ! global atomic index
    ias = idxas(ia,is)

    if (input%gw%debug) then
      write(fdebug,*) 'ncmax,lmaxapw,nlomax:', ncmax, input%groundstate%lmaxapw, nlomax
      write(fdebug,*) 'ias, ncore:', ncore(is)
      write(fdebug,*) 'ias, nlorb:', nlorb(is)
      write(fdebug,*) 'lmaxmb:', input%gw%MixBasis%lmaxmb
      write(fdebug,*) 'maxnup:', maxnup
    end if

    ! radial product functions
    if (allocated(uprod)) deallocate(uprod) 
    allocate(uprod(nrmtmax,maxnup))
    uprod(:,:) = 0.0d0
    
    ! l1 and l2 combination which form L
    if (allocated(eles)) deallocate(eles)
    allocate(eles(maxnup,2))
    eles(:,:) = 0
    
    nup = 0
    nupcore = 0

    !==============
    ! Core states
    !==============
    if (input%gw%coreflag.ne.'vab') then 

      do ist = 1, ncore(is)
        l1 = spl(ist,is)
        if (l1 <= input%gw%MixBasis%lmaxmb) then

          !------------------------------------
          ! Product of core and apw functions
          !------------------------------------
          do l2 = 0, input%gw%MixBasis%lmaxmb
            do io2 = 1, apword(l2,is)
              if (apwdm(io2,l2,is) == 0) then
                nup = nup+1
                nupcore = nupcore+1
                eles(nup,1) = l1
                eles(nup,2) = l2
                do ir = 1, nrmt(is)
                  uprod(ir,nup) = ucore(ir,1,ist,ias)*    &
                  &               apwfr(ir,1,io2,l2,ias)* &
                  &               spr(ir,is)
                end do ! ir
              end if ! apwdm=0
            end do ! io2
          end do ! l2
              
          !----------------------------------------------
          ! Product of core and local orbital functions
          !----------------------------------------------
          do ilo2 = 1, nlorb(is)
            l2 = lorbl(ilo2,is)
            if (l2 <= input%gw%MixBasis%lmaxmb) then
              nup = nup+1
              nupcore = nupcore+1
              eles(nup,1) = l1
              eles(nup,2) = l2 
              do ir = 1, nrmt(is)
                uprod(ir,nup) = ucore(ir,1,ist,ias)* &
                &               lofr(ir,1,ilo2,ias)* &
                &               spr(ir,is)
              end do ! ir
            end if ! l2
          end do! ilo2
              
        end if ! l1
      end do ! ist
          
    end if ! coreflag.ne.vab
    
    !================
    ! Valence states
    !================ 
    do l1 = 0, input%gw%MixBasis%lmaxmb
      do io1 = 1, apword(l1,is)
        if (apwdm(io1,l1,is)==0) then
            
          !-----------------------------------------------
          ! products between valence / conduction states
          !-----------------------------------------------
          do l2 = l1, input%gw%MixBasis%lmaxmb
            do io2 = 1, apword(l2,is)
              if (apwdm(io2,l2,is) == 0) then
                nup = nup+1
                eles(nup,1) = l1
                eles(nup,2) = l2 
                do ir = 1, nrmt(is)
                  uprod(ir,nup) = apwfr(ir,1,io1,l1,ias)* &
                  &               apwfr(ir,1,io2,l2,ias)* &
                  &               spr(ir,is)
                end do ! ir
              end if ! apwdm=0
            end do ! io2
          end do ! l2

          !-------------------------------------------
          ! product between valence / local orbitals
          !-------------------------------------------
          do ilo2 = 1, nlorb(is)
            l2 = lorbl(ilo2,is)
            if (l2 <= input%gw%MixBasis%lmaxmb) then
              nup = nup+1
              eles(nup,1) = l1
              eles(nup,2) = l2 
              do ir = 1, nrmt(is)
                uprod(ir,nup) = apwfr(ir,1,io1,l1,ias)* &
                &               lofr(ir,1,ilo2,ias)*    &
                &               spr(ir,is)
              end do ! ir
            end if ! l2
          end do ! ilo2

        end if ! apwdm=0
      end do ! io1
    end do ! l1

    do ilo1 = 1, nlorb(is)
      l1 = lorbl(ilo1,is)
      if (l1 <= input%gw%MixBasis%lmaxmb) then

        !--------------------------------------------         
        ! products between local and local orbitals
        !--------------------------------------------
        do ilo2 = ilo1, nlorb(is)
          l2 = lorbl(ilo2,is)
          if (l2 <= input%gw%MixBasis%lmaxmb) then
            nup = nup+1
            eles(nup,1) = l1
            eles(nup,2) = l2 
            do ir = 1, nrmt(is)
              uprod(ir,nup) = lofr(ir,1,ilo1,ias)* &
              &               lofr(ir,1,ilo2,ias)* &
              &               spr(ir,is)
            end do ! ir
          end if ! l2
        end do ! ilo2

      end if ! l1
    end do ! ilo1

    if (input%gw%debug) then
      write(fdebug,*) 'setuprod: # of uprod for atom', ias, nup
      write(fdebug,*) '          # from core states ', nupcore
    end if 
    
    !=============================================
    ! allocate and initialize the overlap matrix
    !=============================================
    if (allocated(umat)) deallocate(umat) 
    allocate(umat(maxnup,maxnup))
    umat = 0.0d0
    
    do ipr1 = 1, nup
      do ipr2 = ipr1, nup
        do ir = 1, nrmt(is)
          fr(ir) = uprod(ir,ipr1)*uprod(ir,ipr2)
        end do
        call fderiv(-1,nrmt(is),spr(:,is),fr,gr,cf)
        umat(ipr1,ipr2) = gr(nrmt(is))
      end do ! ipr2
    end do ! nup
     
    if (input%gw%debug) then
      write(fdebug,*) 
      write(fdebug,*) "###umat for atom ", ias
      do ipr1 = 1, nup
        do ipr2 = ipr1, nup
          write(fdebug,'(2i4,f20.6)') ipr1, ipr2, umat(ipr1,ipr2)
        end do ! ipr2
      end do ! ipr1
      write(fdebug,*)        
    end if ! debug
      
    return
end subroutine
!EOC
