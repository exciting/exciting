!BOP
!
!!ROUTINE: kintw
!
!!INTERFACE:
!
subroutine kintw()
!      
!!DESCRIPTION:
!
! This subroutine calculates the occupancy dependent k-point 
! weights for BZ integrations.
!
!!USES:
    use modinput
    use modmain, only: evalsv, evalcr, nstsv, efermi, &
    &                  natmtot, nspecies, natoms, idxas
    use modgw,   only: kset, kqset, kiw, ciw, kwfer, metallic, fdebug
    use mod_core_states, only: ncmax, ncore
     
!!LOCAL VARIABLES:
    implicit none
    integer(4) :: ia
    integer(4) :: ias
    integer(4) :: is
    integer(4) :: ist
    integer(4) :: ik, ikp
    real(8), allocatable :: bandpar(:,:)
    real(8), allocatable :: cwpar(:,:)
!EOP
!BOC

    if (allocated(kiw)) deallocate(kiw)
    allocate(kiw(nstsv,kqset%nkpt))
    kiw = 0.0d0
    
    if (allocated(kwfer)) deallocate(kwfer)
    allocate(kwfer(nstsv,kqset%nkpt))
    kwfer = 0.0d0
    
    metallic = .false.
    
    !-----------------
    ! Valence states
    !-----------------
    allocate(bandpar(nstsv,kqset%nkpt))
    do ik = 1, kqset%nkpt
      ikp = kset%ik2ikp(ik)
      bandpar(:,ik) = evalsv(:,ikp)
    enddo  
    call tetiw(kqset%nkpt,kqset%ntet,nstsv,bandpar, &
    &          kqset%tnodes,kqset%wtet,kqset%tvol,efermi, &
    &          kiw)
    call tetiwsurf(kqset%nkpt,kqset%ntet,nstsv,bandpar, &
    &              kqset%tnodes,kqset%wtet,kqset%tvol,efermi, &
    &              kwfer)
    deallocate(bandpar)
    
    if (dabs(maxval(kwfer))>1.d-6) metallic = .true.
    
    !-------------
    ! Core states
    !-------------
    if ((input%gw%coreflag=='all').or. &
    &   (input%gw%coreflag=='xal')) then
    
      if (allocated(ciw)) deallocate(ciw)
      allocate(ciw(ncmax,natmtot))
      ciw = 0.0d0
    
      allocate(bandpar(1,kqset%nkpt))
      allocate(cwpar(1,kqset%nkpt))
      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas(ia,is)
          do ist = 1, ncore(is)
            bandpar(1,:) = evalcr(ist,ias) 
            call tetiw(kqset%nkpt,kqset%ntet,1,bandpar, &
            &          kqset%tnodes,kqset%wtet,kqset%tvol,efermi, &
            &          cwpar)
            ciw(ist,ias) = cwpar(1,1)
          end do ! ist
        end do ! ia
      end do ! is
      deallocate(bandpar)
      deallocate(cwpar)

    end if ! core
    
    if (input%gw%debug) then
      call linmsg(fdebug,'-','Info(kintw): BZ integration weights')
      write(fdebug,*) 'VALENCE: kiw'
      do ik = 1, kqset%nkpt
        write(fdebug,*) ik, kiw(:,ik)
      enddo
      if (input%gw%coreflag=='all') then
        write(fdebug,*) 'CORE: ciw'
        do is = 1, nspecies
          do ia = 1, natoms(is)
            ias = idxas(ia,is)
            write(fdebug,*) ias, ciw(:,ias)
          end do
        end do
      end if
    endif
      
    return
end subroutine
!EOC          
         
