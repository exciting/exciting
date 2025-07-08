!BOP
!
!!ROUTINE: \verb"expand_evec"
!
!!INTERFACE:
!
    subroutine expand_evec(ik,trans)
!
!!DESCRIPTION:
!
!Calculate the product of an eigenvector with the corresponding matching coefficients
!
!!USES:
    use mod_bands, only: eveckalm, eveckpalm, eveck, eveckp
    use modinput, only: input
    use modmain, only : ngkmax, apwordmax, lmmaxapw, natmtot, &
    &                   nspecies, natoms, idxas, idxlm, apword, nstfv, &
    &                   lsplsymc, isymlat, symlatc
    use modgw,   only : kqset, Gkqset

!!INPUT PARAMETERS:
    implicit none
    integer(4),   intent(in) :: ik    ! index of the k-point      
    character(1), intent(in) :: trans
      
!!LOCAL VARIABLES:
    integer(4) :: ia, is, ias
    integer(4) :: io
    integer(4) :: ist
    integer(4) :: l, m, lm, m1
    integer(4) :: igk, ngk
    integer(4) :: isym, lspl, ilspl, first_state, last_state
    complex(8), allocatable :: apwalm(:,:,:,:)
    
 
!!EXTERNAL ROUTINES: 
    complex(8), external :: zdotc
    complex(8), external :: zdotu
    
!!REVISION HISTORY:
! 
! Created Dic. 2003 by RGA
! Last Modification: July. 20th. 2004 by RGA
! Revision: 5.05.2011 by DIN
!
!EOP  
!BOC
    ! find matching coefficients
    allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
    
    ngk = Gkqset%ngk(1,ik)
    call match(ngk, &
    &          Gkqset%gkc(:,1,ik), &
    &          Gkqset%tpgkc(:,:,1,ik), &
    &          Gkqset%sfacgk(:,:,1,ik),&
    &          apwalm)

    select case (trans)
    
    case ('t','T')
      first_state = lbound( eveck, 2 )
      last_state = ubound( eveck, 2 )
      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas(ia,is)
          do l = 0, input%groundstate%lmaxapw
            do m = -l, l
              lm = idxlm(l,m)
              do io = 1, apword(l,is)
                do ist = first_state, last_state
                  eveckalm(ist,io,lm,ias) = &
                  &  zdotu(ngk,  &
                  &        eveck(1:ngk,ist),1, &
                  &        apwalm(1:ngk,io,lm,ias),1)
                enddo ! ist
              enddo ! io
            enddo ! m1     
          enddo ! l1       
        enddo  !ia
      enddo !is
      
    case ('c','C')
      first_state = lbound( eveckp, 2 )
      last_state = ubound( eveckp, 2 )
      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas(ia,is)
          do l = 0, input%groundstate%lmaxapw
            do m = -l, l
              lm = idxlm(l,m)
              do io = 1, apword(l,is)
                do ist = first_state, last_state
                  eveckpalm(ist,io,lm,ias) = &
                  &  zdotc(ngk, &
                  &        apwalm(1:ngk,io,lm,ias),1, &
                  &        eveckp(1:ngk,ist),1)
                enddo ! ist  
              enddo ! io
            enddo ! m1     
          enddo ! l1       
        enddo  !ia
      enddo !is
      
    case default
      write(*,*)'ERROR in expand_evec'
      write(*,*)'Allowed values of trans are t, T, c or C'
      write(*,'(a19,a1)')'Received trans = ',trans 
      stop 
    end select
    
end subroutine
!EOC      
