
!BOP
! !ROUTINE: gensmallq
! !INTERFACE:

SUBROUTINE gensmallq
! !USES:
    use modinput
    use modmain
    use modgw
!
! !DESCRIPTION:
!
!  Generates the small group of q-vectors, i.e., all R's that $R q = q + G_{R}$.
!  Then the q-dependent BZ(q) is searched.
!  As a last step the subroup of R is determined which regenerate from k=BZ(q) the full BZ.
!
! !LOCAL VARIABLES:

    implicit none
      
    integer :: iqp, ikp  ! q-, k-point indexes for IBZ
    integer :: ik
    integer :: ip, jp, iv(3)
    integer :: isym, lspl, nsym

    real(8) :: s(3,3), v1(3), v2(3), t1
    
    real(8), allocatable :: vklq(:,:)
    integer, allocatable :: scmapq(:)    
    integer, allocatable :: ivwrapq(:,:)
    integer, allocatable :: iwkpq(:)
    
    Real (8), External :: r3taxi    
!
! !REVISION HISTORY:
!
!   Created October 2011 (DIN)
!
!EOP
!BOC

    if (input%gw%reduceq) then
      nqpt=nkpt
      nsym=nsymcrys
    else
      nqpt=nqptnr
      nsym=1
    end if  

    allocate(nsymq(nqpt))
    allocate(nkptq(nqpt))
    allocate(vklq(3,nkptnr))
    allocate(wkpq(nkptnr,nqpt))
    allocate(indkpq(nkptnr,nqpt))
    allocate(iksymq(nkptnr,nqpt))
    allocate(idikpq(nkptnr,nqpt))
    iksymq(:,:)=0

    allocate(nsymkstar(nkptnr,nqpt))
    allocate(isymkstar(nsym,nkptnr,nqpt))
    nsymkstar(:,:)=0
    isymkstar(:,:,:)=0

    allocate(scmapq(nsym))
    allocate(ivwrapq(3,nsym))
    allocate(iwkpq(nkptnr))

!    
!   for each q from IBZ
!
    do iqp = 1, nqpt

!--------------------------------------------------------------------------------
!     find the small group of q
!--------------------------------------------------------------------------------
      if (input%gw%reduceq) then
        v1(:)=vql(:,idikp(iqp))
      else
        v1(:)=vql(:,iqp)
      end if
      
      call findgroupq(.false.,v1,input%structure%epslat, &
     &   bvec,symlat,nsym,lsplsymc,nsymq(iqp),scmapq,ivwrapq)

!--------------------------------------------------------------------------------
!     determine the q-dependent BZ(q) and symmetry operations which regenerate
!     the full (non-reduced) k-grid
!--------------------------------------------------------------------------------      

      ip = 0
      do ik = 1, nkptnr
        v1(:) = vklnr(:,ik)
        ! determine if this point is equivalent to that already in the set
        do isym = 1, nsymq(iqp)
           lspl = lsplsymc(scmapq(isym))
           s(:,:) = dble(symlat(:,:,lspl))
           call r3mtv(s,v1,v2)
           call r3frac(input%structure%epslat,v2,iv)
           do jp = 1, ip
              t1 = r3taxi(vklq(:,jp),v2)
              if (t1.lt.input%structure%epslat) then
                ! equivalent k-point found so add to current weight
                indkpq(ik,iqp) = jp 
                iwkpq(jp) = iwkpq(jp) + 1
                iksymq(ik,iqp) = scmapq(isym)
                goto 10
              end if
           end do !jp
        end do ! isym
        ! add new point to set
        ip = ip + 1
        vklq(:,ip) = v1(:)
        indkpq(ik,iqp) = ip 
        iwkpq(ip) = 1
        iksymq(ik,iqp) = 1
10      continue
      end do ! ik
      nkptq(iqp) = ip ! save the number N(q) of k-points

!     Q-dependent k-point weight      
      do ik = 1, nkptq(iqp)
        wkpq(ik,iqp)=dble(iwkpq(ik))/dble(nkptnr)
      end do
      
!
!     Determine the index of the irreducible point in the non-reduced set
!      
      do ik = 1, nkptnr
        v1(:)=vklnr(:,ik)
        v2(:)=vklq(:,indkpq(ik,iqp))
        t1 = r3taxi(v1,v2)
        if(t1.lt.input%structure%epslat)then
          idikpq(indkpq(ik,iqp),iqp)=ik
        end if
      end do
!
!     Determine the star G_k
!      
      nsymkstar(:,iqp)=0
      do ik = 1, nkptnr
        ikp=indkpq(ik,iqp)
        nsymkstar(ikp,iqp)=nsymkstar(ikp,iqp)+1
        isymkstar(nsymkstar(ikp,iqp),ikp,iqp)=iksymq(ik,iqp)
      end do ! ik
!      
!     Debug information
!
      if (input%gw%reduceq .and. input%gw%debug) then
 
        call boxmsg(99,'-','Generate a small group of q-vectors')
        write(99,*)
        write(99,*) 'iqp = ', iqp
        write(99,*)
        write(99,*) '    nsymq: ', nsymq(iqp)
        write(99,*) '    nkptq: ', nkptq(iqp)
        write(99,*) '    vklq: '
        do ikp = 1, nkptq(iqp)
          write(99,'(i4,4f12.4)') ikp, vklq(:,ikp), wkpq(ikp,iqp)
        end do
        write(99,*) '    idikpq: ', idikpq(1:nkptq(iqp),iqp)      
        write(99,*) '    indkpq: ', indkpq(:,iqp)
        write(99,*) '    iksymq: ', iksymq(:,iqp)
        do ikp = 1, nkptq(iqp)
          write(99,*) 'ikp =', ikp
          write(99,*) 'nsymkstar =', nsymkstar(ikp,iqp)
          write(99,*) 'isymkstar =', isymkstar(1:nsymkstar(ikp,iqp),ikp,iqp)
        enddo
      
!       write(99,*)'iqp=',iqp
!       do ik = 1, nkptnr
!         ikp=indkpq(ik,iqp)
!         isym=iksymq(ik,iqp)
!         lspl=lsplsymc(isym)
!         ilspl=isymlat(lspl)
!         s(:,:)=dble(symlat(:,:,ilspl))
!         call r3mtv(s,vklq(:,ikp),v1)
!         !call r3frac(input%structure%epslat,v1,iv)
!         write(99,*) 'ikp=',ikp,'ik=',ik,'isym=',isym 
!         write(99,*) 'vklq                     R.vklq'
!         write(99,'(6f12.4)') vklq(:,ikp) ,v1
!       end do
!       write(99,*)
      
      end if ! debug
      
    end do ! iqp

    deallocate(iwkpq)
    deallocate(vklq)
    deallocate(scmapq)
    deallocate(ivwrapq)
      
END SUBROUTINE gensmallq
!EOC
