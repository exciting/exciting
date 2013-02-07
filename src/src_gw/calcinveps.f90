!BOP
!
! !ROUTINE: calcinveps
!
! !INTERFACE:
    subroutine calcinveps(iqp)

! !DESCRIPTION:
!
! This subroutine calculates the inverse of the dielectric function
!
! !USES:

    use modmain
    use modgw

! !INPUT PARAMETERS:
      
    implicit none
    integer(4), intent(in) :: iqp
    
    integer(4) :: iom
    integer(4) :: im, jm
    
    real(8)    :: tstart, tend

    integer(4) :: info, lwork
    integer(4), allocatable :: ipiv(:)
    complex(8), allocatable :: eps(:,:)
    complex(8), allocatable :: work(:)
    complex(8), allocatable :: w2b(:), bw1(:)
    
    character(len=10) :: sname="calcinveps"    

! !EXTERNAL ROUTINES: 

    external zgetrf, zgetri
    external zhetrf, zhetri
    complex(8), external :: zdotc, zdotu

! !INTRINSIC ROUTINES: 
      
    intrinsic conjg
    intrinsic cpu_time
        
! !REVISION HISTORY:

! Created 31.01.2007 by JH
! Revisited Dec 2011 by DIN
!
!EOP
!BOC

    call cpu_time(tstart)
    
!----------------------------------------------------------------------!
!                    Non-singular q-points                             !
!----------------------------------------------------------------------!

    lwork=64*matsiz
    allocate(ipiv(mbsiz))
    allocate(work(lwork))
    allocate(eps(mbsiz,mbsiz))
    allocate(bw1(mbsiz),w2b(mbsiz))
    
    do iom = 1, nomeg
    
      eps(1:mbsiz,1:mbsiz)=epsilon(1:mbsiz,1:mbsiz,iom)
      
      if(fflg.eq.2) then  !! real freq
        call zgetrf(mbsiz,mbsiz,eps,mbsiz,ipiv,info)
        call errmsg0(info,sname,"calling zgetrf")
        call zgetri(mbsiz,eps,mbsiz,ipiv,work,lwork,info)
        call errmsg0(info,sname,"calling zgetri")
      else   !! imaginary freq
        call zhetrf('u',mbsiz,eps,mbsiz,ipiv,work,lwork,info)
        call errmsg0(info,sname,"calling zhetrf")
        call zhetri('u',mbsiz,eps,mbsiz,ipiv,work,info)
        call errmsg0(info,sname,"calling zhetri")
      endif

!----------------------------------------------------------------------!
!                    Gamma point                                       !
!----------------------------------------------------------------------!

      if (Gamma) then

        if(fflg.eq.2) then 
          call zgemv('n',mbsiz,mbsiz,zone,eps,mbsiz, &
         &            epsw1(:,iom),1,zzero,bw1,1)
          call zgemv('t',mbsiz,mbsiz,zone,eps,mbsiz, &
         &            epsw2(:,iom),1,zzero,w2b,1)
        else 
          call zhemv('u',mbsiz,zone,eps,mbsiz, &
         &            epsw1(:,iom),1,zzero,bw1,1)
          call zhemv('u',mbsiz,zone,eps,mbsiz, &
         &            epsw2(:,iom),1,zzero,w2b,1)
          w2b=conjg(w2b)
        endif 

        emac(2,iom)=head(iom)
        head(iom)=1.d0/(head(iom)-zdotu(mbsiz,epsw2(:,iom),1,bw1,1))
        emac(1,iom)=1.d0/head(iom)

        epsw1(:,iom)= -head(iom)*bw1(:)
        epsw2(:,iom)= -head(iom)*w2b(:)
        do jm=1,mbsiz
          do im=1,mbsiz
            eps(im,jm)=eps(im,jm)+epsw1(im,iom)*epsw2(jm,iom)/head(iom)
          enddo
        enddo


      endif ! iq.eq.Gamma

      inveps(1:mbsiz,1:mbsiz,iom)=eps(1:mbsiz,1:mbsiz)

    enddo ! iom
    
    deallocate(ipiv,work,w2b,bw1)
    deallocate(eps)

    call cpu_time(tend)
    if(tend.lt.0.0d0)write(fgw,*)'warning, tend < 0'
    call write_cputime(fgw,tend-tstart,'CALCINVEPS')

    end subroutine
!EOC
