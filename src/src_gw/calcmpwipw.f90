!BOP
!
!!ROUTINE: calcmpwipw
!
!!INTERFACE:
!
subroutine calcmpwipw(iq)
!      
! !DESCRIPTION:
!
! This subroutine calculates the matrix elements between PW's and
! orthogonalized IPW's.
!
!!USES:
    use modinput
    use modmain, only : cfunig, zzero, zone, ngrtot, ngrid, igfft
    use modgw,   only : Gset, Gqset, Gqbarc, sgi, sgi_fft, mpwipw, &
    &                   fdebug, time_mpwipw
    use modmpi,  only : rank
    
!!INPUT PARAMETERS:
    implicit none
    integer(4), intent(in) :: iq

!!LOCAL VARIABLES:
    integer(4) :: ig, npw, ngq, ipw, igq
    integer(4), dimension(3) :: igv ! integer coodintates of G-G'
    complex(8), allocatable  :: tmat(:,:)
    real(8) :: tstart, tend
 
!!EXTERNAL ROUTINES: 
    external zgemm

! !REVISION HISTORY:
! 
! Created: May 2006 by RGA
! Revisited 10.05.2011 by DIN

!EOP
!BOC
    call timesec(tstart)

    npw = Gqbarc%ngk(1,iq)
    ngq = Gqset%ngk(1,iq)
        
    if (allocated(mpwipw)) deallocate(mpwipw)
    allocate(mpwipw(ngq,npw))
    mpwipw(:,:) = zzero

    ! Calculate the integral between pw's and IPW's
    allocate(tmat(ngq,npw))
    
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ipw,igq,igv,ig)
!$OMP DO
#endif    
    do ipw = 1, npw
      do igq = 1, ngq
        igv(:) = Gset%ivg(:,Gqbarc%igkig(ipw,1,iq)) - &
        &        Gset%ivg(:,Gqset%igkig(igq,1,iq))
        if ((igv(1).ge.Gset%intgv(1,1)).and.(igv(1).le.Gset%intgv(1,2)).and. &
        &   (igv(2).ge.Gset%intgv(2,1)).and.(igv(2).le.Gset%intgv(2,2)).and. &
        &   (igv(3).ge.Gset%intgv(3,1)).and.(igv(3).le.Gset%intgv(3,2))) then
          ig = Gset%ivgig(igv(1),igv(2),igv(3))
          tmat(igq,ipw) = conjg(cfunig(ig))
        else
          tmat(igq,ipw) = zzero
        end if
      end do ! igq 
    end do ! ipw
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
     
    ! Calculate the overlap PW-IPW integral
    call zgemm( 't','n',ngq,npw,ngq, &
    &           zone,sgi,ngq, &
    &           tmat,ngq, &
    &           zzero,mpwipw,ngq)
    
    ! FFT of S^{*}_{Gi}
    !if (allocated(sgi_fft)) deallocate(sgi_fft)
    !allocate(sgi_fft(ngrtot,ngq))
    !do igq = 1, ngq
    !  sgi_fft(:,igq) = zzero
    !  do ipw = 1, ngq
    !    ig = igfft(Gqset%igkig(ipw,1,iq))
    !    sgi_fft(ig,igq) = sgi(ipw,igq)
    !  end do
    !  call zfftifc(3,ngrid,1,sgi_fft(:,igq))
    !end do

    deallocate(sgi)
        
    if (input%gw%debug) then
      write(fdebug,*) 'CALCMPWIPW, iq = ', iq
      write(fdebug,*) 'CALCMPWIPW, ngq = ', ngq
      write(fdebug,*) 'CALCMPWIPW, npw = ', npw
      write(fdebug,*) "### mpwipw for iq=", iq 
      write(fdebug,*)   
      write(fdebug,*) 'igq  ipw  ivg(ipw)  ivg(igq)  cfunig  mpwipw(jpw,ipw)'
      do ipw = 1, npw, npw/10
        do igq = 1, ngq, ngq/10
          igv(:) = Gset%ivg(:,Gqbarc%igkig(ipw,1,iq))-Gset%ivg(:,Gqset%igkig(igq,1,iq))
          write(fdebug,'(8i5,4f12.6)') igq, ipw, &
          &  Gset%ivg(:,Gqbarc%igkig(ipw,1,iq)), Gset%ivg(:,Gqset%igkig(igq,1,iq)), &
          &  tmat(igq,ipw), mpwipw(igq,ipw)
        end do 
      end do 
    end if ! debug
     
    deallocate(tmat)

    call timesec(tend)
    time_mpwipw = time_mpwipw+tend-tstart

    return
end subroutine
!EOC
