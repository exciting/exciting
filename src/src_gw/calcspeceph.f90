!==================================================================
! Calculates the q-dependent correlation term of the self-energy
! using the frequency convolution
!==================================================================
subroutine calcspeceph(eval2, evalpath, ik)
    use modinput
    use m_getunit
    use modmain, only : pi, zzero, evalsv, idxas, evalcr, efermi, nstfv 
    use modgw 

    ! input variables
    implicit none
    real(8), intent(in)    :: eval2( nstfv, ngridkqtot)
    real(8), intent(in)    :: evalpath( nstfv, nkpt)
    integer(4), intent(in) :: ik 

    ! local variables            
    integer(4) :: iq, ig
    integer(4) :: ie1, ie2
    integer(4) :: iw, fid
    real(8)    :: ekk, ekq, w, eta, wq, wgq, temp, wqf, wgkq, wgkk
    complex(8) :: weight 
    real(8), external :: wgauss
    !real(8)    :: g2eph (ngridkqtot)
    real(8)    :: w0g1, w0g2, esigmar0

    !!! NOTE: I am assuming here that all inputs are in eV and internal units are in Ha
    ! For now, I consider only the coupling to an Einstein phonon (only one mode) with fixed energy 

    wqf   = 2.d0/ngridkqtot  ! this has to be set to the value of the interpolated grid

    ! this has to be adjusted to a model

    wq    = 0.05   / 27.21139 ! this has be taken from a code with linear response. 
    temp  = 0.0008 / 27.21139 ! 300 K -> 0.025 eV 
    
    ! Bose occupation factor for the phonon 
    wgq = wgauss(-wq/temp,-99) 
    wgq = wgq/(1.d0-2.d0*wgq)

    eta = 0.05 / 27.21139 
    ! 
    if (ik .eq. 1 )then 
      write(6,*) ' eta  = ' , eta  * 27.21139
      write(6,*) ' wq   = ' , wq   * 27.21139
      write(6,*) ' temp = ' , temp  *27.21139
      write(6,*) ' wgq  = ' , wgq 
    endif 
    ! 
    !write(6,102), "qset ", qsetd%vkl(:,1), kqsetd%vkl(:,1),vkl(:,1)
    !write(6,102), "qset ", qsetd%vkl(:,2), kqsetd%vkl(:,2),vkl(:,2)
    !write(6,102), "qset ", qsetd%vkl(:,3), kqsetd%vkl(:,3),vkl(:,3)
    !write(6,102), "qset ", qsetd%vkl(:,4), kqsetd%vkl(:,4),vkl(:,4)
    !write(6,102), "qset ", qsetd%vkl(:,5), kqsetd%vkl(:,5),vkl(:,5)
    !102 format(A5,4x,3f6.3,4x,3f6.3,4x,3f6.3)
    !
    !do iq = 1, 2 ! skip loop over dense mesh for BZ integra
    !
    !kqsetd%kset -> q   grid
    !kqsetd%kset -> k+q grid
    !
    ! construct the coupling coefficients 
    if(ik .eq. 1) then
      if(allocated(g2eph)) deallocate (g2eph) 
      allocate(g2eph(ngridkqtot)) 
      g2eph(:) = 0.d0 
      do iq = 1, ngridkqtot ! loop over dense mesh for BZ integral
        do ig = 1, min(100000,ngvec) 
          if (norm2(kqsetd%kset%vkc(:,iq)+gset%vgc(:,ig))**2 .gt.1e-8)then
            g2eph(iq) = g2eph(iq) + 4.d0 * pi / omega * wq / 2.d0 / norm2(kqsetd%kset%vkc(:,iq)+gset%vgc(:,ig))**2 * 0.1 
            !g2eph(iq) = 1.d0
          endif
        enddo
      enddo
    endif

    !g2eph(:) = 1.d0
    g2eph(:) = 0.1/ omega * 4.d0  
    !

    if(ik .eq. 1) then
      call getunit(fid)
      open(fid,file='wgkq'//ik//'.OUT',action='Write',status='Unknown')
      do iq = 1, ngridkqtot ! loop over dense mesh for BZ integral
        do ie1 = ibeph, nbeph ! loop over states
           ekq = eval2(ie1,iq)-efnew
           wgkq = wgauss(-ekq/temp,-99)
           write(fid,*) iq, ie1, ekq, wgkq 
        enddo 
      enddo 
      close(fid)
    endif 
     
    do iq = 1, ngridkqtot ! loop over dense mesh for BZ integral
      !
      !if (norm2(kqsetd%kset%vkc(:,kqsetd%ikqmt2ik_nr(iq))).gt. 1e-15)then
      !if (norm2(kqsetd%kset%vkc(:,iq)).gt.1e-15)then
        !g2eph = 0.4e-04 !/ norm2(kqsetd%vkc(:,iq)-vkc(:,ik))**2 
        !g2eph = 4.d0 * pi / omega * wq / 2.d0 / norm2(kqsetd%kset%vkc(:,kqsetd%ikqmt2ik_nr(iq)))**2  !* (1./epsinf - 1./eps0) 
        !g2eph = 4.d0 * pi / omega * wq / 2.d0 / norm2(kqsetd%kset%vkc(:,iq))**2 !* (1./epsinf - 1./eps0) 
        !g2eph  = 1.d0

      !g2eph   = 1.35e-04 / norm2(kqsetd%vkc(:,iq)-vkc(:,ik))
      !g2eph = 1.35e-05 / norm2(qsetd%vkl(:,iq))
      !else
      !  g2eph    = 0.d0
      !endif
      
      !g2eph = 1.d0 
      !
      !do ie2 = ibsumeph, nbsumeph ! sum over states
        !
        do ie1 = ibeph, nbeph ! loop over states
          !
          ie2  = ie1 ! assuming that polar coupling is only intraband
          ekq = eval2(ie1,iq)-efnew
          !ekq  = eval2(ie1,kqsetd%ik2ikqmt_nr(iq))-efnew
          !ekq  = eval2(ie1,kqsetd%ikqmt2ik_nr(iq))-efnew

          if (.true. .and. (ik .eq. 2) .and. (iq .lt. 31) .and. (ie1 .eq. 24)) then
            write(6,101) iq, 1./(norm2(kqsetd%kset%vkc(:,iq))**2+1.e-6) ,    &
                             1./(norm2(kqsetd%kset%vkc(:,kqsetd%ikqmt2ik_nr(iq)))**2+1.e-6),    &
                             1./(norm2(kqsetd%kset%vkc(:,kqsetd%ik2ikqmt_nr(iq)))**2+1.e-6),    &
                             eval2(ie1,iq)-efnew,    &
                             eval2(ie1,kqsetd%ikqmt2ik_nr(iq))-efnew,    &
                             eval2(ie1,kqsetd%ik2ikqmt_nr(iq))-efnew   
          endif
          101 format(i4,2x,6f15.6)
          !ekq  = eval2(ie1,kqsetd%ik2ikqmt_nr(iq))-efnew
          wgkq = wgauss(-ekq/temp,-99)
          !
          do iw = 1, nomegeph  ! loop over frequency
            !
            w=freq%freqs(iw)
            !
            weight = wqf *                                                              &
                     ((       wgkq + wgq ) / ( w - ( ekq - wq ) - (0.d0,1.d0) * eta ) + &
                      (1.d0 - wgkq + wgq ) / ( w - ( ekq + wq ) - (0.d0,1.d0) * eta ))
            !
            esigmar0 =  g2eph (iq) *  wqf * real (                                          &
                ( (       wgkq + wgq ) / ( -( ekq - wq ) - (0.d0,1.d0) * eta )  +      &
                  (1.d0 - wgkq + wgq ) / ( -( ekq + wq ) - (0.d0,1.d0) * eta )))
            !
            selfeph(ie1,iw,ik) = selfeph(ie1,iw,ik)+ weight * g2eph(iq) - esigmar0
            !
          end do ! iw
          !
        end do ! ie1
        !
      !end do ! ie2
      !
    end do ! iq
    !
    do ie1 = ibeph, nbeph
      ! 
      ekk    = evalpath(ie1,ik)-efnew
      wgkk   = wgauss(-ekk/temp,-99) 
      !
      do iw = 1, nomegeph       ! loop over frequency 
        !
        w=freq%freqs(iw)
        speceph(ie1,iw,ik) = 2.d0/pi*abs(aimag(selfeph(ie1,iw,ik)))/ &
                      ((w-ekk-real(selfeph(ie1,iw,ik)))**2+(aimag(selfeph(ie1,iw,ik))**2))
        !
        !speceph(ie1,iw,ik) = 2.d0/pi*abs(aimag(selfeph(ie1,iw,ik)))/ ( (w-ekk)**2 + (aimag(selfeph(ie1,iw,ik))**2) )
        !
      enddo
      !
    enddo
    !
    return
    !
end subroutine
