!BOP
!!ROUTINE: calchead
!!INTERFACE:
!
subroutine calchead(ik,iomstart,iomend,ndim)
!
!!DESCRIPTION: 
!
! This subroutine calculate the head of the dielectric matrix at the $\Gamma$ point.      
!
!!USES:
    use modinput
    use modmain, only : zzero, zone, pi, evalsv, evalcr, idxas
    use modgw

!!INPUT VARIABLES:
    implicit none
    integer(4), intent(in) :: ik
    integer(4), intent(in) :: iomstart, iomend
    integer(4), intent(in) :: ndim
    
!!LOCAL VARIABLES:
    integer(4) :: ic, icg
    integer(4) :: ia, is, ias
    integer(4) :: ie1, ie2, dimtk
    integer(4) :: ikp, iop, jop
    integer(4) :: iom, fid
    real(8) :: edif    ! energy difference
    real(8) :: edsq    ! edif^2
    real(8) :: tstart, tend
    complex(8) :: coefh
    complex(8) :: pnm, zsum
    
!
!!REVISION HISTORY:
!
! Created 11.02.05 by RGA
! Revisited July 2011 by DIN
!
!EOP
!BOC
    if (input%gw%debug) write(fdebug,*) ' ---- calchead started ----'
    call timesec(tstart)
    
    ! position in the non-reducied grid
    ikp = kset%ik2ikp(ik)     
    
    ! constant prefactor    
    coefh = cmplx(4.d0*pi*vi*occmax,0d0,8)
    
    ! loop over tensor components
    do jop = 1, 3
    do iop = 1, 3
      ! loop over frequencies
      do iom = iomstart, iomend

        !-------------------------
        ! Inter-band contribution
        !-------------------------
        zsum = zzero
        do ie2 = numin, nstdf
          do ie1 = 1, ndim
            if (ie1<=nomax) then
              !==================         
              ! valence-valence
              !==================
              edif = evalsv(ie2,ikp)-evalsv(ie1,ikp)
              if (dabs(edif)>1.d-6) then
                pnm = pmatvv(ie1,ie2,iop)*conjg(pmatvv(ie1,ie2,jop))
                zsum = zsum+fnm(ie1,ie2,iom,ik)*pnm/(edif*edif)
              end if ! edif
            else
              !==================
              ! core-valence
              !==================
              icg = ie1-nomax
              is = corind(icg,1)
              ia = corind(icg,2)
              ias = idxas(ia,is)
              ic = corind(icg,3)
              edif = evalsv(ie2,ikp)-evalcr(ic,ias)
              if (dabs(edif)>1.0d-6) then
                pnm = pmatcv(icg,ie2,iop)*conjg(pmatcv(icg,ie2,jop))
                zsum = zsum+fnm(ie1,ie2,iom,ik)*pnm/(edif*edif)
              end if
            end if
          end do ! ie2
        end do ! ie1
        epsh(iom,iop,jop) = epsh(iom,iop,jop)-coefh*zsum
        
        !-------------------------  
        ! Intra-band contribution
        !-------------------------
        if (metallic) then 
          zsum = zzero
          do ie1 = numin, nomax
            pnm = pmatvv(ie1,ie1,iop)*conjg(pmatvv(ie1,ie1,jop))
            zsum = zsum+kwfer(ie1,ik)*pnm
          enddo ! ie1
          ! for imaginary frequency, a negative sign is needed 
          if (freq%fconv=='imfreq') zsum = -zsum
          epsh(iom,iop,jop) = epsh(iom,iop,jop)- &
          &                   coefh*zsum/(freq%freqs(iom)**2)
        end if ! metallic
        
      end do ! iom
    end do ! iop
    end do ! jop

    ! timing
    call timesec(tend)
    time_dfhead = time_dfhead+tend-tstart
    if (input%gw%debug) write(fdebug,*) ' ---- calchead finished ----'

    return
end subroutine
!EOC
