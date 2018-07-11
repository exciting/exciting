!==================================================================
! Calculates the q-dependent correlation term of the self-energy
! using the frequency convolution
!==================================================================
subroutine calcselfeph(eval2, evalpath, ik)
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
    !
    !g2eph(:) = 1.d0
    g2eph(:) = ( 0.05 )**2  / omega * 4.d0  
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
      do ie2 = ibsumeph, nbsumeph ! sum over states
        !
        do ie1 = ibeph, nbeph ! loop over states
          !
          !ie2  = ie1 ! assuming that polar coupling is only intraband
          !ekq = eval2(ie1,iq)-efnew
          ekq  = eval2(ie2,kqsetd%ik2ikqmt_nr(iq))-efnew
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
          !do iw = 1, nomegeph  ! loop over frequency
          !
          ekk    = evalpath(ie1,ik)-efnew
          !
          weight = wqf *                                                              &
                   ((       wgkq + wgq ) / ( ekk - ( ekq - wq ) - (0.d0,1.d0) * eta ) + &
                    (1.d0 - wgkq + wgq ) / ( ekk - ( ekq + wq ) - (0.d0,1.d0) * eta ))
          !
          !esigmar0 =  g2eph *  wqf * real (                                          &
          !    ( (       wgkq + wgq ) / ( -( ekq - wq ) - (0.d0,1.d0) * eta )  +      &
          !      (1.d0 - wgkq + wgq ) / ( -( ekq + wq ) - (0.d0,1.d0) * eta )))
          !
          selfeph0(ie1,ik) = selfeph0(ie1,ik)+ weight * g2eph(iq) 
          !
          !end do ! iw
          !
        end do ! ie1
        !
      end do ! ie2
      !
    end do ! iq
    !
!    do ie1 = ibeph, nbeph
!      ! 
!      ekk    = evalpath(ie1,ik)-efnew
!      wgkk   = wgauss(-ekk/temp,-99) 
!      !
!      do iw = 1, nomegeph       ! loop over frequency 
!        !
!        w=freq%freqs(iw)
!        speceph(ie1,iw,ik) = 2.d0/pi*abs(aimag(selfeph(ie1,iw,ik)))/ &
!                      ((w-ekk-real(selfeph(ie1,iw,ik)))**2+(aimag(selfeph(ie1,iw,ik))**2))
!        !
!        !speceph(ie1,iw,ik) = 2.d0/pi*abs(aimag(selfeph(ie1,iw,ik)))/ ( (w-ekk)**2 + (aimag(selfeph(ie1,iw,ik))**2) )
!        !
!      enddo
!      !
!    enddo
    !
    return
    !
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FC: Auxiliary functions, for now taken from Quantum Espresso 
!-----------------------------------------------------------------------
function wgauss (x, n)
  !-----------------------------------------------------------------------
  !
  !     this function computes the approximate theta function for the
  !     given order n, at the point x.
  !
  ! --> (n>=0) : Methfessel-Paxton case. See PRB 40, 3616 (1989).
  !
  ! --> (n=-1 ): Cold smearing (Marzari-Vanderbilt). See PRL 82, 3296 (1999)
  !       1/2*erf(x-1/sqrt(2)) + 1/sqrt(2*pi)*exp(-(x-1/sqrt(2))**2) + 1/2
  !
  ! --> (n=-99): Fermi-Dirac case: 1.0/(1.0+exp(-x)).
  !
  !
  use modmain, only : pi 

  implicit none
  real(8) :: wgauss, x
  ! output: the value of the function
  ! input: the argument of the function
  integer :: n
  ! input: the order of the function
  !     
  !    the local variables
  !

  real(8) :: a, hp, arg, hd, xp
  ! the coefficient a_n
  ! the hermitean function
  ! the argument of the exponential
  ! the hermitean function
  ! auxiliary variable (cold smearing)
  integer :: i, ni
  ! counter on the n indices
  ! counter on 2n
  real(8), external :: gauss_freq, qe_erf
  real(8), parameter :: maxarg = 200.d0
  ! maximum value for the argument of the exponential

  ! Fermi-Dirac smearing
  if (n.eq. - 99) then
     if (x.lt. - maxarg) then
        wgauss = 0.d0
     elseif (x.gt.maxarg) then
        wgauss = 1.d0
     else
        wgauss = 1.0d0 / (1.0d0 + exp ( - x) )
     endif
     return

  endif
  ! Cold smearing
  if (n.eq. - 1) then
     xp = x - 1.0d0 / sqrt (2.0d0)
     arg = min (maxarg, xp**2)
     wgauss = 0.5d0 * qe_erf (xp) + 1.0d0 / sqrt (2.0d0 * pi) * exp ( - &
          arg) + 0.5d0
     return

  endif
  ! Methfessel-Paxton
  wgauss = gauss_freq (x * sqrt (2.0d0) )
  if (n.eq.0) return
  hd = 0.d0
  arg = min (maxarg, x**2)
  hp = exp ( - arg)
  ni = 0
  a = 1.d0 / sqrt (pi)
  do i = 1, n
     hd = 2.0d0 * x * hp - 2.0d0 * DBLE (ni) * hd
     ni = ni + 1
     a = - a / (DBLE (i) * 4.0d0)
     wgauss = wgauss - a * hd
     hp = 2.0d0 * x * hd-2.0d0 * DBLE (ni) * hp
     ni = ni + 1
  enddo
  return
end function wgauss
!
!
!
!---------------------------------------------------------------------
function qe_erf (x)  
  !---------------------------------------------------------------------
  !
  !     Error function - computed from the rational approximations of
  !     W. J. Cody, Math. Comp. 22 (1969), pages 631-637.
  !
  !     for abs(x) le 0.47 erf is calculated directly
  !     for abs(x) gt 0.47 erf is calculated via erf(x)=1-erfc(x)
  !
  implicit none  
  real(8), intent(in) :: x
  real(8) :: x2, p1 (4), q1 (4)
  real(8), external :: qe_erfc  
  real(8) :: qe_erf
  data p1 / 2.426679552305318E2, 2.197926161829415E1, &
            6.996383488619136,  -3.560984370181538E-2 /
  data q1 / 2.150588758698612E2, 9.116490540451490E1, &
            1.508279763040779E1, 1.000000000000000 /
  !
  if (abs (x) > 6.d0) then  
     !
     !  erf(6)=1-10^(-17) cannot be distinguished from 1
     !
     qe_erf = sign (1.d8, x)  
  else  
     if (abs (x)  <= 0.47d0) then  
        x2 = x**2  
        qe_erf=x *(p1 (1) + x2 * (p1 (2) + x2 * (p1 (3) + x2 * p1 (4) ) ) ) &
                / (q1 (1) + x2 * (q1 (2) + x2 * (q1 (3) + x2 * q1 (4) ) ) )
     else  
        qe_erf = 1.d0 - qe_erfc (x)  
     endif
  endif
  !
  return  
end function qe_erf
!
!
!---------------------------------------------------------------------
function gauss_freq (x)
  !---------------------------------------------------------------------
  !
  !     gauss_freq(x) = (1+erf(x/sqrt(2)))/2 = erfc(-x/sqrt(2))/2
  !             - See comments in erf
  !
  implicit none
  real(8),intent(in) :: x
  real(8)            :: gauss_freq
  real(8), parameter :: c =  0.7071067811865475
  !        ( c= sqrt(1/2) )
  real(8), external :: qe_erfc
  !
  gauss_freq = 0.5 * qe_erfc ( - x * c)
  !
  return
end function gauss_freq
!
!
function qe_erfc (x)  
  !---------------------------------------------------------------------
  !   
  !     erfc(x) = 1-erf(x)  - See comments in erf
  !
  implicit none  
  real(8),intent(in) :: x
  real(8)            :: qe_erfc
  real(8) :: ax, x2, xm2, p2 (8), q2 (8), p3 (5), q3 (5), pim1
  real(8), external :: qe_erf  
  data p2 / 3.004592610201616E2,  4.519189537118719E2, &
            3.393208167343437E2,  1.529892850469404E2, &
            4.316222722205674E1,  7.211758250883094,   &
            5.641955174789740E-1,-1.368648573827167E-7 /
  data q2 / 3.004592609569833E2,  7.909509253278980E2, &
            9.313540948506096E2,  6.389802644656312E2, &
            2.775854447439876E2,  7.700015293522947E1, &
            1.278272731962942E1,  1.000000000000000 /
  data p3 /-2.996107077035422E-3,-4.947309106232507E-2, &
           -2.269565935396869E-1,-2.786613086096478E-1, &
           -2.231924597341847E-2 /
  data q3 / 1.062092305284679E-2, 1.913089261078298E-1, &
            1.051675107067932,    1.987332018171353,    &
            1.000000000000000 /

  data pim1 / 0.56418958354775629 /  
  !        ( pim1= sqrt(1/pi) )
  ax = abs (x)  
  if (ax > 26.0) then  
     !
     !  erfc(26.0)=10^(-296); erfc( 9.0)=10^(-37);
     !
     qe_erfc = 0.0  
  elseif (ax > 4.0) then  
     x2 = x**2  
     xm2 = (1.0 / ax) **2  
     qe_erfc = (1.0 / ax) * exp ( - x2) * (pim1 + xm2 * (p3 (1) &
          + xm2 * (p3 (2) + xm2 * (p3 (3) + xm2 * (p3 (4) + xm2 * p3 (5) &
          ) ) ) ) / (q3 (1) + xm2 * (q3 (2) + xm2 * (q3 (3) + xm2 * &
          (q3 (4) + xm2 * q3 (5) ) ) ) ) )
  elseif (ax > 0.47) then  
     x2 = x**2  
     qe_erfc = exp ( - x2) * (p2 (1) + ax * (p2 (2) + ax * (p2 (3) &
          + ax * (p2 (4) + ax * (p2 (5) + ax * (p2 (6) + ax * (p2 (7) &
          + ax * p2 (8) ) ) ) ) ) ) ) / (q2 (1) + ax * (q2 (2) + ax * &
          (q2 (3) + ax * (q2 (4) + ax * (q2 (5) + ax * (q2 (6) + ax * &
          (q2 (7) + ax * q2 (8) ) ) ) ) ) ) )
  else  
     qe_erfc = 1.0 - qe_erf (ax)  
  endif
  !
  ! erf(-x)=-erf(x)  =>  erfc(-x) = 2-erfc(x)
  !
  if (x < 0.0) qe_erfc = 2.0 - qe_erfc  
  !
  return  
end function qe_erfc

