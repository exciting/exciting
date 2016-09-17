
module strconst

    real(8), allocatable :: rstr(:,:) ! The lattice vectors used for the
                                      ! calculation of the structure constants
end module strconst

!BOP
!
!!ROUTINE: sigma
!
subroutine sigma(iq,lambdamax)
!
!!DESCRIPTION:
!
!This subroutine calculates the lattice sums $\Sigma^{a,a'}_{\lambda,\mu}(\vec{q})=\sum_{\vec{R}}{\frac{e^{i\vec{q}\cdot\left(\vec{R}+\vec{r}_{aa'}\right)}}%
!{|\vec{R}+\vec{r}_{aa'}|^{(\lambda+1)}}Y_{\lambda\mu}\left(\hat{R}_{aa'}\right)}$ using the method
!described in appendix \ref{ewaldmeth}, in particular equation \ref{strconstdef}.
!
!!USES:
    use modinput
    use modmain, only : nspecies, natoms, idxas, atposc, &
    &                   zzero, bvec, pi, natmtot, zi
    use modgw,   only : kqset, fdebug
    use mod_coulomb_potential, only : sgm
    use mod_misc_gw, only : avec, vi    
    use strconst
      
!!INPUT PARAMETERS:
    implicit none
    integer(4), intent(in) :: iq        ! index of the q-point for which 
                                        ! sigma is calculated
    integer(4), intent(in) :: lambdamax ! Maximum value of lambda
     
!!LOCAL VARIABLES:
    integer(4) :: np ! Number of points for the real space
                     ! summation
    integer(4) :: ng ! Number of points for the reciprocal
                     ! space summation
    integer(4) :: ia, ias, is, ja, jas, js
    integer(4) :: i1
    integer(4) :: lmbd
    integer(4) :: mu
    integer(4) :: lmuind
      
    real(8) :: eta
    real(8) :: rcf
    real(8) :: gcf
    real(8) :: rleng               ! the length of rpaa
    real(8) :: qtraa               ! the scalar product qvec.rpaa
    real(8) :: gausr               ! value of the gaussian function e^(-(rleng/eta)^2)
    real(8) :: pref                ! prefactor for the reciprocal lattice sum = 4 pi^(3/2)/v
    real(8) :: gleng               ! the length of gqv
    real(8) :: gqtra               ! the scalar product gqv.raa
    real(8) :: gausg               ! value of the gaussian function e^(-(gleng*eta/2)^2)
    real(8) :: gtolam              ! value of gleng^(lambda-2)/2^(lambda-0.5)
    real(8) :: gamlam              ! gamlam = Gamma[lambda+1/2]
    real(8) :: erfr
    real(8) :: gammaor             

    real(8), dimension(3) :: qvec  ! cartesian coords. of the q point
    real(8), dimension(3) :: raa   ! Vector going from atom 1 to atom 2.
    real(8), dimension(3) :: rpaa  ! the corresponding sum R+raa
    real(8), dimension(3) :: qtemp
    real(8), dimension(3) :: g     ! vector belonging to the reciprocal space lattice
    real(8), dimension(3) :: gqv   ! the corresponding sum G+qvec

    complex(8) :: ylam((lambdamax+1)*(lambdamax+1)) ! the values of the spherical harmonics
    complex(8) :: stmp1((lambdamax+1)*(lambdamax+1)) ! temporary allocation of the values of sigma
    complex(8) :: stmp2((lambdamax+1)*(lambdamax+1)) ! temporary allocation of the values of sigma
    complex(8) :: expqdr           ! e^(qvec.rpaa)
    complex(8) :: itolam           ! i^lambda
    complex(8) :: term1
      
! !EXTERNAL ROUTINES: 
    real(8), external :: derfc
    real(8), external :: calceta
    real(8), external :: gcutoff
    real(8), external :: rcutoff

! !REVISION HISTORY:
!
! Created 21. Jan. 2004 by RGA
! Last Modified 3. Nov 2005 by RGA
! Revisited: May 2011 by DIN
!
!EOP
!BOC

    if (allocated(sgm)) deallocate(sgm)
    allocate(sgm(natmtot,natmtot,(lambdamax+1)*(lambdamax+1)))
    sgm(:,:,:) = zzero
!
! Lattice sums cutoff parameters
!
    eta = calceta()
    rcf = 2.0d+0*rcutoff(input%gw%barecoul%stctol,eta,10)
    gcf = 2.0d+0*gcutoff(input%gw%barecoul%stctol,eta,10)
    qvec(1:3) = -1.0d0*kqset%vqc(1:3,iq)

    do is = 1, nspecies
      do ia = 1, natoms(is)
        ias = idxas(ia,is)
        do js = 1, nspecies
          do ja = 1, natoms(js)
            jas = idxas(ja,js)
            
            ! Initialize the temporal storage of the lattice sums
            stmp1 = zzero
            stmp2 = zzero
            
            !----------------------------------------------
            ! Calculate all the R's such that R+r_aa < rcf
            !----------------------------------------------
            raa(1:3) = atposc(1:3,ia,is)-atposc(1:3,ja,js)
            call genrstr(rcf,raa,avec,np)
            
            !---------------------------------
            ! Sum over the real space lattice
            !---------------------------------
            do i1 = 1, np
              rpaa(1:3) = rstr(1:3,i1)
              rleng = rstr(4,i1)
              ! calculate the values of the spherical harmonics at rpaa
              call ylm(rpaa,lambdamax,ylam)
              qtraa = qvec(1)*rpaa(1)+qvec(2)*rpaa(2)+qvec(3)*rpaa(3)
              expqdr = cmplx(dcos(qtraa),dsin(qtraa),8)
              gausr = dexp(-1.0d0*rleng*rleng/(eta*eta))
              erfr = derfc(rleng/eta)
              gammaor = dsqrt(pi)*erfr/rleng
              term1 = expqdr*cmplx(gammaor,0.0d0,8)*ylam(1)
              stmp1(1) = stmp1(1)+term1
              do lmbd = 1, lambdamax
                gammaor = (dble(lmbd)-0.5d0)*gammaor/rleng+ &
                &          rleng**(lmbd-2)*gausr/(eta**(2*lmbd-1))
                do mu = -lmbd, lmbd
                  lmuind = lmbd*lmbd+lmbd+mu+1
                  stmp1(lmuind) = stmp1(lmuind)+ &
                  &               cmplx(gammaor,0.0d0,8)*expqdr*ylam(lmuind)
                enddo ! mu
              enddo ! lmbd
            enddo ! i1
            deallocate(rstr)
            
            !----------------------------------------------
            ! Calculate all the G's such that G+q < gcf
            !----------------------------------------------
            qtemp(1:3) = -1.0d0*qvec(1:3)
            call genrstr(gcf,qtemp,bvec,ng)
            
            !---------------------------------------
            ! Sum over the reciprocal space lattice
            !---------------------------------------
            pref = 4.0d0*pi*dsqrt(pi)*vi
            do i1 = 1, ng
              gqv(1:3) = rstr(1:3,i1)
              g(1:3) = gqv(1:3)-qtemp(1:3)
              gleng = rstr(4,i1)
              ! Calculate the values of the spherical harmonics at rpaa
              call ylm(gqv,lambdamax,ylam)
              gqtra = g(1)*raa(1)+g(2)*raa(2)+g(3)*raa(3)
              expqdr = cmplx(dcos(gqtra),dsin(gqtra),8)
              gausg = dexp(-2.5d-1*eta*gleng*eta*gleng)
              gtolam = 1.0d0/(gleng*gleng)
              term1 = cmplx(pref*gtolam*gausg,0.0d0,8)*expqdr*ylam(1)
              stmp2(1) = stmp2(1) + term1
              itolam = cmplx(1.0d0,0.0d0,8)
              do lmbd = 1, lambdamax
                gtolam = -1.0d0*gleng*gtolam/2.0d0
                itolam = itolam*zi
                do mu = -lmbd, lmbd
                  lmuind = lmbd*lmbd+lmbd+mu+1
                  stmp2(lmuind) = stmp2(lmuind)+ &
                  &               cmplx(pref*gtolam*gausg,0.0d0,8)* &
                  &               expqdr*ylam(lmuind)*itolam
                enddo ! mu
              enddo ! lmbd
            enddo !i1
            deallocate(rstr)
              
            gamlam = dsqrt(pi)
            stmp1(1) = stmp1(1)*cmplx(1.0d0/gamlam,0.0d0,8)
            stmp2(1) = stmp2(1)*cmplx(1.0d0/gamlam,0.0d0,8)
            sgm(ias,jas,1) = stmp1(1)+stmp2(1)

            do lmbd = 1, lambdamax
              gamlam = 5.0d-1*dble(2*lmbd-1)*gamlam
              do mu = -lmbd, lmbd
                lmuind = lmbd*lmbd+lmbd+mu+1
                stmp1(lmuind) = stmp1(lmuind)*cmplx(1.0d0/gamlam,0.0d0,8)
                stmp2(lmuind) = stmp2(lmuind)*cmplx(1.0d0/gamlam,0.0d0,8)
                sgm(ias,jas,lmuind) = stmp1(lmuind)+stmp2(lmuind)
              enddo
            enddo
            
            if (ias==jas) sgm(ias,jas,1) = sgm(ias,jas,1)-1.0d0/(eta*pi)
            
          enddo ! ja
        enddo ! js
      enddo ! ia
    enddo ! is
!_______________________________________________________________________________    
!    
    if (input%gw%debug) then
      write(fdebug,*)'### eta, rcutoff, gcutoff', eta, rcf, gcf
      write(fdebug,*)'### qvec', qvec
      write(fdebug,*)'-----------------------------------------------------'
      write(fdebug,*)'                    sigma: iq=', iq 
      write(fdebug,*)'-----------------------------------------------------'
      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas(ia,is)
          do js = 1, nspecies
            do ja =1, natoms(js)
              jas = idxas(ja,js)
              do lmbd = 0, lambdamax
                do mu = -lmbd, lmbd
                  lmuind = lmbd*lmbd+lmbd+mu+1
                  write(fdebug,'(4(i3,1x),2e16.8)') &
                  &  ias,jas,lmbd,mu,sgm(ias,jas,lmuind)
                enddo ! mu
              enddo ! lmbd
            enddo ! ja
          enddo ! js
        enddo ! ia    
      enddo ! is
      write(fdebug,*)'-----------------------------------------------------'
      write(fdebug,*)'                    sigma: END'
      write(fdebug,*)'-----------------------------------------------------'
      write(fdebug,*)
    end if
   
    return
end subroutine sigma
!EOC

!===============================================================================
!
!===============================================================================

!BOP
!
!!ROUTINE: genrstr
!
!!INTERFACE:
subroutine genrstr(rmax,rshift,rbas,nr)
!      
!!DESCRIPTION:
!
! Generates the indexes of the lattice vectors to be included in the
! calculation of the structure constants, under the condition:
!
!\begin{equation}
!|\vec{R}+\vec{r}_{aa'}|\le R_{cutoff}
!\end{equation}
!      
! !USES:      
    use strconst
      
!!INPUT PARAMETERS:      
    implicit none
    real(8), intent(in) :: rmax      ! Maximum radius
    real(8), intent(in) :: rshift(3) ! Shift of the origin
    real(8), intent(in) :: rbas(3,3) ! Bravais lattice basis

!!OUTPUT PARAMETERS:
    integer(4), intent(out) :: nr    ! number of vectors
      
!!LOCAL VARIABLES:
    integer(4) :: i, imax, ir, rdim, i1, i2, i3, gap
    real(8) :: lrmin                 ! minimum length of the basis vectors.
    real(8) :: rleng
    real(8), dimension(3) :: r       ! vector belonging to the real space lattice
    real(8), dimension(3) :: lrbs    ! length of the basis vectors
    real(8), dimension(3) :: rps
    real(8), dimension(4) :: rtmp
    logical :: done

!!REVISION HISTORY:
!
! Created 5. Aug. 2004 by RGA
! Revisited Oct. 2013 by DIN
!
!EOP
!BOC
         
    do i = 1, 3
      lrbs(i) = dsqrt(rbas(1,i)*rbas(1,i)+ &
     &                rbas(2,i)*rbas(2,i)+ &
     &                rbas(3,i)*rbas(3,i))
    enddo
    lrmin = minval(lrbs)
    imax = idint(rmax/lrmin)+1
    rdim = (2*imax+1)*(2*imax+1)*(2*imax+1)
      
    allocate(rstr(4,rdim))
      
    ir = 0
    do i1 = -imax,imax
      do i2 = -imax,imax
        do i3 = -imax,imax
          do i = 1, 3
            r(i) = rbas(i,1)*dble(i1)+rbas(i,2)*dble(i2)+rbas(i,3)*dble(i3)
          end do !i
          rps = r+rshift
          rleng = dsqrt(rps(1)*rps(1)+rps(2)*rps(2)+rps(3)*rps(3))
          if ((rleng <= rmax) .and. (rleng > 1.0d-6)) then
            ir = ir+1
            rstr(1:3,ir) = rps(1:3)
            rstr(4,ir) = dsqrt(rps(1)*rps(1)+rps(2)*rps(2)+rps(3)*rps(3))
          endif
        enddo  
      enddo      
    enddo
    nr = ir
    !
    ! sort by increasing length using shell algorithm
    !
    gap = nr/2
    do while (gap >= 1)
      done = .false.
      do while (.not.done)
        done = .true.
        do i = 1, nr-gap
          if (rstr(4,i) > rstr(4,i+gap)) then
            rtmp(1:4) = rstr(1:4,i)
            rstr(1:4,i) = rstr(1:4,i+gap)
            rstr(1:4,i+gap) = rtmp(1:4)
            done = .false.
          endif
        enddo
      enddo
      gap = gap/2
    enddo
end subroutine genrstr
!EOC
