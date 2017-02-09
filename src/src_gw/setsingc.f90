
subroutine setsingc

    use modinput
    use modmain,  only: omega, pi
    use modgw,    only: fdebug, kqset, singc1, singc2

    real(8) :: beta
    real(8) :: f1
    real(8) :: f2
    real(8) :: intf1   ! BZ integral of the auxiliary function 1      
    real(8) :: intf2   ! BZ integral of the auxiliary function 2      
    real(8) :: sumf1   ! Auxiliary functions for 
    real(8) :: sumf2   ! BZ integrals at singular $\Gamma$ point
    integer :: iq
    
    if ((input%gw%debug).and.(myrank==0)) then
      write(fdebug,*)
      write(fdebug,*) 'Use the auxiliary function by Massida et al. (1993)'
      write(fdebug,*)
    end if
      
    beta = (omega/(6.0d0*pi*pi))**(1.0d0/3.0d0)
    intf1 = omega/(4.0d0*pi*pi*beta)
    intf2 = omega/(4.0d0*pi*pi)*sqrt(pi/beta)

    sumf1 = 0.0d0
    sumf2 = 0.0d0
    do iq = 1, kqset%nkpt
      call genauxf(iq,beta,f1,f2)
      sumf1 = sumf1 + f1
      sumf2 = sumf2 + f2
    enddo  
    sumf1 = sumf1/dble(kqset%nkpt)
    sumf2 = sumf2/dble(kqset%nkpt)

    singc1 = intf1-sumf1
    singc2 = intf2-sumf2

    if ((input%gw%debug).and.(myrank==0)) then
      write(fdebug,*) 'Info(setsingc): Integrals of the auxiliary function'
      write(fdebug,1) beta
      write(fdebug,2) intf1, intf2
      write(fdebug,3) sumf1, sumf2
      write(fdebug,4) singc1, singc2
      1 format('Parameter beta: ',f18.12,/,30x,'q^(-1)',12x,'q^(-2)')
      2 format('Analitic integration: ',2f18.12)
      3 format('Numerical integration: ',2f18.12)
      4 format('Correction factor: ',2f18.12)
    end if

contains

    subroutine genauxf(iq,beta,f1,f2)
!
! Given the index \texttt{iq} of $\vec{q}$, this
! subroutine generates the auxiliary functions $F_1(\vec{q})$ and $F_2(\vec{q})$
! according to the formulas:
!
! \begin{subequations}\label{genauxf-01}
! \begin{align}
! F_1(\vec{q})=&\frac{1}{\Omega}\sum\limits_i{\frac{e^{-\alpha|\vec{q}+\vec{G}_i|^2}}{|\vec{q}+\vec{G}_i|}}\\
! F_2(\vec{q})=&\frac{1}{\Omega}\sum\limits_i{\frac{e^{-\alpha|\vec{q}+\vec{G}_i|^2}}{|\vec{q}+\vec{G}_i|^2}}
! \end{align}
! \end{subequations}
! 
! where $\beta=\left(\frac{\Omega}{6\pi^2}\right)^{\frac{1}{3}}$
!
        use modmain
        use modgw,       only : kqset, Gset, Gqset 
        use mod_misc_gw, only : gammapoint

        ! input parameters
        implicit none
        integer, intent(in) :: iq    ! Index of the q-vector       
        real(8), intent(in) :: beta  ! Parameter of the function
      
        ! output parameters
        real(8), intent(out) :: f1   ! Value of the auxiliary function F_1
        real(8), intent(out) :: f2   ! Value of the auxiliary function F_2
      
        ! local variables
        integer :: ipw, ipwin
        real(8) :: gpq(3), gpq2
        real(8) :: expq

        if (gammapoint(kqset%vqc(1:3,iq))) then
          ipwin = 2
        else
          ipwin = 1
        end if

        f1 = 0.0d0
        f2 = 0.0d0

        ! Loop over G-vectors
        do ipw = ipwin, Gqset%ngk(1,iq)

          gpq(1:3) = Gset%vgc(1:3,Gqset%igkig(ipw,1,iq))+kqset%vqc(1:3,iq)
          
          gpq2 = gpq(1)*gpq(1)+gpq(2)*gpq(2)+gpq(3)*gpq(3)
          expq = dexp(-beta*gpq2)
          
          f1 = f1 + expq/dsqrt(gpq2) ! Notice the error in definition of F_1(q) Eq.(FHIgap-B.3)
          f2 = f2 + expq/gpq2
            
        enddo ! ipw
        return
    end subroutine genauxf
    
end subroutine 
