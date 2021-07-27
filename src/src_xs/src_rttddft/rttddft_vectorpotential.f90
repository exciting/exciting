! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
! Copyright (C) Exciting Code, SOL group. 2020

! HISTORY
! Created July 2019 (Ronaldo)
! Improved documentation: July 2021 (Ronaldo)
! Reference: https://doi.org/10.1088/2516-1075/ac0c26

!> Module to manage the vector potential \( \mathbf{A} \)
module rttddft_VectorPotential

  use precision, only: dp

  implicit none

  private :: Delta_Kick, Cossine_with_Trapezoidal_Envelope, &
    & Cossine_with_Sinsquared_Envelope
  public  :: Calculate_Vector_Potential, Solve_ODE_Vector_Potential

contains
  !> Delta Kick function
  !> IMPORTANT: The delta kick is meant for the electric field
  !> \[
  !>    \mathbf{E} = -\frac{1}{c}\frac{d\mathbf{A}}{dt}
  !> \]
  !> The Vector potential is therefore a step function
  !> If we consider strictly a delta kick, we have
  !> \[
  !>    \mathbf{E} = \mathbf{E}_0 \delta(t-t_0)
  !> \]
  !> this means width \( w = 0 \).
  !> However, it is possible to broaden it as
  !> \[
  !>    \mathbf{E}(t) = \mathbf{E}_0 \frac{15}{16}
  !>     \left( \frac{t-t_0}{w} +1 \right)^2 \left( \frac{t-t_0}{w} -1 \right)^2
  !> \]
  !> for \( t \) between \( t_0 - w \) and \( t_0 + w \), and zero otherwise.
  !> This function is smooth, and has a maximum on \( t_0 \).
  !> The resulting vector potential is
  !> \[
  !>    \mathbf{A}(t) = -c\mathbf{E}_0 \frac{1}{16}
  !>     \left( \frac{t+t_0}{w} +1 \right)^3
  !>     \left[ 3\left( \frac{t-t_0}{w} \right)^2
  !>      - 9\left( \frac{t-t_0}{w} \right) + 8 \right]
  !> \]
  !> for \( t \) between \( t_0 - w \) and \( t_0 + w \),
  !> zero for \( t < t_0 - w \) and \( -c\mathbf{E}_0 \) for \( t > t_0 + w \).
  subroutine Delta_Kick( t, t0, width, amplitude, a )
    implicit none

    !> time \( t \)
    real(dp), intent(in)  :: t
    !> time \( t_0 \) when the kick is applied, or in the case of
    !> nonzero width, it is the center of broadened delta kick
    real(dp), intent(in)  :: t0
    !> broadening of the delta kick
    real(dp), intent(in)  :: width
    !> amplitude \( -cE_0 \)
    real(dp), intent(in)  :: amplitude
    !> the calculated vector potential at time t
    real(dp), intent(out) :: a

    real(dp)              :: tsh

    if ( width /= 0._dp ) then
      tsh = ( t - t0 )/width
      if( tsh >= 1._dp ) then
        a = amplitude
      elseif( ( tsh >= -1._dp) .and. ( tsh .le. 1._dp ) ) then
        a = (1._dp/16._dp) * amplitude * &
          & ( 3._dp*tsh**2 - 9._dp*tsh + 8._dp)*( tsh + 1._dp )**3
      else
        a = 0._dp
      end if
    else
      if ( t >= t0 ) then
        a = amplitude
      else
        a = 0._dp
      end if
    end if
  end subroutine Delta_Kick

  !> Cossine function modulated by a trapezoid
  !> \[
  !>    \mathbf{A}(t) = \mathbf{A}_0  f(t) \cos ( \omega t + \phi )
  !> \]
  !> where \( f(t) \) is the trapezoidal function:
  !> <ul>
  !> <li> \( f(t) = 0 \), if \( t \le t_0 \) or \( t \ge t_0 + w + 2 t_r\) </li>
  !> <li> \( f(t) = 1 \), if \( t_0 + t_r \le t \le t_0 + t_r + w \) </li>
  !> <li> \( f(t) = (t-t_0)/t_r \), if \( t_0 < t < t_0 + t_r \) </li>
  !> <li> \( f(t) = (t_0 + w + 2 t_r - t )/t_r \), if \( t_0 + t_r + w < t < t_0 + w + 2 t_r \) </li>
  !> </ul>
  subroutine Cossine_with_Trapezoidal_Envelope( t, t0, tr, width, omega, phase, &
    & amplitude, a )
    implicit none

    !> time \( t \)
    real(dp), intent(in)  :: t
    !> time \( t_0 \) as defined above
    real(dp), intent(in)  :: t0
    !> rise time \( t_r \) as defined above
    real(dp), intent(in)  :: tr
    !> witdh of the trapezoid
    !> (width of the interval where \( f(t) = 1 \) ).
    real(dp), intent(in)  :: width
    !> angular frequency of the cossine function
    real(dp), intent(in)  :: omega
    !> phase of the cossine function
    real(dp), intent(in)  :: phase
    !> amplitude \( A_0 \)
    real(dp), intent(in)  :: amplitude
    !> the calculated vector potential at time t
    real(dp), intent(out) :: a

    real(dp)              :: tsh

    tsh = (t - t0)
    if ( (tsh >= 0._dp) .and. ( tsh <= 2._dp*tr + width ) ) then
      if( tsh < tr ) then
        a = tsh/tr
      else if ( tsh <= tr + width ) then
        a = 1._dp
      else
        a = ( 2*tr + width - tsh )/tr
      end if
      a = a*amplitude*cos( omega*t + phase )
    else
      a = 0._dp
    end if
  end subroutine Cossine_with_Trapezoidal_Envelope

  !> Cossine function modulated by a sine squared  
  !> \[
  !>    \mathbf{A}(t) = \mathbf{A}_0  f(t) \cos ( \omega t + \phi )
  !> \]
  !> where \( f(t) \) is the following sine squared function:  
  !> <ul>
  !> <li> \( f(t) = 0 \), if \( t \le t_0 \) or \( t \ge t_0 + w\) </li>
  !> <li> \( f(t) = \sin^2( \pi(t-t_0)/w) \), if \( t_0 \le t \le t_0 + w \) </li>
  !> </ul>
  subroutine Cossine_with_Sinsquared_Envelope( t, t0, width, omega, phase, &
    & amplitude, a )
    use constants, only: pi
    implicit none

    !> time \( t \)
    real(dp), intent(in)  :: t
    !> time \( t_0 \) as defined above
    real(dp), intent(in)  :: t0
    !> witdh of the sine squared
    real(dp), intent(in)  :: width
    !> angular frequency \( \omega \) of the cossine function
    real(dp), intent(in)  :: omega
    !> phase of the cossine function
    real(dp), intent(in)  :: phase
    !> amplitude \( A_0 \)
    real(dp), intent(in)  :: amplitude
    !> the calculated vector potential at time \( t \)
    real(dp), intent(out) :: a

    real(dp)              :: tsh

    tsh = (t-t0)
    if ( ( tsh >= 0._dp ) .and. (tsh <= width ) ) then
      a = amplitude*( sin( pi*tsh/width )**2 )*cos( omega*t + phase )
    else
      a = 0._dp
    end if
  end subroutine Cossine_with_Sinsquared_Envelope

  !> Obtain the vector potential \( A(t) \) due to a laser pulse
  subroutine Calculate_Vector_Potential( t, aapplied )
    use rttddft_GlobalVariables, only: nkicks, t0kick, wkick, amplkick, dirkick, &
      & ntrapcos, dirtrapcos, ampltrapcos, omegatrapcos, phasetrapcos, &
      & t0trapcos, trtrapcos, wtrapcos, &
      & nsinsq, dirsinsq, amplsinsq, omegasinsq, phasesinsq, &
      & t0sinsq, tpulsesinsq

    implicit none

    !> time \( t \)
    real(dp),intent(in)         :: t
    !> the vector potential (components `x`, `y`, and `z`)
    real(dp),intent(out)        :: aapplied(3)

    integer                     :: i
    real(dp)                    :: amod

    aapplied(:) = 0._dp

    do i = 1, nkicks
      call Delta_Kick( t, t0kick(i), wkick(i), amplkick(i), amod )
      select case( dirkick(i) )
        case('x')
          aapplied(1) = aapplied(1) + amod
        case('y')
          aapplied(2) = aapplied(2) + amod
        case('z')
          aapplied(3) = aapplied(3) + amod
      end select
    end do

    ! Cossine pulses modulated by a trapezoidal function
    do i = 1, ntrapcos
      call Cossine_with_Trapezoidal_Envelope( t, t0trapcos(i), trtrapcos(i), &
        & wtrapcos(i), omegatrapcos(i), phasetrapcos(i), ampltrapcos(i), amod )
      select case( dirtrapcos(i) )
        case('x')
          aapplied(1) = aapplied(1) + amod
        case('y')
          aapplied(2) = aapplied(2) + amod
        case('z')
          aapplied(3) = aapplied(3) + amod
      end select
    end do

    ! Cossine pulses modulated by a sin squared function
    do i = 1, nsinsq
      call Cossine_with_Sinsquared_Envelope( t, t0sinsq(i), tpulsesinsq(i), &
        & omegasinsq(i), phasesinsq(i), amplsinsq(i), amod )
      select case(dirsinsq(i))
        case('x')
          aapplied(1) = aapplied(1) + amod
        case('y')
          aapplied(2) = aapplied(2) + amod
        case('z')
          aapplied(3) = aapplied(3) + amod
      end select
    end do
  end subroutine Calculate_Vector_Potential


  !> Update the induced vector potential \( \mathbf{A}_{ind} \)
  !> using the current density \( \mathbf{J}(t) \)
  !> We need to solve the differential equation
  !> \[
  !>  \frac{d^2\mathbf{A}_{ind}}{dt^2} = 4 \pi c \mathbf{J}(t).
  !>  \]
  subroutine Solve_ODE_Vector_Potential( TotalFieldIsGiven )
    use rttddft_GlobalVariables, only: tstep, time, jpara, jparanext, &
      & aind, pvec, jind, aext, atot
    use errors_warnings, only: terminate_if_false
    use constants, only: fourpi, pi
    use physical_constants, only: c
    use modmpi, only: mpiglobal
    use mod_lattice, only: omega
    use mod_charge_and_moment, only: chgval
    use modinput, only: input

    implicit none
    !> Tells if the total field (`TotalFieldIsGiven` = true) or if the 
    !> external field (`TotalFieldIsGiven` = false) is given by the laser.  
    !> Reminder: \( \mathbf{A} = \mathbf{A}_{ext} + \mathbf{A}_{ind} \)
    logical, intent(in)   :: TotalFieldIsGiven

    real(dp)              :: beta, fac, den
    real(dp)              :: k1(3,2), k2(3,2), k3(3,2), k4(3,2)
    real(dp)              :: jparamid(3),jindmid(3),jindnext(3)
    real(dp)              :: aauxmid(3),aauxnext(3),smid(3)
    real(dp)              :: asave(3)

    beta = chgval/c/omega
    ! Method of integrating the differential equation
    select case(input%xs%realTimeTDDFT%vectorPotentialSolver)
      case('euler') ! Euler
        aind = aind + fourpi*c*tstep*pvec
        pvec = pvec + tstep*jind
      case('improvedeuler') ! Improved Euler method
        aind = aind + fourpi*c*tstep*(pvec + (0.5_dp)*(tstep)*jind)
        call Calculate_Vector_Potential( time, aauxnext )
        if ( TotalFieldIsGiven ) then
          jindnext = jparanext - beta*(aauxnext)
        else
          jindnext = jparanext - beta*( aind + aauxnext )
        end if
        pvec = pvec + (0.5_dp)*tstep*( jind + jindnext )
      case('midpoint')
        call Calculate_Vector_Potential( time, aauxnext )
        if ( TotalFieldIsGiven ) then
          jindnext = jparanext - beta*( aauxnext )
          jindmid = 0.5_dp*( jind + jindnext )
          aind = aind + fourpi*c*tstep*( pvec + 0.5_dp*tstep*jindmid )
          pvec = pvec + tstep*jindmid
        else
          asave = aind
          jparamid = (0.5_dp)*( jparanext + jpara )
          smid = jparamid(:) - 0.5_dp*beta*( aext + aauxnext )
          fac = pi*beta*c*(tstep**2)
          den = 1_dp + fac
          fac = (1_dp - fac)/den
          aind = (fourpi*c*tstep/den)*( pvec + 0.5d0*tstep*smid ) + &
            & fac*aind(:)
          pvec = (tstep/den)*(smid - beta*asave ) + fac*pvec
        end if
      case('rk4') ! Runge-Kutta 4th order
        ! Before we begin with rk4, we need to extrapolate jpara and aext
        jparamid = (0.5_dp)*( jparanext + jpara )
        call Calculate_Vector_Potential( time-(0.5_dp)*tstep, aauxmid )
        call Calculate_Vector_Potential( time, aauxnext )
        ! Now, we apply Runge Kutta of 4th order
        if ( TotalFieldIsGiven ) then
          k1(:,1) = jpara(:) - beta*atot(:)
          k1(:,2) = fourpi*c*pvec(:)
          k2(:,1) = jparamid(:) - beta*aauxmid(:)
          k2(:,2) = fourpi*c*(pvec(:) + (tstep/2._dp)*k1(:,1))
          k3(:,1) = jparamid(:) - beta*aauxmid(:)
          k3(:,2) = fourpi*c*(pvec(:) + (tstep/2._dp)*k2(:,1))
          k4(:,1) = jparanext(:) - beta*aauxnext(:)
          k4(:,2) = fourpi*c*(pvec(:) + (tstep)*k3(:,1))
        else
          k1(:,1) = jpara(:) - beta*(aext(:) + aind(:))
          k1(:,2) = fourpi*c*pvec(:)
          k2(:,1) = jparamid(:) - beta*(aauxmid(:) + aind(:) + (tstep/2._dp)*k1(:,2) )
          k2(:,2) = fourpi*c*(pvec(:) + (tstep/2._dp)*k1(:,1))
          k3(:,1) = jparamid(:) - beta*(aauxmid(:) + aind(:) + (tstep/2._dp)*k2(:,2) )
          k3(:,2) = fourpi*c*(pvec(:) + (tstep/2._dp)*k2(:,1))
          k4(:,1) = jparanext(:) - beta*(aauxnext(:) + aind(:) + (tstep)*k3(:,2) )
          k4(:,2) = fourpi*c*(pvec(:) + (tstep)*k3(:,1))
        end if
        pvec(:) = pvec(:) + (tstep/6._dp)*( &
          & k1(:,1) + 2._dp*k2(:,1) + 2._dp*k3(:,1) + k4(:,1) )
        aind(:) = aind(:) + (tstep/6._dp)*( &
          & k1(:,2) + 2._dp*k2(:,2) + 2._dp*k3(:,2) + k4(:,2) )
      case default
        ! Method not recognized
        ! We need to stop the code
        call terminate_if_false( mpiglobal, .false., &
          & 'Error(Solve_ODE_Vector_Potential): method given in &
          & input%xs%rt_tddft%updateAind is not recognized.' )

    end select
  end subroutine Solve_ODE_Vector_Potential

end module rttddft_VectorPotential
