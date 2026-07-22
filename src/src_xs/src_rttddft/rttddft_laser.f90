!> Module with the several laser pulses considered
module rttddft_laser
#include "asserts.fpp"
  use constants, only: pi, twopi
  use modinput, only: kick_type_array, trapCos_type_array, sinSq_type_array
  use physical_constants, only: c
  use precision, only: dp, i32
  use rttddft_VectorField, only: direction, cartesian_direction
  
  implicit none

  private

  character(len=*), parameter :: field_total = 'total'
  character(len=*), parameter :: field_external = 'external'
  real(dp), parameter :: eps_kick_width = 1.e-14_dp

  !> Enum with the type of applied field
  !> There are 2 possibilities that the applied field can assume: "total" or "external"
  enum, bind(C)
    enumerator :: applied_field
    enumerator :: total, external
  end enum 

  !> Abstract type to model a generic laser pulse
  type, abstract :: Laser_Pulse
    private 
    !> Time \( t_0 \) when the laser pulse is applied.
    !> For the kick with nonzero width, it is the center of broadened delta kick
    real(dp) :: t_0
    !> Broadening
    real(dp) :: width
    !> Amplitude 
    real(dp) :: amplitude
    !> Direction: `x`, `y`, `z`
    integer(kind(direction)) :: direction
  contains
    procedure(evaluate_), private, deferred :: evaluate
    procedure(evaluate_), private, deferred :: evaluate_derivative
  end type

  !> Laser pulse whose shape resembles an impulsive electric field. The vector potential,
  !> however, has the shape of a (smoothed) step function
  type, extends(Laser_Pulse) :: Delta_Kick
  contains
    procedure, public :: initialize => initialize_delta_kick
    procedure, private :: evaluate => evaluate_delta
    procedure, private :: evaluate_derivative => evaluate_derivative_delta
  end type

  !> Laser pulse whose shape is given by a cossine function modulated by a 
  !> trapezoidal envelope
  type, extends(Laser_Pulse) :: Cossine_with_Trapezoidal_Envelope
    private
    !> Rise time \( t_r \)
    real(dp) :: t_r
    !> Angular frequency of the cossine function
    real(dp) :: omega
    !> phase of the cossine function
    real(dp) :: phase
  contains
    procedure, public :: initialize => initialize_cos_trapez
    procedure, private :: evaluate => evaluate_cos_trapez
    procedure, private :: evaluate_derivative => evaluate_derivative_cos_trapez
  end type

  !> Laser pulse whose shape is given by a cossine function modulated by a 
  !> sin squared envelope (resembles a Gaussian envelope)
  type, extends(Laser_Pulse) :: Cossine_with_Sinsquared_Envelope
    private
    !> angular frequency \( \omega \) of the cossine function
    real(dp) :: omega
    !> phase of the cossine function
    real(dp) :: phase
  contains
    procedure, public :: initialize => initialize_cos_sinSquared
    procedure, private :: evaluate => evaluate_cos_sinSquared
    procedure, private :: evaluate_derivative => evaluate_derivative_cos_sinSquared
  end type

  type, public :: Set_of_Laser_Pulses
    private
    !> There are 2 possibilities for `given_field`: "total" or "external"
    integer(kind(applied_field)) :: given_field = total
    !> Delta kicks used to build the vector potential
    type(Delta_Kick), allocatable :: delta_kicks(:)
    !> Laser pulses with the shape cossine with trapezoidal envelope
    !> used to build the vector potential
    type(Cossine_with_Trapezoidal_Envelope), allocatable :: cos_trapezoidal(:)
    !> Laser pulses with the shape cossine with sin squared envelope
    !> used to build the vector potential
    type(Cossine_with_Sinsquared_Envelope), allocatable :: cos_sinsquared(:)
  contains
    procedure, public :: set_given_field, initialize_pulses
    procedure, public :: is_external_field_given, is_total_field_given
    procedure, public :: applied_vector_potential, get_dA_dt
    final :: destructor_set_laser_pulses
  end type

  abstract interface
    pure function evaluate_( this, t ) result( r )
      import :: Laser_Pulse, dp
      class(Laser_Pulse), intent(in) :: this
      !> time \( t \)
      real(dp), intent(in)  :: t
      !> the vector potential (or its derivative) at time \( t \)
      real(dp) :: r
    end function
  end interface

contains

  !> Return an enum with the field type
  function field_type( name ) result( field )
    !> String with the name of field type
    character(len=*), intent(in) :: name
    !> resulting enum
    integer(kind(applied_field)) :: field

    CALL_ASSERT( trim( name ) == field_total .or. trim( name ) == field_external ,  'name must be ' // field_total // ' or ' // field_external )
  
    select case( trim( name ) )
      case( field_total )
        field = total
      case( field_external )
        field = external
    end select
  end function

  !> Initialize a delta kick. See [[evaluate_delta]] for a complete documentation
  subroutine initialize_delta_kick( this, t_0, width, amplitude, dir )
    class(Delta_Kick), intent(inout) :: this
    !> Time \(t_0\) with the center of the delta pulse
    real(dp), intent(in) :: t_0
    !> Broadening width of the delta pulse
    real(dp), intent(in) :: width
    !> Amplitude \(E_0\)
    real(dp), intent(in) :: amplitude
    !> Character containing the cartesian direction: must be `x`, `y`, or `z`
    character, intent(in) :: dir

    this%t_0 = t_0
    this%width = width
    this%amplitude = -c*amplitude
    this%direction = cartesian_direction( dir )
  end subroutine

  !> Delta Kick function
  !> IMPORTANT: The delta kick is meant for the electric field
  !> \[
  !>    \mathbf{E} = -\frac{1}{c}\frac{d\mathbf{A}}{dt}
  !> \]
  !> The Vector potential is, therefore, a step function.
  !> If we consider strictly a delta kick, we have
  !> \[
  !>    \mathbf{E} = \mathbf{E}_0 \delta(t - t_0).
  !> \]
  !> It is possible to broaden it as
  !> \[
  !>    \mathbf{E}(t) = \mathbf{E}_0 \frac{15}{16 w}
  !>     \left( \frac{t-t_0}{w} +1 \right)^2 \left( \frac{t-t_0}{w} -1 \right)^2
  !> \]
  !> for \( t \) between \( t_0 - w \) and \( t_0 + w \), and zero otherwise.
  !> This function is smooth, and has a maximum value of \( \frac{15 \mathbf{E}_0 }{16 w} \)  on \( t_0 \).
  !> The resulting vector potential is
  !> \[
  !>    \mathbf{A}(t) = -c\mathbf{E}_0 \frac{1}{16}
  !>     \left( \frac{t+t_0}{w} +1 \right)^3
  !>     \left[ 3\left( \frac{t-t_0}{w} \right)^2
  !>      - 9\left( \frac{t-t_0}{w} \right) + 8 \right]
  !> \]
  !> for \( t \) between \( t_0 - w \) and \( t_0 + w \),
  !> zero for \( t < t_0 - w \) and \( -c\mathbf{E}_0 \) for \( t > t_0 + w \).
  !> If \(w = 0 \), the vector potential is zero for \( t < t_0 \) and \( -c\mathbf{E}_0 \) for \( t > t_0 \)
  pure function evaluate_delta( this, t ) result(a)
    class(Delta_Kick), intent(in) :: this
    !> time \( t \)
    real(dp), intent(in)  :: t
    !> the calculated vector potential at time \( t \)
    real(dp) :: a

    real(dp) :: tsh

    if ( abs( this%width ) > eps_kick_width ) then
      tsh = ( t - this%t_0 )/this%width
      if( tsh >= 1._dp ) then
        a = this%amplitude
      elseif( ( tsh >= -1._dp) .and. ( tsh <= 1._dp ) ) then
        a = (1._dp/16._dp) * this%amplitude * &
          & ( 3._dp*tsh**2 - 9._dp*tsh + 8._dp)*( tsh + 1._dp )**3
      else
        a = 0._dp
      end if
    else
      if ( t >= this%t_0 ) then
        a = this%amplitude
      else
        a = 0._dp
      end if
    end if
  end function

  pure function evaluate_derivative_delta( this, t ) result( a )
    class(Delta_Kick), intent(in) :: this
    !> time \( t \)
    real(dp), intent(in) :: t
    !> the derivative of the vector potential at time \( t \)
    real(dp) :: a

    real(dp) :: t_aux

    a = 0._dp
    if ( abs( this%width ) > eps_kick_width ) then
      t_aux = ( t - this%t_0 ) / this%width
      if( ( t_aux >= -1._dp) .and. ( t_aux <= 1._dp ) ) then
        a = (15._dp/16._dp) * this%amplitude * (t_aux + 1._dp)**2*(t_aux - 1._dp)**2 / this%width
      end if
    else
      if ( t == this%t_0 ) a = huge( 1._dp )
    end if
  end function

  !> Initialize cossine function with trapezoidal envelope. See [[evaluate_cos_trapez]] for a complete documentation
  subroutine initialize_cos_trapez( this, t_0, width, amplitude, dir, t_rise, omega, phase )
    class(Cossine_with_Trapezoidal_Envelope), intent(inout) :: this
    !> Time \(t_0\) when the pulse is applied 
    real(dp), intent(in) :: t_0
    !> Width of the laser pulse
    real(dp), intent(in) :: width
    !> Amplitude \(A_0\)
    real(dp), intent(in) :: amplitude
    !> Character containing the cartesian direction: must be `x`, `y`, or `z`
    character, intent(in) :: dir
    !> Rise time \(t_r\)
    real(dp), intent(in) :: t_rise
    !> Angular frequency \(\omega\)
    real(dp), intent(in) :: omega
    !> Phase \(\phi\)
    real(dp), intent(in) :: phase
    
    this%t_0 = t_0
    this%width = width
    this%amplitude = amplitude
    this%direction = cartesian_direction( dir )
    this%t_r = t_rise
    this%omega = omega
    this%phase = phase
  end subroutine

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
  pure function evaluate_cos_trapez( this, t ) result( a )
    class(Cossine_with_Trapezoidal_Envelope), intent(in) :: this
    !> Time \( t \)
    real(dp), intent(in)  :: t
    !> the calculated vector potential at time t
    real(dp) :: a

    real(dp)              :: t_aux

    associate( t_0=>this%t_0, width=>this%width, tr=>this%t_r, &
               amplitude=>this%amplitude, omega=>this%omega, phase=>this%phase )
      t_aux = (t - t_0)
      if ( (t_aux >= 0._dp) .and. ( t_aux <= 2._dp*tr + width ) ) then
        if( t_aux < tr ) then
          a = t_aux/tr
        else if ( t_aux <= tr + width ) then
          a = 1._dp
        else
          a = ( 2._dp * tr + width - t_aux )/tr
        end if
        a = a*amplitude*cos( omega*t + phase )
      else
        a = 0._dp
      end if
    end associate
  end function

  pure function evaluate_derivative_cos_trapez( this, t ) result( a )
    class(Cossine_with_Trapezoidal_Envelope), intent(in) :: this
    !> time \( t \)
    real(dp), intent(in)  :: t
    !> the derivative of the vector potential at time \( t \)
    real(dp) :: a

    real(dp) :: t_aux, env, env_dot

    a = 0._dp
    associate( t_0=>this%t_0, width=>this%width, tr=>this%t_r, &
               amplitude=>this%amplitude, omega=>this%omega, phase=>this%phase )
      t_aux = t - t_0
      if ( (t_aux >= 0._dp) .and. (t_aux <= 2._dp * tr + width) ) then
        if( t_aux < tr ) then
          env = t_aux / tr
          env_dot = 1._dp / tr
        else if ( t_aux >= tr + width ) then
          env = (width + 2._dp * tr - t_aux) / tr
          env_dot = -1._dp / tr
        else
          env = 1._dp
          env_dot = 0._dp
        end if
        a = amplitude * (env_dot * cos( omega * t + phase ) - env * omega * sin( omega * t + phase ))
      end if
    end associate
  end function

  !> Initialize cossine function with sin squared envelope. See [[evaluate_cos_sinSquared]] for a complete documentation
  subroutine initialize_cos_sinSquared( this, t_0, width, amplitude, dir, omega, phase )
    class(Cossine_with_Sinsquared_Envelope), intent(inout) :: this
    !> Time \(t_0\) when the pulse is applied
    real(dp), intent(in) :: t_0
    !> Width of the laser pulse
    real(dp), intent(in) :: width
    !> Amplitude \(A_0\)
    real(dp), intent(in) :: amplitude
    !> Character containing the cartesian direction: must be `x`, `y`, or `z`
    character, intent(in) :: dir
    !> Angular frequency \(\omega\)
    real(dp), intent(in) :: omega
    !> Phase \(\phi\)
    real(dp), intent(in) :: phase

    this%t_0 = t_0
    this%width = width
    this%amplitude = amplitude
    this%direction = cartesian_direction( dir )
    this%omega = omega
    this%phase = phase
  end subroutine

  !> Cossine function modulated by a sine squared  
  !> \[
  !>    \mathbf{A}(t) = \mathbf{A}_0  f(t) \cos ( \omega t + \phi )
  !> \]
  !> where \( f(t) \) is the following sine squared function:  
  !> <ul>
  !> <li> \( f(t) = 0 \), if \( t \le t_0 \) or \( t \ge t_0 + w\) </li>
  !> <li> \( f(t) = \sin^2( \pi(t-t_0)/w) \), if \( t_0 \le t \le t_0 + w \) </li>
  !> </ul>
  pure function evaluate_cos_sinSquared( this, t ) result( a )
    class(Cossine_with_Sinsquared_Envelope), intent(in) :: this
    !> time \( t \)
    real(dp), intent(in)  :: t
    !> vector potential at time \( t \)
    real(dp) :: a

    real(dp)              :: t_aux

    associate( t_0=>this%t_0, width=>this%width, amplitude=>this%amplitude, &
               omega=>this%omega, phase=>this%phase )
      t_aux = (t-t_0)
      if ( ( t_aux >= 0._dp ) .and. (t_aux <= width ) ) then
        a = amplitude*( sin( pi*t_aux/width )**2 )*cos( omega*t + phase )
      else
        a = 0._dp
      end if
    end associate
  end function

  pure function evaluate_derivative_cos_sinSquared( this, t ) result( a )
    class(Cossine_with_Sinsquared_Envelope), intent(in) :: this
    !> time \( t \)
    real(dp), intent(in)  :: t
    !> the calculated vector potential at time \( t \)
    real(dp) :: a

    real(dp)              :: t_aux

    associate( t_0=>this%t_0, width=>this%width, amplitude=>this%amplitude, &
               omega=>this%omega, phase=>this%phase )
      t_aux = (t-t_0)
      if ( ( t_aux >= 0._dp ) .and. (t_aux <= width ) ) then
        a = (amplitude)*( -omega*(sin( pi*t_aux/width )**2)*sin( omega*t + phase ) + &
          (pi/width)*sin( twopi*t_aux/width )*cos( omega*t + phase ) )
      else
        a = 0._dp
      end if
    end associate
  end function

  !> Obtain the total vector potential at time \(t\) due to an array of laser pulses
  pure function get_vector_potential_from_lasers( lasers, t ) result(a)
    !> Array of laser pulses
    class(Laser_Pulse), contiguous, intent(in) :: lasers(:)
    !> Time \(t\)
    real(dp), intent(in) :: t
    !> Cartesian components `x`, `y`, and `z` of the vector potential
    real(dp) :: a(3)

    integer(i32) :: i
    
    a = 0._dp
    do i = 1, size( lasers )
      a( lasers(i)%direction ) = a( lasers(i)%direction ) + lasers(i)%evaluate( t )
    end do
  end function

  !> Obtain the total vector potential at time \(t\) due to a set of laser pulses
  pure function applied_vector_potential( this, t ) result( a )
    !> Set of laser pulses
    class(Set_of_Laser_Pulses), intent(in) :: this
    !> Time \(t\)
    real(dp), intent(in) :: t
    !> Cartesian components `x`, `y`, and `z` of the vector potential
    real(dp) :: a(3)

    a = get_vector_potential_from_lasers( this%delta_kicks, t )
    a = a + get_vector_potential_from_lasers( this%cos_trapezoidal, t )
    a = a + get_vector_potential_from_lasers( this%cos_sinsquared, t )
  end function

  !> Obtain the the derivative of the vector potential \(dA/dt\) due to an array of laser pulses
  pure function get_dA_dt_from_lasers( lasers, t ) result( da_dt )
    !> Array of laser pulses
    class(Laser_Pulse), contiguous, intent(in) :: lasers(:)
    !> time \(t\)
    real(dp), intent(in) :: t
    !> Cartesian components `x`, `y`, and `z` of the \(dA/dt\)
    real(dp) :: da_dt(3)

    integer(i32) :: i, m
    
    m = size( lasers )
    da_dt = 0._dp
    do i = 1, m
      da_dt( lasers(i)%direction ) = da_dt( lasers(i)%direction ) + lasers(i)%evaluate_derivative( t )
    end do
  end function

  !> Obtain the vector potential \( dA(t)/dt \) due to all applied lasers pulses
  pure function get_dA_dt( this, t ) result( da_dt )
    !> Set of laser pulses
    class(Set_of_Laser_Pulses), intent(in) :: this
    !> time \(t\)
    real(dp), intent(in) :: t
    !> Cartesian components `x`, `y`, and `z` of the \(dA/dt\)
    real(dp) :: da_dt(3)

    da_dt = get_dA_dt_from_lasers( this%delta_kicks, t )
    da_dt = da_dt + get_dA_dt_from_lasers( this%cos_trapezoidal, t )
    da_dt = da_dt + get_dA_dt_from_lasers( this%cos_sinsquared, t )
  end function

  pure logical function is_external_field_given( this )
    class(Set_of_Laser_Pulses), intent(in) :: this
    is_external_field_given = (this%given_field == external)
  end function

  pure logical function is_total_field_given( this )
    class(Set_of_Laser_Pulses), intent(in) :: this
    is_total_field_given = (this%given_field == total)
  end function

  !> Set the component `given_field` of a `Set_of_Laser_Pulses` object as the 
  !> enum corresponding to the string `given_field` as dummy argument
  subroutine set_given_field( this, given_field )
    class(Set_of_Laser_Pulses), intent(inout) :: this
    !> String spelling out the name of the given field. It must be equal to `field_external` or `field_total`
    character(len=*), intent(in) :: given_field

    this%given_field = field_type( given_field )

  end subroutine

  !> Initialize the object of class `Set_of_Laser_Pulses`
  subroutine initialize_pulses( this, kick_array, trapCos_array, sinSq_array )
    class(Set_of_Laser_Pulses), intent(inout) :: this
    !> Array with kick pulses
    type(kick_type_array), intent(in) :: kick_array(:)
    !> Array with trapCos pulses
    type(trapCos_type_array), intent(in) :: trapCos_array(:)
    !> Array with sinSq pulses
    type(sinSq_type_array), intent(in) :: sinSq_array(:)

    integer(i32) :: i, m

    m = size( kick_array )
    allocate( this%delta_kicks(m) )
    do i = 1, m
      associate( kick => kick_array(i)%kick )
        call this%delta_kicks(i)%initialize( kick%t0, kick%width, kick%amplitude, kick%direction )
      end associate
    end do

    m = size( trapCos_array )
    allocate( this%cos_trapezoidal(m) )
    do i = 1, m
      associate( trap => trapCos_array(i)%trapCos )
        call this%cos_trapezoidal(i)%initialize( trap%t0, trap%width, trap%amplitude, &
          trap%direction, trap%riseTime, trap%omega, trap%phase )
      end associate
    end do

    m = size( sinSq_array )
    allocate( this%cos_sinsquared(m) )
    do i = 1, m
      associate( sinsq => sinSq_array(i)%sinSq )
        call this%cos_sinsquared(i)%initialize( sinsq%t0, sinsq%pulseLength, sinsq%amplitude, &
          sinsq%direction, sinsq%omega, sinsq%phase )
      end associate
    end do

  end subroutine

  elemental impure subroutine destructor_set_laser_pulses( this )
    type(Set_of_Laser_Pulses), intent(inout) :: this
    if( allocated(this%delta_kicks) ) deallocate( this%delta_kicks )
    if( allocated(this%cos_trapezoidal) ) deallocate( this%cos_trapezoidal )
    if( allocated(this%cos_sinsquared) ) deallocate( this%cos_sinsquared )
  end subroutine


end module