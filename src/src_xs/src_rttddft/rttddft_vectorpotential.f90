! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
! Copyright (C) Exciting Code, SOL group. 2020

! HISTORY
! Created July 2019 (Ronaldo)
! Improved documentation: July 2021 (Ronaldo)
! Refactored May 2024 (Ronaldo)
! Reference: https://doi.org/10.1088/2516-1075/ac0c26

!> Module to manage the vector potential \( \mathbf{A} \)
module rttddft_VectorPotential
  use constants, only: fourpi, pi
  use mod_lattice, only: omega
  use mod_charge_and_moment, only: chgval
  use modinput, only: laser_type
  use modmpi, only: terminate
  use physical_constants, only: c
  use precision, only: dp, i32
  use rttddft_laser, only: Set_of_Laser_Pulses
  use rttddft_VectorField, only: Uniform_Vector_Field

  implicit none

  private

  public  :: update_a_ind_and_p_vec

  character(len=*), parameter :: solver_euler = 'euler'
  character(len=*), parameter :: solver_improved_euler = 'improvedeuler'
  character(len=*), parameter :: solver_midpoint = 'midpoint'
  character(len=*), parameter :: solver_rk4 = 'rk4'

  !> Enum with the solver type for the vector potential
  enum, bind(C)
    enumerator :: solver_types
    enumerator :: euler, improved_euler, midpoint, rk4
  end enum

  type, public, extends(Uniform_Vector_Field) :: Vector_Potential_Field
  end type

  type, public, extends(Set_of_Laser_Pulses) :: Vector_Potential
    !> \(A_{ind}\): vector potential (induced part) with its `x`, `y`, and `z` components
    type(Vector_Potential_Field) :: a_ind
    !> \(A = A_{ind} + A_{ext}\): vector potential (total = induced + external) with its `x`, `y`, and `z` components
    type(Vector_Potential_Field) :: a_tot
    !> Type of solver used for the vector potential
    integer(kind(solver_types)), private  :: vector_potential_solver
  contains
    procedure :: initialize => initialize_from_input
    procedure :: set_a_tot_a_ind
    procedure :: evaluate_a_tot
    procedure :: is_solver_euler => a_vec_is_solver_euler
  end type

  interface Vector_Potential_Field
    module procedure :: Vector_Potential_Field_Constructor
  end interface

contains
  pure function Vector_Potential_Field_Constructor( vector ) result(a)
    real(dp), intent(in) :: vector(3)
    type(Vector_Potential_Field) :: a

    a%components = vector
  end function  

  function solver_type( solver_name ) result(solver)
    character(len=*), intent(in) :: solver_name
    integer(kind(solver_types)) :: solver 

    select case( trim( solver_name ) )
      case( solver_euler )
        solver = euler
      case( solver_improved_euler )
        solver = improved_euler
      case( solver_midpoint )
        solver = midpoint
      case( solver_rk4 )
        solver = rk4
      case default
        call terminate('unknown solver_type')
    end select
  end function

  !> Initialize the interface to the input variables that define the laser pulses
  subroutine initialize_from_input( this, laser, vectorPotentialSolver )
    class(Vector_Potential), intent(inout) :: this
    !> Type with the variables given in the input file
    type(laser_type), intent(in) :: laser
    !> Method used to update the vector potential
    character(len=*), intent(in) :: vectorPotentialSolver

    this%vector_potential_solver = solver_type( vectorPotentialSolver )
    call this%set_given_field( laser%fieldType )
    call this%initialize_pulses( laser%kickarray, laser%trapCosarray, laser%sinSqarray )
  end subroutine

  pure subroutine evaluate_a_tot( this, t )
    class(Vector_Potential), intent(inout) :: this
    !> time \(t\)
    real(dp), intent(in) :: t

    this%a_tot%components = this%applied_vector_potential( t ) 
    if( this%is_external_field_given() ) this%a_tot%components = this%a_tot%components + this%a_ind%components 
  end subroutine

  pure subroutine set_a_tot_a_ind( this, a_tot, a_ind )
    class(Vector_Potential), intent(inout) :: this
    class(Vector_Potential_Field), intent(in) :: a_tot, a_ind
    
    this%a_tot%components = a_tot%components
    this%a_ind%components = a_ind%components
  end subroutine

  pure logical function a_vec_is_solver_euler( this )
    class(Vector_Potential), intent(in) :: this

    a_vec_is_solver_euler = ( this%vector_potential_solver == euler )
  end function

  !> Update the induced vector potential \( \mathbf{A}_{ind} \)
  !> using the current density \( \mathbf{J}(t) \)
  !> We need to solve the differential equation
  !> \[
  !>  \frac{d^2\mathbf{A}_{ind}}{dt^2} = 4 \pi c \mathbf{J}(t).
  !>  \]
  subroutine update_a_ind_and_p_vec( this, time, dt, jind, jpara, jparanext, pvec )
    class(Vector_Potential), intent(inout) :: this
    !> Time \( t \)
    real(dp), intent(in)      :: time
    !> Time step \( \Delta t \)
    real(dp), intent(in)      :: dt
    !> Total (induced) current density
    real(dp), intent(in) :: jind(3)
    !> Paramagnetic component of the current density
    real(dp), intent(in) :: jpara(3)
    !> Auxiliary variable, used to evolve `jpara`
    real(dp), intent(in) :: jparanext(3)
    !> Polarization field
    real(dp), intent(inout) :: pvec(3)

    real(dp)              :: beta, fac, den
    real(dp)              :: k1(3,2), k2(3,2), k3(3,2), k4(3,2)
    real(dp)              :: jparamid(3),jindmid(3),jindnext(3)
    real(dp)              :: aauxmid(3),aauxnext(3),smid(3)
    real(dp)              :: asave(3)

    beta = chgval/c/omega
    ! Method of integrating the differential equation
    select case( this%vector_potential_solver )
      case( euler ) ! Euler
        call this%a_ind%add_vector( fourpi*c*dt*pvec )
        pvec = pvec + dt*jind
      case( improved_euler ) ! Improved Euler method
        call this%a_ind%add_vector( fourpi*c*dt*(pvec + (0.5_dp)*(dt)*jind) )
        aauxnext = this%applied_vector_potential( time )
        if ( this%is_total_field_given() ) then
          jindnext = jparanext - beta*(aauxnext)
        else
          jindnext = jparanext - beta*( this%a_ind%components + aauxnext )
        end if
        pvec = pvec + (0.5_dp)*dt*( jind + jindnext )
      case( midpoint )
        aauxnext = this%applied_vector_potential( time )
        if ( this%is_total_field_given() ) then
          jindnext = jparanext - beta*( aauxnext )
          jindmid = 0.5_dp*( jind + jindnext )
          call this%a_ind%add_vector( fourpi*c*dt*( pvec + 0.5_dp*dt*jindmid ) )
          pvec = pvec + dt*jindmid
        else
          asave = this%a_ind%components
          jparamid = (0.5_dp)*( jparanext + jpara )
          smid = jparamid(:) - 0.5_dp*beta*( this%a_tot%components - this%a_ind%components + aauxnext )
          fac = pi*beta*c*(dt**2)
          den = 1_dp + fac
          fac = (1_dp - fac)/den
          this%a_ind%components = (fourpi*c*dt/den)*( pvec + 0.5_dp*dt*smid ) + fac*this%a_ind%components
          pvec = (dt/den)*(smid - beta*asave ) + fac*pvec
        end if
      case( rk4 ) ! Runge-Kutta 4th order
        ! Before we begin with rk4, we need to extrapolate jpara and aext
        jparamid = (0.5_dp)*( jparanext + jpara )
        aauxmid = this%applied_vector_potential( time-(0.5_dp)*dt )
        aauxnext = this%applied_vector_potential( time )
        ! Now, we apply Runge Kutta of 4th order
        if ( this%is_total_field_given() ) then
          k1(:,1) = jpara(:) - beta*this%a_tot%components
          k1(:,2) = fourpi*c*pvec(:)
          k2(:,1) = jparamid(:) - beta*aauxmid(:)
          k2(:,2) = fourpi*c*(pvec(:) + (dt/2._dp)*k1(:,1))
          k3(:,1) = jparamid(:) - beta*aauxmid(:)
          k3(:,2) = fourpi*c*(pvec(:) + (dt/2._dp)*k2(:,1))
          k4(:,1) = jparanext(:) - beta*aauxnext(:)
          k4(:,2) = fourpi*c*(pvec(:) + (dt)*k3(:,1))
        else
          k1(:,1) = jpara(:) - beta*( this%a_tot%components )
          k1(:,2) = fourpi*c*pvec(:)
          k2(:,1) = jparamid(:) - beta*(aauxmid(:) + this%a_ind%components + (dt/2._dp)*k1(:,2) )
          k2(:,2) = fourpi*c*(pvec(:) + (dt/2._dp)*k1(:,1))
          k3(:,1) = jparamid(:) - beta*(aauxmid(:) + this%a_ind%components + (dt/2._dp)*k2(:,2) )
          k3(:,2) = fourpi*c*(pvec(:) + (dt/2._dp)*k2(:,1))
          k4(:,1) = jparanext(:) - beta*(aauxnext(:) + this%a_ind%components + (dt)*k3(:,2) )
          k4(:,2) = fourpi*c*(pvec(:) + (dt)*k3(:,1))
        end if
        pvec(:) = pvec(:) + (dt/6._dp)*( k1(:,1) + 2._dp*k2(:,1) + 2._dp*k3(:,1) + k4(:,1) )
        call this%a_ind%add_vector( (dt/6._dp)*( k1(:,2) + 2._dp*k2(:,2) + 2._dp*k3(:,2) + k4(:,2) ) )
      case default
        ! Method not recognized
        ! We need to stop the code
        call terminate( 'Error(Solve_ODE_Vector_Potential): method is not recognized.' )

    end select
  end subroutine

end module rttddft_VectorPotential
