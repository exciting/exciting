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
  use modinput, only: laser_type
  use modmpi, only: terminate
  use precision, only: dp, i32
  use rttddft_laser, only: Set_of_Laser_Pulses
  use rttddft_VectorField, only: Uniform_Vector_Field

  implicit none

  private

  public  :: euler, improved_euler, midpoint, rk4

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
    integer(kind(solver_types))  :: vector_potential_solver
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

end module rttddft_VectorPotential
