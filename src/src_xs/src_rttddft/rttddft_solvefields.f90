module rttddft_solve_fields

  use constants, only: fourpi, pi
  use modmpi, only: terminate
  use mod_lattice, only: omega
  use physical_constants, only: c
  use precision, only: dp
  use rttddft_CurrentDensity, only: Current_Density, Current_Density_Field
  use rttddft_Polarization, only: Polarization
  use rttddft_VectorPotential, only: Vector_Potential, euler, improved_euler, midpoint, rk4

  implicit none

  private

  public :: update_a_ind_and_p_vec

contains 

!> Update the induced vector potential \( \mathbf{A}_{ind} \)
!> using the current density \( \mathbf{J}(t) \)
!> We need to solve the differential equation
!> \[
!>  \frac{d^2\mathbf{A}_{ind}}{dt^2} = 4 \pi c \mathbf{J}(t).
!>  \]
subroutine update_a_ind_and_p_vec( t, dt, j_t_minus_dt, j_para_t, vec_pot, p_vec, active_charge )
  !> Time \( t \)
  real(dp), intent(in) :: t
  !> Time step \( \Delta t \)
  real(dp), intent(in) :: dt
  !> Current density at time \(t-\Delta t\)
  class(Current_Density), intent(in) :: j_t_minus_dt
  !> Current density at time \(t\)
  class(Current_Density_Field), intent(in) :: j_para_t
  !> In: vector potential at time \(t-\Delta t\)
  !> Out: vector potential with `a_ind` at time \(t\)
  class(Vector_Potential), intent(inout) :: vec_pot
  !> In: Polarization vector at time \(t-\Delta t\)
  !> Out: Polarization vector at time \(t\)
  class(Polarization), intent(inout) :: p_vec
  !> Total charge of active electrons
  real(dp), intent(in) :: active_charge

  real(dp) :: beta, fac, den, k1(3,2), k2(3,2), k3(3,2), k4(3,2), j_para_mid(3), &
    j_ind_mid(3), j_ind_t(3), a_mid(3), smid(3), a_applied_t(3), a_save(3), a_ext_t_minus_t(3)

  beta = active_charge / c / omega
  select case( vec_pot%vector_potential_solver )
    case( euler ) ! Euler
      call vec_pot%a_ind%add_vector( fourpi*c*dt*p_vec%components )
      call p_vec%add_vector( dt*j_t_minus_dt%total_components() )
    case( improved_euler ) ! Improved Euler method
      call vec_pot%a_ind%add_vector( fourpi*c*dt*(p_vec%components + (0.5_dp)*(dt)*j_t_minus_dt%total_components()) )
      j_ind_t = j_para_t%components - beta*( vec_pot%applied_vector_potential( t ) )
      if ( .not. vec_pot%is_total_field_given() ) j_ind_t = j_ind_t - beta*( vec_pot%a_ind%components )
      call p_vec%add_vector( 0.5_dp*dt*( j_t_minus_dt%total_components() + j_ind_t ) )
    case( midpoint )
      a_applied_t = vec_pot%applied_vector_potential( t )
      if ( vec_pot%is_total_field_given() ) then
        j_ind_t = j_para_t%components - beta*( a_applied_t )
        j_ind_mid = 0.5_dp*( j_t_minus_dt%total_components() + j_ind_t )
        call vec_pot%a_ind%add_vector( fourpi*c*dt*( p_vec%components + 0.5_dp*dt*j_ind_mid ) )
        call p_vec%add_vector( dt*j_ind_mid )
      else
        j_para_mid = (0.5_dp)*( j_para_t%components + j_t_minus_dt%paramagnetic%components )
        a_applied_t = vec_pot%applied_vector_potential( t )
        a_ext_t_minus_t = vec_pot%a_tot%components - vec_pot%a_ind%components
        smid = j_para_mid - 0.5_dp*beta*( a_applied_t + a_ext_t_minus_t )
        fac = pi*beta*c*(dt**2)
        den = 1_dp + fac
        fac = (1_dp - fac)/den
        a_save = vec_pot%a_ind%components
        vec_pot%a_ind%components = (fourpi*c*dt/den)*( p_vec%components + 0.5_dp*dt*smid ) + fac*vec_pot%a_ind%components
        p_vec%components = (dt/den)*(smid - beta*a_save ) + fac*p_vec%components
      end if
    case( rk4 ) ! Runge-Kutta 4th order
      ! Before we begin with rk4, we need to extrapolate jpara and aext
      j_para_mid = (0.5_dp)*( j_para_t%components + j_t_minus_dt%paramagnetic%components )
      a_mid = vec_pot%applied_vector_potential( t-0.5_dp*dt )
      a_applied_t = vec_pot%applied_vector_potential( t )
      k1(:,1) = j_t_minus_dt%paramagnetic%components - beta*vec_pot%a_tot%components
      k1(:,2) = fourpi*c*p_vec%components
      ! Now, we apply Runge Kutta of 4th order
      if ( vec_pot%is_total_field_given() ) then
        k2(:,1) = j_para_mid - beta*a_mid
        k2(:,2) = fourpi*c*(p_vec%components + (dt/2._dp)*k1(:,1))
        k3(:,1) = j_para_mid - beta*a_mid
        k3(:,2) = fourpi*c*(p_vec%components + (dt/2._dp)*k2(:,1))
        k4(:,1) = j_para_t%components - beta*a_applied_t
        k4(:,2) = fourpi*c*(p_vec%components + (dt)*k3(:,1))
      else
        k2(:,1) = j_para_mid - beta*(a_mid(:) + vec_pot%a_ind%components + (dt/2._dp)*k1(:,2) )
        k2(:,2) = fourpi*c*(p_vec%components + (dt/2._dp)*k1(:,1))
        k3(:,1) = j_para_mid - beta*(a_mid(:) + vec_pot%a_ind%components + (dt/2._dp)*k2(:,2) )
        k3(:,2) = fourpi*c*(p_vec%components + (dt/2._dp)*k2(:,1))
        k4(:,1) = j_para_t%components - beta*(a_applied_t + vec_pot%a_ind%components + (dt)*k3(:,2) )
        k4(:,2) = fourpi*c*(p_vec%components + (dt)*k3(:,1))
      end if
      call p_vec%add_vector( (dt/6._dp)*( k1(:,1) + 2._dp*k2(:,1) + 2._dp*k3(:,1) + k4(:,1) ) )
      call vec_pot%a_ind%add_vector( (dt/6._dp)*( k1(:,2) + 2._dp*k2(:,2) + 2._dp*k3(:,2) + k4(:,2) ) )
    case default
      call terminate( 'Error(Solve_ODE_Vector_Potential): method given in &
        & input%xs%rt_tddft%updateAind is not recognized.' )
  end select

end subroutine

end module 