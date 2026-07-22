!! Module to calculate the van-der-Waals contribution due to semiempirical corrections to the Dynamical Matrix

module mod_semiempirical_vdw_dynmat
       use precision, only: i32, dp
       use modmain, only: pi
       implicit none

       private
       public :: vdW_dynmat_TS, vdW_dynmat_DFTD2

 contains

 !! NAME
 !! vdw_dynmat
 !!
 !! FUNCTION
 !! Compute the second derivative of the van der Waals interaction. 
 !! The used semiemperical formalism follows the describtion of Tkatschenko and Scheffler (TSvdW)
 !! DOI: https://doi.org/10.1103/PhysRevLett.102.073005
 !! and Grimme (DFT-D2)
 !! DOI: https://doi.org/10.1002/jcc.20495
 !! INPUTS
 !!     s6, sr6, damping_const, cutoff: parameter to calculate the pairwise vdW-energy
 !!     nq = number of q points
 !!     qvec = list of g-vectors
 !!
 !! OUTPUT
 !!     vdW_dynmat: contribution to dynamical matrix from TS-vdW dispersion potential

 subroutine  vdW_dynmat_TS(s6,sr6,damping_const,cutoff,nq,qvec,dynmat)
       use mod_atoms, only: sprmax, atposc, nspecies, natoms, spzn, idxas, natmtot
       use vdw_general_routines, only: vdw_energy_pairwiseC6, set_default_vdW_parameters, getlatticerepetition
       use TS_vdW_module, only: get_TS_parameters, C6ab, R0_eff_ab 
       use modinput, only: input
       use, intrinsic :: ieee_arithmetic
       implicit none
 !! ==== Input Variables ====
       real(dp), intent(in)         :: s6, sr6, damping_const, cutoff 
       integer, intent(in)          :: nq 
       real(dp), intent(in)         :: qvec(3,nq) 

 !! ==== Output Variables ====
       complex(dp), intent(out)     :: dynmat(3*natmtot, 3*natmtot, nq) 

 !! ==== Local Variables ====
       real(dp)   :: xyz(3, natmtot)
       integer    :: latrep(3)
       integer    :: iat, jat, ia, is, iq, alpha, beta, im, i, j
       integer    :: tau_a, tau_b, tau_c
       real(dp)   :: tau(3), tau_coeff(3)
       real(dp)   :: dx, dy, dz, r
       real(dp)   :: dyn_help
       real(dp)   :: c_denom
       real(dp)   :: r_vec(3)
       complex(dp):: phase 
       real(dp)   :: dotqR
       real(dp)   :: ifc !interatomic force constant
       real(dp), allocatable :: asr_sum(:,:,:) !accoustic sum rule correction
       real(dp)   :: test_block(3,3)
       integer    :: i_idx, j_idx
       real(dp)   :: arg, dyn_help_ratio, term1, term2 
       real(dp) :: term2_part1, term2_part2 
       real(dp) :: xyz_frac(3, natmtot)
       allocate(asr_sum(3,3,natmtot))
       dynmat = (0.0_dp, 0.0_dp) 
       asr_sum = 0.0_dp

       if (.not. allocated(C6ab)) then
            allocate(C6ab(natmtot,natmtot), R0_eff_ab(natmtot,natmtot))
            C6ab = 0.0_dp
            R0_eff_ab = 0.0_dp      
           call get_TS_parameters()       
       end if
 
       ! copy all atomic positions in an (3, natmtot)-array 
       do is = 1, nspecies
            do ia = 1, natoms(is)
                  xyz(:, idxas(ia,is)) = atposc(:, ia, is)
                  im = idxas(ia,is) 
                  xyz_frac(:,im) = input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:)
            end do
       end do 

       call getlatticerepetition(latrep, cutoff)

       dynmat = cmplx(0.0_dp, 0.0_dp, kind=dp)

       !! main loop regarding atom i and atom j
       do iat = 1, natmtot
            do jat = 1, natmtot
                  !! loop over all periodic pictures
                  do tau_a = -latrep(1), latrep(1)
                        do tau_b = -latrep(2), latrep(2)
                              do tau_c = -latrep(3), latrep(3)
                                    !! avoid self-interaction 
                                    if (iat == jat .and. tau_a == 0 .and. tau_b == 0 .and. tau_c == 0) Cycle

                                    tau_coeff = dble((/tau_a, tau_b, tau_c/))
                                    tau = matmul(input%structure%crystal%basevect, tau_coeff)
                                    
                                    !! calculate distances
                                    dx = xyz(1, jat) - xyz(1, iat) + tau(1)
                                    dy = xyz(2, jat) - xyz(2, iat) + tau(2)
                                    dz = xyz(3, jat) - xyz(3, iat) + tau(3)

                                    r_vec(1) = dx
                                    r_vec(2) = dy
                                    r_vec(3) = dz
                                    r = norm2(r_vec)

                                    if ((r .gt. cutoff) .or. (r .lt. 1.0d-6)) cycle

                                    c_denom = sr6 * R0_eff_ab(iat,jat) !constant part of denominator
                                    dyn_help = exp(damping_const * (r / c_denom - 1))

                                    if (c_denom < 1.0d-6) cycle  
                                    arg = damping_const * (r / c_denom - 1.0_dp)
                                    if (arg > 100.0_dp) then
                                          dyn_help_ratio = 1.0_dp
                                    else 
                                          dyn_help = exp(arg)
                                          dyn_help_ratio = dyn_help / (1.0_dp +dyn_help)
                                    end if
           
                                    !! ==== analytical formula of the vdW-contribution to the interatomic force constant ====
                                    term1 = -s6 * C6ab(iat,jat) / (r**10) * (dyn_help_ratio / (1.0_dp + dyn_help)**2) * &
                                          ( 48.0_dp * (1.0_dp + dyn_help)**2 - &
                                          13.0_dp * (damping_const * r / c_denom) * (1.0_dp + dyn_help) - &
                                          (damping_const * r / c_denom)**2 * (dyn_help - 1.0_dp) )

                                    term2_part1 = damping_const * r / (1.0_dp + dyn_help)**2
                                    term2_part2 = 6.0_dp * c_denom / (1.0_dp + dyn_help)     
                                    term2 = s6 * C6ab(iat,jat) * dyn_help * (term2_part1 - term2_part2) / ( (r**8) * c_denom )   

                                          !! loop over cartesian directions
                                          do alpha = 1, 3
                                                do beta = 1, 3
                              
                                                      !! calculate the contribution to the interatomic force constants                                         
                                                      ifc = - term1 * r_vec(alpha) * r_vec(beta)     
                                                      if (alpha == beta) then
                                                          ifc = ifc + term2
                                                      end if
                                          
                                                      !! q=0 correction
                                                      asr_sum(alpha,beta,iat) = asr_sum(alpha,beta,iat) - ifc 

                                                      !! Fourier-Transformation of the interatomic force constant
                                                      do iq = 1, nq
                                                            dotqR = 2.0_dp * pi * ( &
                                                                        qvec(1,iq)* (real(tau_a, dp) + xyz_frac(1, jat) - xyz_frac(1, iat)) + &
                                                                        qvec(2,iq)* (real(tau_b, dp) + xyz_frac(2, jat) - xyz_frac(2, iat)) + &
                                                                        qvec(3,iq)* (real(tau_c, dp) + xyz_frac(3, jat) - xyz_frac(3, iat)) )
                                                            phase = cmplx(cos(dotqR), sin(dotqR), kind=dp)

                                                            dynmat(alpha + (iat -1) * 3, beta + (jat -1) * 3, iq) = &
                                                            dynmat(alpha + (iat -1) * 3, beta + (jat -1) * 3, iq) &
                                                            + ifc * phase 
                                                      end do !iq        

                                                end do !beta
                                          end do !alpha
                              end do !tau_c
                        end do !tau_b
                  end do !tau_a
            end do !jat
       end do !iat

      ! =========================================================================
      ! Acoustic Sum Rule (ASR) Correction
      ! =========================================================================
      ! Enforce the acoustic sum rule to restore translational invariance, which
      ! is slightly broken due to numerical inaccuracies.
      ! This ensures that the acoustic modes correctly vanish at the Gamma point.
      ! =========================================================================
       do iq = 1, nq
            do iat = 1, natmtot
                  do alpha = 1, 3
                        do beta = 1, 3
                              dynmat(alpha + (iat -1) * 3, beta + (iat -1) * 3, iq) = & ! Masses are globally added at the diagonalization-routine
                              dynmat(alpha + (iat -1) * 3, beta + (iat -1) * 3, iq) + asr_sum(alpha, beta, iat) 
                        end do
                  end do
            end do
      end do 

      if (allocated(asr_sum)) deallocate(asr_sum)

      if ( any(isnan(real(dynmat))) .or. any(isnan(aimag(dynmat))) ) then
            print*, "DEBUG: NaNs in dynmat detected!"                        
      end if

      if ( any(.not. ieee_is_finite(real(dynmat))) ) then
            print*, "DEBUG: Matrix contains infinity"                               
      end if

      end subroutine vdW_dynmat_TS

      subroutine  vdW_dynmat_DFTD2(s6,sr6,damping_const,cutoff,nq,qvec,dynmat_DFTD2)
            use mod_atoms, only: sprmax, atposc, nspecies, natoms, spzn, idxas, natmtot
            use vdw_general_routines, only: vdw_energy_pairwiseC6, set_default_vdW_parameters, getlatticerepetition
            Use DFT_D2_module, Only : loadoldpar
            use modinput, only: input
            use, intrinsic :: ieee_arithmetic
            implicit none

            !! ==== Input Variables ====
            real(dp), intent(in)         :: s6, sr6, damping_const, cutoff 
            integer, intent(in)          :: nq 
            real(dp), intent(in)         :: qvec(3,nq) 

             !! ==== Output Variables ====
            complex(dp), intent(out)     :: dynmat_DFTD2(3*natmtot, 3*natmtot, nq) 

            !! ==== Local Variables ====
                  real(dp)   :: xyz(3, natmtot)
                  integer    :: latrep(3)
                  integer    :: iat, jat, ia, is, iq, alpha, beta, im, i, j
                  integer    :: tau_a, tau_b, tau_c
                  real(dp)   :: tau(3), tau_coeff(3)
                  real(dp)   :: dx, dy, dz, r
                  real(dp)   :: dyn_help
                  real(dp)   :: c_denom
                  real(dp)   :: r_vec(3)
                  complex(dp):: phase 
                  real(dp)   :: dotqR
                  real(dp)   :: ifc !interatomic force constant
                  real(dp), allocatable :: asr_sum(:,:,:) !accoustic sum rule correction
                  real(dp)   :: test_block(3,3)
                  integer    :: i_idx, j_idx
                  real(dp)   :: arg, dyn_help_ratio, term1, term2 !! NEW
                  real(dp) :: term2_part1, term2_part2    
                  real(dp), allocatable :: C6ab(:,:), R0_eff_ab(:,:)        
                  real(dp) :: xyz_frac(3, natmtot)
                  allocate(asr_sum(3,3,natmtot))
                  dynmat_DFTD2 = (0.0_dp, 0.0_dp) 
                  asr_sum = 0.0_dp
                  
                  if (.not. allocated(C6ab)) then
                        allocate(C6ab(natmtot,natmtot), R0_eff_ab(natmtot,natmtot))
                        C6ab = 0.0_dp
                        R0_eff_ab = 0.0_dp      
                  end if

                  if ( maxval(R0_eff_ab) < 1.0d-8 ) then
                        call loadoldpar(C6ab, R0_eff_ab) 
                   end if

                   do is = 1, nspecies
                        do ia = 1, natoms(is)
                              xyz(:, idxas(ia,is)) = atposc(:, ia, is)
                              im = idxas(ia,is) 
                              xyz_frac(:,im) = input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:)
                        end do
                   end do 
                   call getlatticerepetition(latrep, cutoff)
                   dynmat_DFTD2 = cmplx(0.0_dp, 0.0_dp, kind=dp)
              
       !! main loop regarding atom i and atom j
                   do iat = 1, natmtot
                        do jat = 1, natmtot
            
                        !! loop over all periodic pictures
                              do tau_a = -latrep(1), latrep(1)
                                    do tau_b = -latrep(2), latrep(2)
                                          do tau_c = -latrep(3), latrep(3)
                                                !! avoid self-interaction 
                                                if (iat == jat .and. tau_a == 0 .and. tau_b == 0 .and. tau_c == 0) Cycle

                                                tau_coeff = dble((/tau_a, tau_b, tau_c/))
                                                tau = matmul(input%structure%crystal%basevect, tau_coeff)
                                                
                                                !! calculate distances
                                                dx = xyz(1, jat) - xyz(1, iat) + tau(1)
                                                dy = xyz(2, jat) - xyz(2, iat) + tau(2)
                                                dz = xyz(3, jat) - xyz(3, iat) + tau(3)
            
                                                r_vec(1) = dx
                                                r_vec(2) = dy
                                                r_vec(3) = dz
                                                r = norm2(r_vec)

                                                if ((r .gt. cutoff) .or. (r .lt. 1.0d-6)) cycle

                                                c_denom = sr6 * R0_eff_ab(iat,jat) !constant part of denominator

                                                dyn_help = exp(damping_const * (r / c_denom - 1))
            
                                                if (c_denom < 1.0d-6) cycle 
                                                arg = damping_const * (r / c_denom - 1.0_dp)
                                                if (arg > 100.0_dp) then
                                                      dyn_help_ratio = 1.0_dp
                                                else 
                                                      dyn_help = exp(arg)
                                                      dyn_help_ratio = dyn_help / (1.0_dp +dyn_help)
                                                end if

                                                term1 = -s6 * C6ab(iat,jat) / (r**10) * (dyn_help_ratio / (1.0_dp + dyn_help)**2) * &
                                                ( 48.0_dp * (1.0_dp + dyn_help)**2 - &
                                                13.0_dp * (damping_const * r / c_denom) * (1.0_dp + dyn_help) - &
                                                (damping_const * r / c_denom)**2 * (dyn_help - 1.0_dp) ) 
      
                                                term2_part1 = damping_const * r / (1.0_dp + dyn_help)**2
                                                term2_part2 = 6.0_dp * c_denom / (1.0_dp + dyn_help)     
                                                term2 = s6 * C6ab(iat,jat) * dyn_help * (term2_part1 - term2_part2) / ( (r**8) * c_denom ) 
                                                
                                                do alpha = 1, 3
                                                      do beta = 1, 3
                                                                  
                                                            ifc = - term1 * r_vec(alpha) * r_vec(beta)  
                                                            if (alpha == beta) then
                                                                ifc = ifc + term2 
                                                            end if
      
                                                            asr_sum(alpha,beta,iat) = asr_sum(alpha,beta,iat) - ifc
                                                            do iq = 1, nq
                                                                  dotqR = 2.0_dp * pi * ( &
                                                                        qvec(1,iq)* (real(tau_a, dp) + xyz_frac(1, jat) - xyz_frac(1, iat)) + &
                                                                        qvec(2,iq)* (real(tau_b, dp) + xyz_frac(2, jat) - xyz_frac(2, iat)) + &
                                                                        qvec(3,iq)* (real(tau_c, dp) + xyz_frac(3, jat) - xyz_frac(3, iat)) )

                                                                  phase = cmplx(cos(dotqR), sin(dotqR), kind=dp)
                                                                  dynmat_DFTD2(alpha + (iat -1) * 3, beta + (jat -1) * 3, iq) = &
                                                                  dynmat_DFTD2(alpha + (iat -1) * 3, beta + (jat -1) * 3, iq) &
                                                                  + ifc * phase
                                                            end do !iq        
                                                      end do !beta
                                                end do !alpha
                              end do !tau_c
                        end do !tau_b
                  end do !tau_a
            end do !jat
       end do !iat
       do iq = 1, nq
            do iat = 1, natmtot
                  do alpha = 1, 3
                        do beta = 1, 3
                              dynmat_DFTD2(alpha + (iat -1) * 3, beta + (iat -1) * 3, iq) = & 
                              dynmat_DFTD2(alpha + (iat -1) * 3, beta + (iat -1) * 3, iq) + asr_sum(alpha, beta, iat) 
                        end do
                  end do
            end do
      end do 
      if (allocated(asr_sum)) deallocate(asr_sum)

      end subroutine vdW_dynmat_DFTD2

end module mod_semiempirical_vdw_dynmat 
 