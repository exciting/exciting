! Copyright (C) 2013 exciting team
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: mt_pot
! !INTERFACE:
!
!
Subroutine mt_pot(pot,basis,mt_h,bxcmt_case)
  ! !USES:
  Use modinput
  Use modmain
  Use mod_compute_muffin_tin_potential, only: integrate_mt_potential, build_mt_potential_matrix_element
  Use mod_lattice_harmonics, only: transform_mt_potential_lattice_harmonics, get_lattice_harmonics, &
       lattice_harmonics_type
  ! !DESCRIPTION:
  !   Calculates the potential energy contribution to the muffin-tin Hamiltonian.
  !
  !EOP
  !BOC
  Implicit None
  Real(8), intent(in) :: pot(lmmaxvr,nrmtmax,natmtot)
  type(apw_lo_basis_type) :: basis
  Type (MTHamiltonianList) :: mt_h
  ! Checks if the routine is used to compute the MT potential or the MT XC magnetic field.
  logical, intent(in) :: bxcmt_case

  ! local variables
  Integer :: is, ia, ias, nr, ir, if1,if3
  Integer :: l1, l2, l3, m2, lm2, m1, m3, lm1, lm3
  Integer :: ilo, ilo1, ilo2, io, io1, io2, maxnlo, maxaa, num_sph_harm
  Real (8) :: t1,t2,angular
  Real (8), allocatable :: haaintegrals(:,:,:,:,:),halointegrals(:,:,:,:),hlolointegrals(:,:,:), &
       radial_integral(:), mt_potential(:,:,:)
  complex (8) :: zsum
  complex (8), allocatable :: gaunt_ryy_general(:,:,:,:)
  type(lattice_harmonics_type) :: lattice_harmonics

  ! Check if lattice harmonics (symmetrized spherical harmonics) should be used in the computation of the MT potential
  ! and, if so, transform the potential coefficients and the gaunt coefficients to the lattice-harmonic representation.
  if (input%groundstate%LatticeHarmonics .and. .not. bxcmt_case) then
     lattice_harmonics = get_lattice_harmonics()
     num_sph_harm = maxval(sum(lattice_harmonics%number, dim=1)) + 1

     allocate(mt_potential(num_sph_harm, nrmtmax, natmtot))
     call transform_mt_potential_lattice_harmonics(pot, nspecies, natoms, idxas, nrmt, nrmtmax, natmtot, &
          lattice_harmonics%coefficients, lattice_harmonics%number, idxlm, mt_potential)

     allocate(gaunt_ryy_general(natmtot, num_sph_harm, lmmaxapw, lmmaxapw))
     gaunt_ryy_general = lattice_harmonics%gaunt_coefficients
  else
     num_sph_harm = lmmaxvr

     allocate(mt_potential(num_sph_harm, nrmtmax, natmtot))
     mt_potential = pot

     allocate(gaunt_ryy_general(natmtot, num_sph_harm, lmmaxapw, lmmaxapw))
     do ias = 1, natmtot
        gaunt_ryy_general(ias, :, :, :) = gntryy
     end do
  end if

  ! APW-APW storage initialisation
  haaijSize=0
  Do is = 1, nspecies
     if1=0
     Do l1 = 0, input%groundstate%lmaxmat
        Do m1 = - l1, l1
           lm1 = idxlm (l1, m1)
           Do io1 = 1, apword (l1, is)
              if1=if1+1
           End Do
        End Do
     End Do
     if (if1.gt.haaijSize) haaijSize=if1
  Enddo

  Allocate (haaintegrals(num_sph_harm, apwordmax, 0:input%groundstate%lmaxapw, apwordmax, 0:input%groundstate%lmaxmat))
  haaintegrals (:, :, :, :, :)=1d100
  Allocate (radial_integral(num_sph_harm))

  maxnlo=mt_h%maxnlo
  if (maxnlo.gt.0) then
     Allocate (halointegrals(num_sph_harm, apwordmax, 0:input%groundstate%lmaxmat, nlomax))
     ! LO-LO storage initialisation
     allocate(hlolointegrals(num_sph_harm,nlomax,nlomax))
  endif

  ! begin loops over atoms and species
  Do is = 1, nspecies
     nr = nrmt (is)
     Do ia = 1, natoms (is)
        ias = idxas (ia, is)
        !---------------------------!
        !     APW-APW integrals     !
        !---------------------------!
        ! Radial integrals first
#ifdef USEOMP
        !$OMP PARALLEL DEFAULT(NONE) SHARED(lorbl,nlorb,input,apword,num_sph_harm,apwfr,lofr,mt_potential,spr,nr,haaintegrals,hlolointegrals,halointegrals,is,ias) PRIVATE(lm2,m2,l2,ir,t1,l1,l3,t2,angular,io1,io2,ilo1,ilo2,io,ilo,radial_integral)
#endif
        Do l1 = 0, input%groundstate%lmaxmat
           Do io1 = 1, apword (l1, is)
              Do l3 = 0, input%groundstate%lmaxmat
                 Do io2 = 1, apword (l3, is)
                    call integrate_mt_potential(num_sph_harm, apwfr(:, 1, io1, l1, ias), &
                         apwfr(:, 1, io2, l3, ias), nr, mt_potential(1:num_sph_harm, 1:nr, ias), spr(:, is), radial_integral)
                    haaintegrals(:, io2, l3, io1, l1) = radial_integral
                 End Do
              End Do
           End Do
        End Do

        !--------------------------------------!
        !     local-orbital-APW integrals      !
        !--------------------------------------!
        Do ilo = 1, nlorb (is)
           l1 = lorbl (ilo, is)
           Do l3 = 0, input%groundstate%lmaxmat
              Do io = 1, apword (l3, is)
                 call integrate_mt_potential(num_sph_harm, lofr(:, 1, ilo, ias), &
                      apwfr(:, 1, io, l3, ias), nr, mt_potential(1:num_sph_harm, 1:nr, ias), spr(:, is), radial_integral)
                 halointegrals(:, io, l3, ilo) = radial_integral
              End Do
           End Do
        End Do

        !-----------------------------------------------!
        !     local-orbital-local-orbital integrals     !
        !-----------------------------------------------!
        Do ilo1 = 1, nlorb (is)
           l1 = lorbl (ilo1, is)
           Do ilo2 = 1, nlorb (is)
              l3 = lorbl (ilo2, is)
              call integrate_mt_potential(num_sph_harm, lofr(:, 1, ilo1, ias), &
                   lofr(:, 1, ilo2, ias), nr, mt_potential(1:num_sph_harm, 1:nr, ias), spr(:, is), radial_integral)
              hlolointegrals(:, ilo1, ilo2) = radial_integral
           End Do
        End Do

#ifdef USEOMP
        !$OMP END PARALLEL
#endif

        ! Now the angular integrals
        t1 = 0.5d0 * rmt (is) ** 2
        if1=0
        Do l1 = 0, input%groundstate%lmaxmat
           Do m1 = - l1, l1
              lm1 = idxlm (l1, m1)
              Do io1 = 1, apword (l1, is)
                 if1=if1+1
                 if3=0
                 Do l3 = 0, input%groundstate%lmaxmat
                    Do m3 = - l3, l3
                       lm3 = idxlm (l3, m3)
                       Do io2 = 1, apword (l3, is)
                          if3=if3+1
                          call build_mt_potential_matrix_element(num_sph_harm, gaunt_ryy_general(ias, :, lm3, lm1), &
                               haaintegrals (:, io2, l3, io1, l1), zsum)
                          mt_h%main%aa(if1,if3,ias)=mt_h%main%aa(if1,if3,ias)+zsum
                       End Do
                    End Do
                 End Do
              End Do
           End Do
        End Do

        if1=0
        Do ilo = 1, nlorb (is)
           l1 = lorbl (ilo, is)
           Do m1 = - l1, l1
              lm1 = idxlm (l1, m1)
              if1=if1+1
              if3=0
              Do l3 = 0, input%groundstate%lmaxmat
                 Do m3 = - l3, l3
                    lm3 = idxlm (l3, m3)
                    Do io = 1, apword (l3, is)
                       if3=if3+1
                       call build_mt_potential_matrix_element(num_sph_harm, gaunt_ryy_general(ias, :, lm3, lm1), &
                            halointegrals(:, io, l3, ilo), zsum)
                       mt_h%main%loa(if1,if3,ias)=mt_h%main%loa(if1,if3,ias)+zsum
                       mt_h%main%alo(if3,if1,ias)=mt_h%main%loa(if1,if3,ias)
                    End Do
                 End Do
              End Do
           End Do
        End Do

        if1=0
        Do ilo1 = 1, nlorb (is)
           l1 = lorbl (ilo1, is)
           Do m1 = - l1, l1
              lm1 = idxlm (l1, m1)
              if1=if1+1
              if3=0
              Do ilo2 = 1, nlorb (is)
                 l3 = lorbl (ilo2, is)
                 Do m3 = - l3, l3
                    lm3 = idxlm (l3, m3)
                    if3=if3+1
                    call build_mt_potential_matrix_element(num_sph_harm, gaunt_ryy_general(ias, :, lm3, lm1), &
                         hlolointegrals(:,ilo1,ilo2), zsum)
                    mt_h%main%lolo(if1,if3,ias)=mt_h%main%lolo(if1,if3,ias)+zsum
                 End Do
              End Do
           End Do
        End Do
        ! end loops over atoms and species
     End Do
  End Do
  ! cleaning up
  deallocate(haaintegrals)
  if (allocated(halointegrals)) deallocate(halointegrals)
  if (allocated(hlolointegrals)) deallocate(hlolointegrals)
  Return
End Subroutine mt_pot
!EOC
