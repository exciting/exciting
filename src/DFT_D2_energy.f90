subroutine DFT_D2_energy(e_disp)
  use modinput
  use DFT_D2_parameters, only : max_elem, cutoff, damping_const, s6, rs6
  use DFT_D2_subroutines, only : loadoldpar, getnumofatoms, getatomdata, getlatticerepetition
  implicit none
  real(8) :: e_disp
  real(8) :: r0ab(max_elem,max_elem)
  real(8) :: c6ab(max_elem,max_elem)
  real(8), dimension(:,:), allocatable :: xyz
  integer, dimension(:), allocatable :: iz
  integer :: numofatoms
  integer :: latrep(3)
  integer :: iat, jat, tau_a, tau_b, tau_c, tau_coeff(3)
  real(8) :: c6, tau(3), dx, dy, dz, r, r6, damp6
  
  call loadoldpar(c6ab,r0ab)
  call getnumofatoms(numofatoms)
  allocate(xyz(3,numofatoms))
  allocate(iz(numofatoms))
  call getatomdata(numofatoms,xyz,iz) 
  call getlatticerepetition(latrep, cutoff)
!!!!!!! calculation of e_disp !!!!!!!
  e_disp = 0
  do iat = 1,numofatoms-1
     do jat = iat+1,numofatoms
        c6 = c6ab(iz(jat),iz(iat))
        do tau_a = -latrep(1),latrep(1)
           do tau_b = -latrep(2),latrep(2)
              do tau_c = -latrep(3),latrep(3)
                 tau_coeff = (/tau_a, tau_b, tau_c/)
                 tau = matmul(input%structure%crystal%basevect,tau_coeff)
                 dx=xyz(1,iat)-xyz(1,jat)+tau(1)
                 dy=xyz(2,iat)-xyz(2,jat)+tau(2)
                 dz=xyz(3,iat)-xyz(3,jat)+tau(3)
                 r=sqrt(dx*dx+dy*dy+dz*dz)
                 if(r.gt.cutoff) cycle
                 damp6=1./(1.+exp(-damping_const*(r/(rs6*r0ab(iz(jat),iz(iat)))-1.)))
                 r6=r**6
                 e_disp =e_disp+c6*damp6/r6
              enddo
           enddo
        enddo
     enddo
  enddo

  do iat = 1,numofatoms
     jat = iat
     c6 = c6ab(iz(jat),iz(iat))
     do tau_a = -latrep(1),latrep(1)
        do tau_b = -latrep(2),latrep(2)
           do tau_c = -latrep(3),latrep(3)
              if (tau_a.eq.0 .and. tau_b.eq.0 .and. tau_c.eq.0) cycle
              tau_coeff = (/tau_a, tau_b, tau_c/)
              tau = matmul(input%structure%crystal%basevect,tau_coeff)
              dx = tau(1)
              dy = tau(2)
              dz = tau(3)
              r = sqrt(dx*dx+dy*dy+dz*dz)
              if(r.gt.cutoff) cycle
              damp6 = 1./(1.+exp(-damping_const*(r/(rs6*r0ab(iz(jat),iz(iat)))-1.)))
              r6 = r**6
              e_disp = e_disp+c6*damp6/r6*0.50d0
           enddo
        enddo
     enddo
  enddo  
  e_disp = -s6*e_disp
end subroutine DFT_D2_energy
