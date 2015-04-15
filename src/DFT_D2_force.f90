subroutine DFT_D2_force(force_disp,numofatoms)
  use modinput
  use DFT_D2_parameters, only : max_elem, cutoff, damping_const, s6, rs6
  use DFT_D2_subroutines
  implicit none
  real(8) :: force_disp (3, numofatoms)
  real(8) :: r0ab(max_elem,max_elem)
  real(8) :: c6ab(max_elem,max_elem)
  real(8) :: xyz(3,numofatoms)
  integer :: iz(numofatoms)
  integer :: numofatoms
  integer :: latrep(3)
  integer :: iat, jat, tau_a, tau_b, tau_c, tau_coeff(3)
  real(8) :: c6, tau(3), dx, dy, dz, r, r8, r0ab_, help1_exp, help2

  call loadoldpar(c6ab,r0ab)
!!$  call getnumofatoms(numofatoms)
!!$  allocate(force_disp(3,numofatoms))
!!$  allocate(xyz(3,numofatoms))
!!$  allocate(iz(numofatoms))
  call getatomdata(numofatoms,xyz,iz) 
  call getlatticerepetition(latrep, cutoff)
!!!!!!! calculation of force_disp !!!!!!!
  force_disp = 0
   do iat = 1,numofatoms
      do jat = 1,numofatoms
         if (iat .eq. jat) cycle
         c6 = c6ab(iz(jat),iz(iat))
         do tau_a = -latrep(1),latrep(1)
            do tau_b = -latrep(2),latrep(2)
               do tau_c = -latrep(3),latrep(3)
                  tau_coeff = (/tau_a, tau_b, tau_c/)
                  tau = matmul(input%structure%crystal%basevect,tau_coeff)
                  dx=xyz(1,jat)-xyz(1,iat)+tau(1)
                  dy=xyz(2,jat)-xyz(2,iat)+tau(2)
                  dz=xyz(3,jat)-xyz(3,iat)+tau(3)
                  r=sqrt(dx*dx+dy*dy+dz*dz)
                  if(r.gt.cutoff) cycle
                  !damp6=1./(1.+exp(-damping_const*(r/(rs6*r0ab(iz(jat),iz(iat)))-1.)))
                  r8=r**8
                  r0ab_=rs6*r0ab(iz(jat),iz(iat))
                  help1_exp=exp(damping_const*(r/r0ab_-1.))
                  help2=c6*help1_exp*(damping_const*r-6*r0ab_*(1+help1_exp))/r8/r0ab_/(1+help1_exp)**2
                  force_disp(1,iat) =force_disp(1,iat)+help2*dx
                  force_disp(2,iat) =force_disp(2,iat)+help2*dy
                  force_disp(3,iat) =force_disp(3,iat)+help2*dz
               enddo
            enddo
         enddo
      enddo
   enddo
 
!!$   do iat = 1,numofatoms
!!$      jat = iat
!!$      c6 = c6ab(iz(jat),iz(iat))
!!$      do tau_a = -latrep(1),latrep(1)
!!$         do tau_b = -latrep(2),latrep(2)
!!$            do tau_c = -latrep(3),latrep(3)
!!$               if (tau_a.eq.0 .and. tau_b.eq.0 .and. tau_c.eq.0) cycle
!!$               tau_coeff = (/tau_a, tau_b, tau_c/)
!!$               tau = matmul(input%structure%crystal%basevect,tau_coeff)
!!$               dx = tau(1)
!!$               dy = tau(2)
!!$               dz = tau(3)
!!$               r = sqrt(dx*dx+dy*dy+dz*dz)
!!$               if(r.gt.cutoff) cycle
!!$               damp6 = 1./(1.+exp(-damping_const*(r/(rs6*r0ab(iz(jat),iz(iat)))-1.)))
!!$               r6 = r**6
!!$               force_disp = force_disp+c6*damp6/r6*0.50d0
!!$            enddo
!!$         enddo
!!$      enddo
!!$   enddo  
   force_disp = s6*force_disp
end subroutine DFT_D2_force
