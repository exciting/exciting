!
!
!
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: zpotcoul
! !INTERFACE:
!
!
Subroutine zpotcoul2 (nr, nrmax, ld, r, igp0, gpc, jlgpr, ylmgp, sfacgp, &
& zn, zrhomt, zrhoir, hybrid, zvclmt, zvclir, zrho0)
      Use modinput
      Use mod_kpoint, only: nkptnr
! !USES:
      Use modmain
#ifdef USEOMP
      use omp_lib
#endif
#ifdef PSOLVER
   use dictionaries
   use Poisson_Solver
   use at_domain
   use numerics, only: onehalf
#endif
! !INPUT/OUTPUT PARAMETERS:
!   nr     : number of radial points for each species (in,integer(nspecies))
!   nrmax  : maximum nr over all species (in,integer)
!   ld     : leading dimension of r (in,integer)
!   r      : radial mesh for each species (in,real(ld,nspecies))
!   igp0   : index of the shortest G+p-vector (in,integer)
!   gpc    : G+p-vector lengths (in,real(ngvec))
!   jlgpr  : spherical Bessel functions for evergy G+p-vector and muffin-tin
!            radius (in,real(0:lmaxvr+npsden+1,ngvec,nspecies))
!   ylmgp  : spherical harmonics of the G+p-vectors (in,complex(lmmaxvr,ngvec))
!   sfacgp : structure factors of the G+p-vectors (in,complex(ngvec,natmtot))
!   zn     : nuclear charges at the atomic centers (in,real(nspecies))
!   zrhomt : muffin-tin charge density (in,complex(lmmaxvr,nrmax,natmtot))
!   zrhoir : interstitial charge density (in,complex(ngrtot))
!   hybrid : hybrid calculation option 0 or 1 (in, integer)
!   zvclmt : muffin-tin Coulomb potential (out,complex(lmmaxvr,nrmax,natmtot))
!   zvclir : interstitial Coulomb potential (out,complex(ngrtot))
!   zrho0  : G+p=0 term of the pseudocharge density (out,complex)
! !DESCRIPTION:
!   Calculates the Coulomb potential of a complex charge density by solving
!   Poisson's equation. First, the multipole moments of the muffin-tin charge
!   are determined for the $j$th atom of the $i$th species by
!   $$ q_{ij;lm}^{\rm MT}=\int_0^{R_i}r^{l+2}\rho_{ij;lm}(r)dr+z_{ij}Y_{00}
!   \,\delta_{l,0}\;, $$
!   where $R_i$ is the muffin-tin radius and $z_{ij}$ is a point charge located
!   at the atom center (usually the nuclear charge, which should be taken as
!   {\bf negative}). Next, the multipole moments of the continuation of the
!   interstitial density, $\rho^{\rm I}$, into the muffin-tin are found with
!   $$ q_{ij;lm}^{\rm I}=4\pi i^l R_i^{l+3}\sum_{\bf G}\frac{j_{l+1}(GR_i)}
!    {GR_i}\rho^{\rm I}({\bf G})\exp(i{\bf G}\cdot{\bf r}_{ij})Y_{lm}^*
!    (\hat{\bf G}), $$
!   remembering that
!   $$ \lim_{x\rightarrow 0}\frac{j_{l+n}(x)}{x^n}=\frac{1}{(2n+1)!!}
!    \delta_{l,0} $$
!   should be used for the case ${\bf G}=0$. A pseudocharge is now constructed
!   which is equal to the real density in the interstitial region and whose
!   multipoles are the difference between the real and interstitial muffin-tin
!   multipoles. This pseudocharge density is smooth in the sense that it can be
!   expanded in terms of the finite set of ${\bf G}$-vectors. In each muffin-tin
!   the pseudocharge has the form
!   $$ \rho_{ij}^{\rm P}({\bf r})=\rho^{\rm I}({\bf r}-{\bf r}_{ij})+\sum_{lm}
!    \rho_{ij;lm}^{\rm P}\frac{1}{R_i^{l+3}}\left(\frac{r}{R_i}\right)^l\left(1-
!    \frac{r^2}{R_i^2}\right)^{N_i}Y_{lm}(\hat{\bf r}) $$
!   where
!   $$ \rho_{ij;lm}^{\rm P}=\frac{(2l+2N_i+3)!!}{2^N_iN_i!(2l+1)!!}\left(
!    q_{ij;lm}^{\rm MT}-q_{ij;lm}^{\rm I}\right) $$
!   and $N_i\approx\frac{1}{2}R_iG_{\rm max}$ is generally a good choice.
!   The pseudocharge in reciprocal-space is given by
!   $$ \rho^{\rm P}({\bf G})=\rho^{\rm I}({\bf G})+\sum_{ij;lm}2^{N_i}N_i!
!    \frac{4\pi(-i)^l}{\Omega R_i^l}\frac{j_{l+N_i+1}(GR_i)}{(GR_i)^{N_i+1}}
!    \rho_{ij;lm}^{\rm P}\exp(-i{\bf G}\cdot{\bf r}_{ij})Y_{lm}(\hat{\bf G}) $$
!   which may be used for solving Poisson's equation directly
!   $$ V^{\rm P}({\bf G})=\begin{cases}
!     4\pi\frac{\rho^{\rm P}({\bf G})}{G^2} & G>0 \\
!     0 & G=0 \end{cases}\;. $$
!   The usual Green's function approach is then employed to determine the
!   potential in the muffin-tin sphere due to charge in the sphere. In other
!   words
!   $$ V_{ij;lm}^{\rm MT}(r)=\frac{4\pi}{2l+1}\left(\frac{1}{r^{l+1}}\int_0^r
!    \rho_{ij;lm}^{\rm MT}(r'){r'}^{l+2}dr'+r^l\int_r^{R_i}\frac{
!    \rho_{ij;lm}^{\rm MT}(r')}{{r'}^{l-1}}dr'\right)+\frac{1}{Y_{00}}
!    \frac{z_{ij}}{r}\delta_{l,0} $$
!   where the last term is the monopole arising from the point charge. All that
!   remains is to add the homogenous solution of Poisson's equation,
!   $$ V_{ij}^{\rm H}({\bf r})=\sum_{lm}V_{ij;lm}^{\rm H}\left(\frac{r}
!    {R_i}\right)^lY_{lm}(\hat{\bf r}), $$
!   to the muffin-tin potential so that it is continuous at the muffin-tin
!   boundary. Therefore the coefficients, $\rho_{ij;lm}^{\rm H}$, are given by
!   $$ V_{ij;lm}^{\rm H}=4\pi i^l\sum_{\bf G}j_{l}(Gr)V^{\rm P}({\bf G})
!    \exp(i{\bf G}\cdot{\bf r}_{ij})Y_{lm}^*(\hat{\bf G})-V_{ij;lm}^{\rm MT}
!    (R_i). $$
!   Finally note that the ${\bf G}$-vectors passed to the routine can represent
!   vectors with a non-zero offset, ${\bf G}+{\bf p}$ say, which is required for
!   calculating Coulomb matrix elements.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
      Implicit None
! arguments
      Integer, Intent (In) :: nr (nspecies)
      Integer, Intent (In) :: nrmax
      Integer, Intent (In) :: ld
      Real (8), Intent (In) :: r (ld, nspecies)
      Integer, Intent (In) :: igp0
      Real (8), Intent (In) :: gpc (ngvec)
      Real (8), Intent (In) :: jlgpr(0:input%groundstate%lmaxvr+ &
      & input%groundstate%npsden+1, ngvec, nspecies)
      Complex (8), Intent (In) :: ylmgp (lmmaxvr, ngvec)
      Complex (8), Intent (In) :: sfacgp (ngvec, natmtot)
      Real (8), Intent (In) :: zn (nspecies)
      Complex (8), Intent (In) :: zrhomt (lmmaxvr, nrmax, natmtot)
      Complex (8), Intent (In) :: zrhoir (ngrtot)
      Integer , Intent(In) :: hybrid
      Complex (8), Intent (Out) :: zvclmt (lmmaxvr, nrmax, natmtot)
      Complex (8), Intent (Out) :: zvclir (ngrtot)
      Complex (8), Intent (Out) :: zrho0
! local variables
      Integer :: is, ia, ias, l, m, lm
      Integer :: ir, ig, ifg
      Real (8) :: fpo, t1, t2, t3
      Complex (8) zsum, zt1, zt2
! automatic arrays
      Real (8) :: rmtl (0:input%groundstate%lmaxvr+3, nspecies)
      Real (8) :: rl (nrmax, 0:input%groundstate%lmaxvr)
      Complex (8) vilm (lmmaxvr)
      Complex (8) qmt (lmmaxvr, natmtot)
      Complex (8) qi (lmmaxvr, natmtot)
      Complex (8) qilocal (lmmaxvr)
      Complex (8) zrp (lmmaxvr)
      real(8) :: vn(nrmax)
      real(8) :: third 
      parameter (third=0.3333333333333333333333d0)
#ifdef USEOMP
      integer ithr,nthreads,whichthread
#endif
      real(8), allocatable :: vdplmt(:,:,:), vdplir(:) 

! external functions
      Real (8) :: factnm
      External factnm
!for compiling without PSOLVER library
   logical psolver0d, psolver1d, psolver2d, psolver3d, exciting, exciting0d, excite
   real(8) :: r_c
#ifdef PSOLVER

! PSolver related
   character(len=1) :: solvertype
   character(len=64) :: chain
   integer :: i1,i2,i3,n1,n2,n3,i1_max,i2_max,i3_max,isf_order, i, bc_type
   real(kind=8) :: max_diff, alpha_bc, beta_ac, gamma_ab, bc, ac, ab
   real(kind=8) :: sigma,length,hgrid,mu,energy,offset,acell,epot
   real(kind=8), dimension(:), allocatable :: fake_arr
   real(kind=8), allocatable :: psi, r_v(:), c_v(:)
   type(coulomb_operator) :: kernel
   type(dictionary), pointer :: dict
   type(domain) :: dom
   integer, dimension(3) :: ndims
   real(kind=8), dimension(3) :: hgrids, acel
   real(8), dimension(3,3) :: n_avec
   
#endif

! Initialise PSolver
   psolver0d=(input%groundstate%vha.eq."psolver0d")
   psolver1d=(input%groundstate%vha.eq."psolver1d")
   psolver2d=(input%groundstate%vha.eq."psolver2d")
   psolver3d=(input%groundstate%vha.eq."psolver3d")
   exciting =(input%groundstate%vha.eq."exciting")
   exciting0d=(input%groundstate%vha.eq."exciting0d")
   if (exciting.or.exciting0d) then
      excite = .True.
   end if

   if (hybrid.eq.1) then! for hybrid calculations, using methods for Fock operator calculations
     exciting0d=(input%groundstate%hybrid%singularity.eq."exc0d")
     psolver0d=(input%groundstate%hybrid%singularity.eq."ps0d")
     if (psolver0d) then
            excite = .False.
     end if
   end if
   

#ifdef PSOLVER
if (.not.excite) then
!call f_lib_initialize()

   n1=ngrid(1)
   n2=ngrid(2)
   n3=ngrid(3)

isf_order = 16

   if (psolver0d) then! free in all directions
       !write(*,*)"type Free"
       solvertype='F'
       bc_type = 0
   else if ((psolver3d)) then! periodic x, y, z
       solvertype='P'
       bc_type = 3
       !write(*,*)"type P"
   else if (psolver1d) then!free in x,y, periodic in z
       solvertype = 'W'
       bc_type = 1
       isf_order=60
       !write(*,*) "type W"
   else if (psolver2d) then! free in y, periodic x, z
       solvertype = 'S'
       bc_type = 2
       isf_order=60
       !write(*,*) "type S"
   end if 

   




do i=1,3 
   hgrids(i) = sqrt(dot_product(avec(:,i), avec(:,i)))/ngrid(i)
   
   n_avec(:, i) = avec(:,i)/sqrt(dot_product(avec(:,i), avec(:,i)))
   !write(*,*)i, hgrids(i)
end do
  

!   bc = acos(DOT_PRODUCT(avec(:,2), avec(:,3)) / (SQRT(DOT_PRODUCT(avec(:,2), avec(:,2))) * SQRT(DOT_PRODUCT(avec(:,3), avec(:,3)) ) ) )!y, z
!   ac = acos(DOT_PRODUCT(avec(:,1), avec(:,3)) / (SQRT(DOT_PRODUCT(avec(:,1), avec(:,1))) * SQRT(DOT_PRODUCT(avec(:,3), avec(:,3)) ) ) )!x, z
!   ab = acos(DOT_PRODUCT(avec(:,1), avec(:,2)) / (SQRT(DOT_PRODUCT(avec(:,1), avec(:,1))) * SQRT(DOT_PRODUCT(avec(:,2), avec(:,2)) ) ) )!z, y
!check angles
!write(*,*)  "bc23", bc*57.3
!write(*,*)  "ac13", ac*57.3
!write(*,*)  "ab12", ab*57.3



   dom=domain_new(units=ATOMIC_UNITS,bc=geocode_to_bc_enum(solvertype), abc = n_avec)!alpha_bc=bc, beta_ac=ac, gamma_ab=ab, acell = ngrid*hgrids
   dict=>dict_new('kernel' .is. dict_new('isf_order' .is. isf_order))
   kernel=pkernel_init(0,1,dict,dom,ngrid,hgrids)
   call dict_free(dict)
   call pkernel_set(kernel,verbose=.true.)




end if

#endif


      fpo = fourpi / omega
! solve Poisson's equation for the isolated charge in the muffin-tin
      Do is = 1, nspecies
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ias)
!$OMP DO
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            Call zpotclmt (input%groundstate%ptnucl, &
            & input%groundstate%lmaxvr, nr(is), r(:, is), zn(is), &
            & lmmaxvr, zrhomt(:, :, ias), zvclmt(:, :, ias))
            ias = idxas (ia, is)
         End Do
!$OMP END DO
!$OMP END PARALLEL
      End Do

! add nuclear contributions
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
           vmad(ias)=zvclmt (1, 1,ias)*y00
           If (zn(is) .Ne. 0.d0) Then
             Call potnucl (input%groundstate%ptnucl, nr(is), r(:,is), zn(is), vn)
             t1 = 1.d0 / y00
             Do ir = 1, nr(is)
               zvclmt (1, ir,ias) = zvclmt (1, ir,ias) + t1 * vn (ir)
             End Do
           End If
         End Do
      Enddo

! add dipole correction
      !if ((iscl>0).and.(input%groundstate%dipolecorrection)) then
      !  allocate(vdplmt(lmmaxvr,nrmax,natmtot),vdplir(ngrtot))
      !  call dipole_correction(vdplmt,vdplir)
      !  do is = 1, nspecies
      !    do ia = 1, natoms(is)
      !      ias = idxas(ia,is)
      !      do ir = 1, nr(is)
      !        zvclmt(1,ir,ias) = zvclmt(1,ir,ias)+vdplmt(1,ir,ias)
      !        zvclmt(3,ir,ias) = zvclmt(3,ir,ias)+vdplmt(3,ir,ias)
      !      end do
      !    end do
      !  end do
      !end if

! compute (R_mt)^l
      Do is = 1, nspecies
         rmtl (0, is) = 1.d0
         Do l = 1, input%groundstate%lmaxvr + 3
            rmtl (l, is) = rmtl (l-1, is) * rmt (is)
         End Do
      End Do
! compute the multipole moments from the muffin-tin potentials
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            lm = 0
            Do l = 0, input%groundstate%lmaxvr
               t1 = dble (2*l+1) * rmtl (l+1, is) / fourpi
               Do m = - l, l
                  lm = lm + 1
                  qmt (lm, ias) = t1 * zvclmt (lm, nr(is), ias)
               End Do
            End Do
!            write(*,*) qmt (:, ias)
!            write(*,*)
         End Do
      End Do

! Fourier transform density to G-space and store in zvclir
      zvclir (:) = zrhoir (:)
      Call zfftifc (3, ngrid,-1, zvclir)
! find the multipole moments of the interstitial charge density
      qi (:, :) = 0.d0
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(qilocal,ig,ifg,zt1,t1,lm,t2,zt2,m,l,nthreads,whichthread,ithr) SHARED(input,gpc,sfacgp,zvclir,rmt,qi,ias,ngvec,is,ylmgp,jlgpr,zil,rmtl,igfft)
            qilocal=0d0
!$OMP DO    
#else
            qilocal=0d0
#endif
            Do ig = 1, ngvec
               ifg = igfft (ig)
               If (gpc(ig) .Gt. input%structure%epslat) Then
                  zt1 = fourpi * zvclir (ifg) * sfacgp (ig, ias)
                  t1 = 1.d0 / (gpc(ig)*rmt(is))
                  lm = 0
                  Do l = 0, input%groundstate%lmaxvr
                     t2 = t1 * rmtl (l+3, is) * jlgpr (l+1, ig, is)
                     zt2 = t2 * zt1 * zil (l)
                     Do m = - l, l
                        lm = lm + 1
                        qilocal (lm) = qilocal (lm) + zt2 * conjg (ylmgp(lm, ig))
                     End Do
                  End Do
               Else
                  t1 = fourpi * y00 * rmtl (3, is) / 3.d0
                  qilocal (1) = qilocal (1) + t1 * zvclir (ifg)
               End If
            End Do
#ifdef USEOMP
!$OMP END DO
            nthreads=omp_get_num_threads()
            whichthread=omp_get_thread_num()
            do ithr=0,nthreads-1
              if (ithr.eq.whichthread) then
                qi(:,ias)=qi(:,ias)+qilocal(:) 
              endif
!$OMP BARRIER
            enddo
!$OMP END PARALLEL
#else
            qi(:,ias)=qi(:,ias)+qilocal(:)
#endif
         End Do
      End Do

! find the smooth pseudocharge within the muffin-tin whose multipoles are the
! difference between the real muffin-tin and interstitial multipoles
      Do is = 1, nspecies
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            lm = 0
            Do l = 0, input%groundstate%lmaxvr
! note the factor 2^N*N! is omitted because of reciprocal term in the
! form-factor
               t1 = factnm (2*(l+input%groundstate%npsden)+3, 2) / &
              & factnm (2*l+1, 2)
               Do m = - l, l
                  lm = lm + 1
                  zrp (lm) = (qmt(lm, ias)-qi(lm, ias)) * t1
               End Do
            End Do
!            write(*,*) zrp (:)
!            write(*,*)
! add the pseudocharge and real interstitial densities in G-space
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(ig,ifg,zt1,t1,t2,zsum,m,l,lm,t3) SHARED(input,gpc,zvclir,rmt,qi,ias,ngvec,is,ylmgp,jlgpr,zil,rmtl,igfft,zrp,fpo,sfacgp)
!$OMP DO    
#endif
            Do ig = 1, ngvec
               ifg = igfft (ig)
               If (gpc(ig) .Gt. input%structure%epslat) Then
                  zt1 = fpo * conjg (sfacgp(ig, ias))
                  t1 = gpc (ig) * rmt (is)
                  t2 = t1 ** (input%groundstate%npsden+1)
                  lm = 0
                  Do l = 0, input%groundstate%lmaxvr
                     lm = lm + 1
                     zsum = zrp (lm) * ylmgp (lm, ig)
                     Do m = - l + 1, l
                        lm = lm + 1
                        zsum = zsum + zrp (lm) * ylmgp (lm, ig)
                     End Do
                     t3 = jlgpr (input%groundstate%npsden+l+1, ig, is) / (t2*rmtl(l, is))
                     zvclir (ifg) = zvclir (ifg) + t3 * zt1 * zsum * conjg (zil(l))
                  End Do
               Else
                  t1 = fpo * y00 / factnm (2*input%groundstate%npsden+3, 2)
                  zvclir (ifg) = zvclir (ifg) + t1 * zrp (1)
               End If
            End Do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL 
#endif

         End Do
      End Do





!------------------------------

#ifdef PSOLVER
if (.not.excite) then
!write(*,*)"psolver"
! Fourier transform interstitial potential to real space
      Call zfftifc (3, ngrid, 1, zvclir)
      allocate(fake_arr(1),r_v(n1*n2*n3), c_v(n1*n2*n3))


      !if (psolver0d) then
          !write(*,*)"reorder"
      !    call reorder(bc_type, zvclir,ngrid, r_v, c_v)     
      !else 
          !write(*,*)"else order"
          r_v=dble(zvclir)
          c_v=dimag(zvclir)
      !end if

      offset=0d0


      call H_potential('G',kernel,r_v,fake_arr,energy,offset,.false.,quiet='no')

      call H_potential('G',kernel,c_v,fake_arr,energy,offset,.false.,quiet='yes')


      zvclir=dcmplx(r_v,c_v)!combine complex and real solutions


     !if (bc_type.eq.0) then
        !write(*,*)"deorder"
     !   call reorder(bc_type, zvclir,ngrid, r_v, c_v)
        !zvclir=dcmplx(r_v,c_v)
      !end if
      deallocate(fake_arr, r_v, c_v)


! Fourier transform interstitial potential to reciprocal space
      Call zfftifc (3, ngrid, -1, zvclir)
      !zvclir(1)=0d0
      



end if!psolver

#endif
if (excite) then!exciting or exciting0d
!write(*,*)"exciting0"
! set zrho0 (pseudocharge density coefficient of the smallest G+p vector)
      ifg = igfft (igp0)
      zrho0 = zvclir (ifg)
      !zvclir (ifg) = 0.d0
      


!r_c = (3d0*omega*nkptnr/(fourpi))**(1d0/3d0)
r_c = (omega*nkptnr)**(1d0/3d0)*0.50d0!0.171d0


! solve Poissons's equation in G-space for the pseudocharge
If (associated(input%groundstate%hybrid)) Then 

   If ((input%groundstate%vha.eq."exciting0d").or.(input%groundstate%hybrid%singularity.eq."exc0d")) Then

      Do ig = 1, ngvec
         ifg = igfft (ig)
         If (gpc(ig) .Gt. input%structure%epslat) Then 
            zvclir (ifg) = fourpi * zvclir (ifg)*(1d0-cos(gpc(ig) * r_c )) / (gpc(ig)**2)
         Else
            zvclir (ifg) = zvclir(ifg)*(fourpi*0.5d0)*r_c**2
         End If
      End Do

   Else
   
      ifg = igfft (igp0)
      zrho0 = zvclir (ifg)
      zvclir (ifg) = 0.d0

      Do ig = 1, ngvec
         ifg = igfft (ig)
         If (gpc(ig) .Gt. input%structure%epslat) Then
            zvclir (ifg) = fourpi * zvclir (ifg) / (gpc(ig)**2)
         Else
            zvclir (ifg) = 0.d0
         End If
      End Do

   End If

Else If (input%groundstate%vha.eq."exciting0d") Then

   Do ig = 1, ngvec
      ifg = igfft (ig)
      If (gpc(ig) .Gt. input%structure%epslat) Then
         zvclir (ifg) = fourpi * zvclir (ifg)*(1d0-cos(gpc(ig) * r_c )) / (gpc(ig)**2)
      Else
         zvclir (ifg) = zvclir(ifg)*(fourpi*0.5d0)*r_c**2
      End If
   End Do

Else

   ifg = igfft (igp0)
   zrho0 = zvclir (ifg)
   zvclir (ifg) = 0.d0
   Do ig = 1, ngvec
      ifg = igfft (ig)
      If (gpc(ig) .Gt. input%structure%epslat) Then
         zvclir (ifg) = fourpi * zvclir (ifg) / (gpc(ig)**2)
      Else
         zvclir (ifg) = 0.d0
      End If
   End Do

End If
End If ! excite check



!------------------------------

! match potentials at muffin-tin boundary by adding homogeneous solution
      Do is = 1, nspecies
! compute (r/R_mt)^l
         Do ir = 1, nr (is)
            t1 = r (ir, is) / rmt (is)
            rl (ir, 0) = 1.d0
            Do l = 1, input%groundstate%lmaxvr
               rl (ir, l) = rl (ir, l-1) * t1
            End Do
         End Do
         Do ia = 1, natoms (is)
            ias = idxas (ia, is)
! find the spherical harmonic expansion of the interstitial potential at the
! muffin-tin radius
            vilm (:) = 0.d0

#ifdef USEOMP
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(ig,ifg,zt1,zt2,zsum,m,l,lm,t3,qilocal,nthreads,whichthread,ithr) SHARED(input,zvclir,ias,ngvec,is,ylmgp,jlgpr,zil,rmtl,igfft,sfacgp,vilm)
            qilocal=0d0
!$OMP DO  
#else
            qilocal=0d0  
#endif
            Do ig = 1, ngvec
               ifg = igfft (ig)
               zt1 = fourpi * zvclir (ifg) * sfacgp (ig, ias)
               lm = 0
               Do l = 0, input%groundstate%lmaxvr
                  zt2 = jlgpr (l, ig, is) * zt1 * zil (l)
                  Do m = - l, l
                     lm = lm + 1
                     qilocal (lm) = qilocal (lm) + zt2 * conjg (ylmgp(lm, ig))
!                     vilm (lm) = vilm (lm) + zt2 * conjg (ylmgp(lm, ig))
                  End Do
               End Do
            End Do
#ifdef USEOMP
!$OMP END DO
            nthreads=omp_get_num_threads()
            whichthread=omp_get_thread_num()
            do ithr=0,nthreads-1
              if (ithr.eq.whichthread) then
                vilm  = vilm +qilocal
              endif
!$OMP BARRIER
            enddo
!$OMP END PARALLEL 
#else
            vilm  = vilm +qilocal 
#endif

! add homogenous solution

            zt1 = vilm (1) - zvclmt (1, nr(is), ias)
            vmad(ias)=vmad(ias)+zt1*rl(1,0)*y00

            lm = 0
            Do l = 0, input%groundstate%lmaxvr
               Do m = - l, l
                  lm = lm + 1
                  zt1 = vilm (lm) - zvclmt (lm, nr(is), ias)
                  Do ir = 1, nr (is)
                     zvclmt (lm, ir, ias) = zvclmt (lm, ir, ias) + zt1 &
                    & * rl (ir, l)
                  End Do
               End Do
            End Do
         End Do
      End Do
! Fourier transform interstitial potential to real-space
      Call zfftifc (3, ngrid, 1, zvclir)

! add dipole correction
      if ((iscl>0).and.(input%groundstate%dipolecorrection)) then
        allocate(vdplmt(lmmaxvr,nrmax,natmtot),vdplir(ngrtot))
        call dipole_correction(vdplmt,vdplir)
        do is = 1, nspecies
          do ia = 1, natoms(is)
            ias = idxas(ia,is)
            do ir = 1, nr(is)
              zvclmt(1,ir,ias) = zvclmt(1,ir,ias)+vdplmt(1,ir,ias)
              zvclmt(3,ir,ias) = zvclmt(3,ir,ias)+vdplmt(3,ir,ias)
            end do
          end do
        end do
        zvclir(:) =  zvclir(:)+vdplir(:)
        deallocate(vdplmt,vdplir)
      end if

#ifdef PSOLVER
if (.not.excite) then
   !call f_lib_finalize()
   call pkernel_free(kernel)
end if
#endif

      Return
End Subroutine
!EOC
