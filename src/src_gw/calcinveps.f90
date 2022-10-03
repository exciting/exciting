!> Compute the inverse of the dielectric function, including
!> applying spherical averaging.
subroutine calcinveps(iomstart,iomend)
       use modmain
       use modgw
   
       implicit none
       integer(4), intent(in) :: iomstart, iomend
       integer(4) :: iom
       integer(4) :: im
       integer(4) :: info, lwork
       real(8)    :: tstart, tend
       complex(8), allocatable :: eps(:,:)
       complex(8), allocatable :: w2b(:), bw1(:)
       integer(4), allocatable :: ipiv(:)
       complex(8), allocatable :: work(:)
   
       character(len=10) :: sname="calcinveps"
   
       external zgetrf, zgetri
       external zhetrf, zhetri
       complex(8), external :: zdotc, zdotu
   
       call timesec(tstart)
   
       ! local arrays for body and wings
       allocate(eps(mbsiz,mbsiz))
       if (Gamma) allocate(bw1(mbsiz),w2b(mbsiz))
   
       ! LAPACK working arrays
       lwork = 64*mbsiz
       allocate(ipiv(mbsiz))
       allocate(work(lwork))
   
       ! lopp over frequencies
       do iom = iomstart, iomend
   
         ! array for body and its inverse
         eps(1:mbsiz,1:mbsiz) = epsilon(1:mbsiz,1:mbsiz,iom)
   
         select case (freq%fconv)
           case('refreq')
             call zgetrf(mbsiz,mbsiz,eps,mbsiz,ipiv,info)
             call errmsg0(info,sname,"calling zgetrf")
             call zgetri(mbsiz,eps,mbsiz,ipiv,work,lwork,info)
             call errmsg0(info,sname,"calling zgetri")

           ! TODO(Alex) Issue #132. Test replacing LU factorisation with Cholesky, for inversion
           ! It should be faster.  
           case('imfreq')
             call zgetrf(mbsiz,mbsiz,eps,mbsiz,ipiv,info)
             call errmsg0(info,sname,"calling zgetrf")
             call zgetri(mbsiz,eps,mbsiz,ipiv,work,lwork,info)
             call errmsg0(info,sname,"calling zgetri")
         end select
      
         if (Gamma) then
           call angular_averaging(iom, mbsiz, eps)
         endif 
   
         !--------------------------------------------
         ! store into the same array
         !--------------------------------------------
         epsilon(1:mbsiz,1:mbsiz,iom) = eps(1:mbsiz,1:mbsiz)
   
       enddo ! iom
   
       deallocate(ipiv, work)
       deallocate(eps)
       if (Gamma) deallocate(bw1, w2b)
   
       !===========================================
       ! \epsilon^{-1}_{ij}-\delta_{ij}
       !===========================================
       do iom = iomstart, iomend
         if (Gamma) epsh(1,1,iom) = epsh(1,1,iom)-zone
         do im = 1, mbsiz
           epsilon(im,im,iom) = epsilon(im,im,iom)-zone
         end do
       end do ! iom
   
       call timesec(tend)
       time_dfinv = time_dfinv+tend-tstart
   
       contains 
        ! TODO(Alex). Issue  141. Test restoring epsilon spherical averaging from exciting oxygen in GW
         subroutine angular_averaging(iom,n,bi)
           use modxs,      only : symt2
           implicit none
           ! input variables
           integer, intent(in) :: iom
           integer, intent(in) :: n
           complex(8), intent(InOut) :: bi(n,n)
   
           ! local variables
           integer, parameter :: nsphcov = 5810, iq0 = 1
           integer :: j1, j2, i
           integer :: iop, jop
           real(8) :: vomega, qsz
           real(8) :: q0eps(3), modq0
           complex (8) :: L(3,3), L_diag(3), dtns(3,3), w1, w2
           complex(8), allocatable :: u(:,:), v(:,:), s(:,:), t(:,:)
           !> Tolerance for zero
           real(8), parameter :: tol = 1.e-8
             
           ! crystal volume
           vomega = omega*product(input%gw%ngridq)
           ! Wigner-Seitz radius and spherical approximation to 1/q^2 average
           qsz = (6*pi**2/vomega)**(1.d0/3.d0)
           ! weight for sqrt(4pi)/q based on Wigner-Seitz radius
           w1 = qsz**2 * vomega / (4.d0*pi**2) * dsqrt(fourpi)
           ! weight for 4pi/q^2 based on Wigner-Seitz radius
           w2 = 2*qsz*vomega/pi
           
           ! definition of U_{\alpha} (column wing) (B.7)
           allocate(u(n,3))
           u(:,:) = epsw1(:,:,iom)
           ! row wing
           allocate(v(3,n))
           v(:,:) = transpose(epsw2(:,:,iom))
           
           ! definition of S{\alpha} (B.13)
           allocate(s(n,3))
           s(:,:) = matmul(bi,u)
           allocate(t(3,n))
           t(:,:) = matmul(v,bi)
   
           ! definition of L = H-conjg(U)*S (B.14) = \epsilon_{00}+LFE
           L(:,:) = epsh(:,:,iom) - matmul(v,s)
               
           ! symmetrize the dielectric tensor
           dtns(:,:) = L(:,:)
           do iop = 1, 3
             do jop = 1, 3
               call symt2app(iop, jop, 1, symt2, dtns, L(iop,jop))
             end do
           end do
   
           !-------------------------------------------
           ! Store the symmetrized macroscopic tensor
           !-------------------------------------------
           eps00(:,:,iom) = L(:,:)
           L_diag = [(L(i, i), i = 1, 3)]
   
           !====================================================
           ! calculate the averaged inverse dielectric function
           !====================================================
           q0eps = [1.0d0, 1.0d0, 1.0d0]
           modq0 = dot_product(q0eps, q0eps)

           ! TODO(Alex) Note, sometimes L is zero, making the denominator also zero
           ! hence initialise as 0 and only set for finite values
           epsh(1, 1, iom) = cmplx(0.d0, 0.d0)
           if (real(dot_product(L_diag, L_diag)) > tol) then
               epsh(1, 1, iom) = modq0 / dot_product(L_diag, q0eps)
           endif
               
           epsw1(:,1,iom) = -epsh(1,1,iom)* &
           &                (s(:,1)*q0eps(1)+s(:,2)*q0eps(2)+s(:,3)*q0eps(3))/dsqrt(modq0)
   
           epsw2(:,1,iom) = -epsh(1,1,iom)* &
           &                (t(1,:)*q0eps(1)+t(2,:)*q0eps(2)+t(3,:)*q0eps(3))/dsqrt(modq0)
           
           do j2 = 1, n
             do j1 = 1, n
               bi(j1,j2) = bi(j1,j2) + &
               &           epsw1(j1,1,iom)*epsw2(j2,1,iom)/epsh(1,1,iom)
             end do
           end do
             
           deallocate(u, v, s, t)
         end subroutine
   
   end subroutine
   