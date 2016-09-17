
!============================================================
! Based on the routine src_xs/angavsc0.f90 by S. Sagmeister
!============================================================

subroutine angular_averaging(iom,n,bi)

    use modinput
    use modmain
    use modgw
    use mod_mpi_gw, only : myrank
    use modxs,      only : symt2
    use invert
    
    implicit none
    ! input variables
    integer, intent(in) :: iom
    integer, intent(in) :: n
    complex(8), intent(InOut) :: bi(n,n)
    ! local varibles
    integer, parameter :: nsphcov = 5810, iq0 = 1
    integer :: i1, i2, j1, j2, ji
    integer :: lmmaxdielt
    integer :: iop, jop, itp, lm, ntpsph
    real(8) :: vomega, t00, r, qsz
    real(8) :: q0eps(3), modq0
    complex (8) :: L(3,3), dtns(3,3), w1, w2
    
    real(8), allocatable :: plat(:,:), p(:), tp(:,:), spc(:,:), w(:)
    complex(8), allocatable :: m00lm(:), mx0lm(:), mxxlm(:)
    complex(8), allocatable :: ei00(:), eix0(:), ei0x(:), eixx(:)
    complex(8), allocatable :: ei00lm(:), eix0lm(:), ei0xlm(:), eixxlm(:)
    complex(8), allocatable :: ylm(:), zylm(:, :)
    complex(8), allocatable :: u(:,:), v(:,:), s(:,:), t(:,:) 
    complex(8), allocatable :: e3(:,:), ie3(:,:)
      
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
    u(:,:) = epsw1(:,iom,:)
    ! row wing
    allocate(v(3,n))
    v(:,:) = transpose(epsw2(:,iom,:))
    
    ! definition of S{\alpha} (B.13)
    allocate(s(n,3))
    s(:,:) = matmul(bi,u)
    allocate(t(3,n))
    t(:,:) = matmul(v,bi)

    ! definition of L = H-conjg(U)*S (B.14) = \epsilon_{00}+LFE
    L(:,:) = epsh(iom,:,:)-matmul(v,s)
        
    ! symmetrize the dielectric tensor
    dtns(:,:) = L(:,:)
    do iop = 1, 3
      do jop = 1, 3
        call symt2app(iop,jop,1,symt2,dtns,L(iop,jop))
      end do
    end do

    !-------------------------------------------
    ! Store the symmetrized macroscopic tensor
    !-------------------------------------------
    eps00(iom,:,:) = L(:,:)

    !====================================================
    ! calculate the averaged inverse dielectric function
    !====================================================
    select case (trim(input%gw%scrcoul%sciavtype))
    
    case('isotropic')
    
      q0eps(:) = input%gw%scrcoul%q0eps(:)
      modq0 = q0eps(1)**2+q0eps(2)**2+q0eps(3)**2

      epsh(iom,1,1) = modq0 / (L(1,1)*q0eps(1)+L(2,2)*q0eps(2)+L(3,3)*q0eps(3))
      
      epsw1(:,iom,1) = -epsh(iom,1,1)* &
      &                (s(:,1)*q0eps(1)+s(:,2)*q0eps(2)+s(:,3)*q0eps(3))/dsqrt(modq0)

      epsw2(:,iom,1) = -epsh(iom,1,1)* &
      &                (t(1,:)*q0eps(1)+t(2,:)*q0eps(2)+t(3,:)*q0eps(3))/dsqrt(modq0)
      
      do j2 = 1, n
        do j1 = 1, n
          bi(j1,j2) = bi(j1,j2) + &
          &           epsw1(j1,iom,1)*epsw2(j2,iom,1)/epsh(iom,1,1)
        end do
      end do
      
    case ('anisotropic')
    
      ! scaling factor
      t00 = vomega/(twopi**3)
      
      ! number of points on a sphere
      ntpsph = input%gw%scrcoul%nleblaik
      
      ! spherical harmonics parameters
      lmmaxdielt = (input%gw%scrcoul%lmaxdielt+1)**2
      if (lmmaxdielt > ntpsph) Then
        write(*,*)
        write(*,'("Error(angular_averaging): lmmaxdielt > ntpsph: ", 2i6)') &
        &  lmmaxdielt, ntpsph
        write(*,*)
        stop
      end if
      
      allocate(plat(3,ntpsph), p(ntpsph))
      allocate(m00lm(lmmaxdielt), mx0lm(lmmaxdielt), mxxlm(lmmaxdielt))
      allocate(ei00(ntpsph), eix0(ntpsph), ei0x(ntpsph), eixx(ntpsph))
      allocate(ei00lm(lmmaxdielt), eix0lm(lmmaxdielt), &
      &        ei0xlm(lmmaxdielt), eixxlm(lmmaxdielt))
      allocate(ylm(lmmaxdielt), zylm(ntpsph,lmmaxdielt))
      
      ! generate Lebedev-Laikov grid
      allocate(spc(3,ntpsph))
      allocate(w(ntpsph))
      call leblaik(ntpsph, spc, w)
      
      ! generate theta and phi angles
      allocate(tp(2,ntpsph))
      do itp = 1, ntpsph
        call sphcrd(spc(:,itp), r, tp(:,itp))
      end do
      
      ! generate spherical harmonics on covering set
      do itp = 1, ntpsph
        call genylm(input%gw%scrcoul%lmaxdielt, tp(:, itp), ylm)
        zylm(itp, :) = ylm(:)
      end do
      
      ! unit vectors of spherical covering set in lattice coordinates
      plat = matmul(binv, spc)
      
      ! distances to subcell cell boundaries in reciprocal space
      do itp = 1, ntpsph
        p(itp) = 1.d0 / (2.d0*maxval(abs(input%gw%ngridq(:)*plat(:,itp)), 1))
      end do
      
      !========
      ! HEAD
      !========
      
      !------------
      ! 1/(p*L*p)
      !------------
      do itp = 1, ntpsph
        ei00(itp) = 1.d0 / dot_product(spc(:, itp), matmul(L,spc(:,itp)))
      end do
      
      ! calculate lm-expansion coefficients
      do lm = 1, lmmaxdielt
        ei00lm(lm) = fourpi * dot_product(zylm(:,lm), ei00*w)
        m00lm(lm) = fourpi * dot_product(zylm(:,lm), p*w)
        mx0lm(lm) = fourpi * dot_product(zylm(:,lm), p**2/2.d0*w)
        mxxlm(lm) = fourpi * dot_product(zylm(:,lm), p**3/3.d0*w)
      end do
      
      !--------------------------------
      ! subcell average (inverse head)
      !--------------------------------
      epsh(iom,1,1) = t00 * dot_product(mxxlm, ei00lm)
     
      !========
      ! WINGS
      !========
      
      do j1 = 1, n
      
        do itp = 1, ntpsph
          !----------------------------
          ! Inverse wing, -p*S/(p*L*p)
          !----------------------------
          eix0(itp) = -dot_product(spc(:,itp), s(j1,:)) * ei00(itp)
          !----------------------------
          ! Inverse wing, -p*T/(p*L*p)
          !----------------------------
          ei0x(itp) = -dot_product(spc(:,itp), t(:,j1)) * ei00(itp)
        end do ! itp
        do lm = 1, lmmaxdielt
          eix0lm (lm) = fourpi * dot_product(zylm(:,lm), eix0*w)
          ei0xlm (lm) = fourpi * dot_product(zylm(:,lm), ei0x*w)
        end do
        
        ! subcell average (wings)
        epsw1(j1,iom,1) = t00 * dot_product(mxxlm, eix0lm)
        epsw2(j1,iom,1) = t00 * dot_product(mxxlm, ei0xlm)
        
      end do ! j1
      
      !========
      ! BODY
      !========
      
      !--------------------------------------
      ! Inverse body, B^-1 + p*S p*T/(p*L*p)
      !--------------------------------------
      if (input%gw%scrcoul%sciavbd) then
        do j2 = 1, n
        do j1 = 1, n
          do itp = 1, ntpsph
            eixx(itp) = bi(j1,j2) + &
            &           dot_product(spc(:,itp), s(j1,:)) * &
            &           dot_product(spc(:,itp), t(:,j2)) / &
            &           ei00(itp)
          end do
          do lm = 1, lmmaxdielt
            eixxlm(lm) = fourpi * dot_product(zylm(:,lm), eixx*w)
          end do
          ! subcell average (body)
          bi(j1,j2) = t00 * dot_product(mxxlm, eixxlm)
        end do ! j1
        end do ! j2
      else
        ! In the original version bi is not changed!
        ! I decided to take into account the found averaged head and wings
        ! via block matrix inversion
        do j2 = 1, n
        do j1 = 1, n
          bi(j1,j2) = bi(j1,j2) + &
          &           epsw1(j1,iom,1) * &
          &           epsw2(j2,iom,1) / &
          &           epsh(iom,1,1)
        end do ! j1
        end do ! j2
      end if
      
      deallocate(ei00, eix0, ei0x, eixx, ei00lm, eix0lm, ei0xlm)
      deallocate(m00lm, mx0lm, mxxlm)
      deallocate(ylm, zylm, tp, spc, w, plat, p)
      
    case ('screendiag', 'invscreendiag')
    
      !
      ! NOT YET TESTED
      !
    
      if (input%gw%scrcoul%sciavbd) Then
        write(*,*)
        write(*, '("Error(angular_averaging): (inv)screendiag-method does not &
        &allow for averaging the body of W")')
        write(*,*)
        stop
      end if
      
      !--------------------------------------------------
      ! invert dielectric matrix including 3 times G=0 
      ! according to the limits q->0_x, q->0_y, q->0_z
      !--------------------------------------------------
      allocate(e3(n+3,n+3), ie3(n+3,n+3))
      
      ! G=0, G'=0 elements
      e3(1:3,1:3) = epsh(iom,:,:)
      
      ! G!=0, G'=0 components and vice versa
      do i1 = 1, 3
        do j2 = 1, n
          e3(i1,j2+3) = epsw1(j2,iom,i1)
        end do
      end do
      do j1 = 1, n
        do i2 = 1, 3
          e3(j1+3,i2) = epsw2(j1,iom,i2)
        end do
      end do
      do j1 = 1, n
        do j2 = 1, n
          e3 (j1+3, j2+3) = epsilon(j1,j2,iom)
        end do
      end do
            
      call zinvert_hermitian(0, e3, ie3)
         
      ! select again
      select case (trim(input%gw%scrcoul%sciavtype))
     
      case ('screendiag')
        ! head
        epsh(iom,1,1) = w2 / ((e3(1,1)+ie3(2,2)+ie3(3,3))/3.d0)
        ! wings, set to zero in this approximation
        epsw1(:,iom,1) = zzero
        epsw2(:,iom,1) = zzero
        ! body, only diagonal is assigned
        bi(:,:) = zzero
        forall (j1=1:n)
          bi(j1,j1) = 1.d0 / e3(j1+3,j1+3)
        end forall
               
      case ('invscreendiag')
        ! head
        epsh(iom,1,1) = w2 * (ie3(1,1)+ie3(2,2)+ie3(3,3))/3.d0
        ! wings
        forall (j1 = 1:n)
          epsw1(j1,iom,1) = w1 * (ie3(j1+3,1)+ie3(j1+3,2)+ie3(j1+1,3))/3.d0
          epsw2(j1,iom,1) = conjg(epsw1(j1,iom,1))
        end forall
        ! body
        forall (j1=1:n, j2=1:n)
          bi(j1,j2) = ie3(j1+3,j2+3)
        end forall
        
      end select
      deallocate (e3, ie3)
      
    case default
      write (*,*)
      write (*, '("Error(angular_averaging): invalid averaging method")')
      write (*,*)
      stop
      
    end select

    deallocate(u, v, s, t)
    
end subroutine
