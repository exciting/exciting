! Parts of this code are copied/adopted from Mitsuaki Kawamura's
! library 'libtetrabz'. See below copyright notice.
!
! Copyright (c) 2014 Mitsuaki Kawamura
!
! Permission is hereby granted, free of charge, to any person obtaining a
! copy of this software and associated documentation files (the
! "Software"), to deal in the Software without restriction, including
! without limitation the rights to use, copy, modify, merge, publish,
! distribute, sublicense, and/or sell copies of the Software, and to
! permit persons to whom the Software is furnished to do so, subject to
! the following conditions:
! 
! The above copyright notice and this permission notice shall be included
! in all copies or substantial portions of the Software.
  
module mod_opt_tetra
  use mod_kpointset
  use modmpi
  implicit none
  private
  
  type t_set
      integer :: ntetra
      integer :: ttype
      logical :: isreduced
      real(8) :: wlsm( 4, 20)
      integer, allocatable :: tetra(:,:)
      real(8), allocatable :: tetwgt(:)
      logical :: initialized = .false.
  end type t_set

  real(8), save :: eps = 1.d-10

  public :: t_set, &
            opt_tetra_init, opt_tetra_destroy, &
            opt_tetra_wgt_theta,     opt_tetra_int_theta, &
            opt_tetra_wgt_delta,     opt_tetra_int_delta, &
            opt_tetra_wgt_deltadiff, opt_tetra_int_deltadiff, &
            opt_tetra_wgt_dbldelta,  opt_tetra_int_dbldelta, &
            opt_tetra_efermi
  
  contains

  ! NOTE:
  ! When an integral should be evaluated for many energies/frequencies
  ! the weight array easily explodes for dense integration grids.
  ! Therefore, for each integral, two routines are available:
  ! opt_tetra_wgt_xxx and opt_tetra_int_xxx
  ! The _wgt_ routines provide the integration weights for each k-point for 
  ! a single energie/frequency.
  ! The _int_ routines directly computes the integral (sum over k) for multiple
  ! energies/frequencies. The matrix elements need to be provided.

!##############################################################
!# INTEGRATION ROUTINES
!##############################################################

!**************************************************************
!* Heaviside theta integral
!**************************************************************
  subroutine opt_tetra_wgt_theta( self, nk, nb, eb, ne, e, wgt)
    type( t_set), intent( in)     :: self             ! tetrahedron set
    integer, intent( in)          :: nk               ! number of k-points
    integer, intent( in)          :: nb               ! number of bands
    real(8), intent( in)          :: eb( nb, nk)      ! band energies
    integer, intent( in)          :: ne               ! number of energies/frequencies
    real(8), intent( in)          :: e(*)             ! energy/frequency to evaluate at
    real(8), intent( out)         :: wgt( nb, nk, *)  ! integration weights

    integer :: nn, i, it, ib, ie, idx(4)
    real(8) :: et(4), wt(4), add(20)

    nn = 20
    if( self%ttype .eq. 1) nn = 4
    wgt(:,:,1:ne) = 0.d0
#ifdef USEOMP
!$omp parallel default( shared) private( i, ib, it, ie, et, idx, wt, add)
!$omp do
#endif
    do ib = 1, nb
      do it = 1, self%ntetra
        if( self%ttype .eq. 1) then
          et = eb( ib, self%tetra(1:4,it))
        else
          et = matmul( self%wlsm, eb( ib, self%tetra(:,it)))
        end if
        call opt_tetra_sort4( et, idx)
        et = et( idx)
        do ie = 1, ne
          call opt_tetra_getwgt_theta( et-e( ie), wt)
          wt( idx) = self%tetwgt( it)*wt
          if( maxval( abs( wt)) .gt. eps) then
            if( self%ttype .eq. 1) then
              add(1:4) = wt
            else
              add = matmul( wt, self%wlsm)
            end if
            do i = 1, nn
              wgt( ib, self%tetra( i, it), ie) = wgt( ib, self%tetra( i, it), ie) + add(i)
            end do
          end if
        end do
      end do
    end do
#ifdef USEOMP
!$omp end do
!$omp end parallel
#endif
    return
  end subroutine opt_tetra_wgt_theta

  subroutine opt_tetra_int_theta( self, nk, nb, eb, ne, e, ld, resr, matr, resc, matc)
    type( t_set), intent( in)          :: self              ! tetrahedron set
    integer, intent( in)               :: nk                ! number of k-points
    integer, intent( in)               :: nb                ! number of bands
    real(8), intent( in)               :: eb( nb, nk)       ! band energies
    integer, intent( in)               :: ne                ! number of energies/frequencies
    real(8), intent( in)               :: e( ne)            ! energies/frequencies to evaluate at
    integer, intent( in)               :: ld                ! leading dimension of matrix elements
    real(8), optional, intent( out)    :: resr( ne, ld, nb) ! the resulting integral (real)
    real(8), optional, intent( in)     :: matr( ld, nb, nk) ! the matrix elements (real)
    complex(8), optional, intent( out) :: resc( ne, ld, nb) ! the resulting integral (complex)
    complex(8), optional, intent( in)  :: matc( ld, nb, nk) ! the matrix elements (complex)

    integer :: nn, it, ib, ie, idx(4)
    real(8) :: et(nb,4), eti(4), wti(4)
    real(8) :: mtir(ld,4)
    complex(8) :: mtic(ld,4)
    logical :: r
    real(8), allocatable    :: mtr(:,:,:), mttr(:,:,:)
    complex(8), allocatable :: mtc(:,:,:), mttc(:,:,:)

    r = .false.
    if( present( resr)) r = .true.
    if( (.not. r) .and. (.not. present( resc))) then
      write(*,*)
      write(*,'("Error (opt_tetra_int_theta): Either a real or a complex result array must be provided.")')
      stop
    end if

    nn = 20
    if( self%ttype .eq. 1) nn = 4
    if( r) then
      resr = 0.d0
    else
      resc = cmplx( 0.d0, 0.d0, 8)
    end if
#ifdef USEOMP
!$omp parallel default( shared) private( ib, it, ie, et, mtr, mttr, mtc, mttc, mtir, mtic, idx, eti, wti)
#endif
    if( present( matr)) then
      allocate( mtr( ld, nb, 4))
      if( self%ttype .eq. 2) allocate( mttr( ld, nb, nn))
      mtr = 1.d0
    else if( present( matc)) then
      allocate( mtc( ld, nb, 4))
      if( self%ttype .eq. 2) allocate( mttc( ld, nb, nn))
      mtc = cmplx( 1.d0, 0.d0, 8)
    end if
#ifdef USEOMP
!$omp do
#endif
    do it = 1, self%ntetra
      if( self%ttype .eq. 1) then
        et = eb( :, self%tetra(1:4,it))
        if(         r .and. present( matr)) mtr = matr(:,:,self%tetra(1:4,it))
        if( (.not. r) .and. present( matc)) mtc = matc(:,:,self%tetra(1:4,it))
      else
        et = matmul( eb( :, self%tetra(:,it)), transpose( self%wlsm))
        if(         r .and. present( matr)) then
          mttr = matr(:,:,self%tetra(:,it))
          call dgemm( 'n', 't', ld*nb, 4, nn, 1.d0, &
                 mttr, ld*nb, &
                 self%wlsm, 4, 0.d0, &
                 mtr, ld*nb)
        end if
        if( (.not. r) .and. present( matc)) then
          mttc = matc(:,:,self%tetra(:,it))
          call zgemm( 'n', 't', ld*nb, 4, nn, cmplx( 1.d0, 0.d0, 8), &
                 mttc, ld*nb, &
                 cmplx( self%wlsm, 0.d0, 8), 4, cmplx( 0.d0, 0.d0, 8), &
                 mtc, ld*nb)
        end if
      end if
      do ib = 1, nb
        eti = et( ib, :)
        call opt_tetra_sort4( eti, idx)
        eti = eti( idx)
        if( r) then
          mtir = mtr( :, ib, idx)
        else
          mtic = mtc( :, ib, idx)
        end if
        do ie = 1, ne
          call opt_tetra_getwgt_theta( eti-e( ie), wti)
          if( maxval( abs( wti)) .gt. eps) then
            if( r) then
#ifdef USEOMP
!$omp critical
#endif
              call dgemv( 'n', ld, 4, self%tetwgt( it), &
                     mtir, ld, &
                     wti, 1, 1.d0, &
                     resr( ie, :, ib), 1)
#ifdef USEOMP
!$omp end critical
#endif
            else
#ifdef USEOMP
!$omp critical
#endif
              call zgemv( 'n', ld, 4, cmplx( self%tetwgt( it), 0.d0, 8), &
                     mtic, ld, &
                     cmplx( wti, 0.d0, 8), 1, cmplx( 1.d0, 0.d0, 8), &
                     resc( ie, :, ib), 1)
#ifdef USEOMP
!$omp end critical
#endif
            end if
          end if
        end do
      end do
    end do
#ifdef USEOMP
!$omp end do
#endif
    if( present( matr)) then
      deallocate( mtr)
      if( self%ttype .eq. 2) deallocate( mttr)
    else if( present( matc)) then
      deallocate( mtc)
      if( self%ttype .eq. 2) deallocate( mttc)
    end if
#ifdef USEOMP
!$omp end parallel
#endif

    return
  end subroutine opt_tetra_int_theta

  subroutine opt_tetra_getwgt_theta( e, wt)
    real(8), intent( in)  :: e(4)
    real(8), intent( out) :: wt(4)
    
    real(8) :: c, w(4,4), a(4,4)

    wt = 0.d0 
    if( 0.d0 .le. e(1)) return
    if( 0.d0 .ge. e(4)) then
      wt = 0.25d0
      return
    end if
    a = coeff(4,e)
    if( (e(1) .lt. 0.d0) .and. (0.d0 .lt. e(2))) then
      call tetra_a1( e, a, c, w)
      wt = c*0.25d0*sum( w, 1)
      return
    else if( (e(2) .le. 0.d0) .and. (0.d0 .lt. e(3))) then
      call tetra_b1( e, a, c, w)
      wt = c*sum( w, 1)
      call tetra_b2( e, a, c, w)
      wt = wt + c*sum( w, 1)
      call tetra_b3( e, a, c, w)
      wt = wt + c*sum( w, 1)
      wt = 0.25d0*wt
      return
    else if( (e(3) .le. 0.d0) .and. (0.d0 .lt. e(4))) then
      call tetra_c1( e, a, c, w)
      wt = c*sum( w, 1)
      call tetra_c2( e, a, c, w)
      wt = wt + c*sum( w, 1)
      call tetra_c3( e, a, c, w)
      wt = wt + c*sum( w, 1)
      wt = 0.25d0*wt
      return
    end if
  end subroutine opt_tetra_getwgt_theta
  
!**************************************************************
!* Dirac delta integral
!**************************************************************
  subroutine opt_tetra_wgt_delta( self, nk, nb, eb, ne, e, wgt)
    type( t_set), intent( in)     :: self             ! tetrahedron set
    integer, intent( in)          :: nk               ! number of k-points
    integer, intent( in)          :: nb               ! number of bands
    real(8), intent( in)          :: eb( nb, nk)      ! band energies
    integer, intent( in)          :: ne               ! number of energies/frequencies
    real(8), intent( in)          :: e(*)             ! energy/frequency to evaluate at
    real(8), intent( out)         :: wgt( nb, nk, *)  ! integration weights

    integer :: nn, i, it, ib, ie, idx(4)
    real(8) :: et(4), wt(4), add(20)

    nn = 20
    if( self%ttype .eq. 1) nn = 4
    wgt(:,:,1:ne) = 0.d0
#ifdef USEOMP
!$omp parallel default( shared) private( i, ib, it, ie, et, idx, wt, add)
!$omp do
#endif
    do ib = 1, nb
      do it = 1, self%ntetra
        if( self%ttype .eq. 1) then
          et = eb( ib, self%tetra(1:4,it))
        else
          et = matmul( self%wlsm, eb( ib, self%tetra(:,it)))
        end if
        call opt_tetra_sort4( et, idx)
        et = et( idx)
        do ie = 1, ne
          call opt_tetra_getwgt_delta( et-e( ie), wt)
          wt( idx) = self%tetwgt( it)*wt
          if( maxval( abs( wt)) .gt. eps) then
            if( self%ttype .eq. 1) then
              add(1:4) = wt
            else
              add = matmul( wt, self%wlsm)
            end if
            do i = 1, nn
              wgt( ib, self%tetra( i, it), ie) = wgt( ib, self%tetra( i, it), ie) + add(i)
            end do
          end if
        end do
      end do
    end do
#ifdef USEOMP
!$omp end do
!$omp end parallel
#endif
    return
  end subroutine opt_tetra_wgt_delta

  subroutine opt_tetra_int_delta( self, nk, nb, eb, ne, e, ld, resr, matr, resc, matc)
    type( t_set), intent( in)          :: self              ! tetrahedron set
    integer, intent( in)               :: nk                ! number of k-points
    integer, intent( in)               :: nb                ! number of bands
    real(8), intent( in)               :: eb( nb, nk)       ! band energies
    integer, intent( in)               :: ne                ! number of energies/frequencies
    real(8), intent( in)               :: e( ne)            ! energies/frequencies to evaluate at
    integer, intent( in)               :: ld                ! leading dimension of matrix elements
    real(8), optional, intent( out)    :: resr( ne, ld, nb) ! the resulting integral (real)
    real(8), optional, intent( in)     :: matr( ld, nb, nk) ! the matrix elements (real)
    complex(8), optional, intent( out) :: resc( ne, ld, nb) ! the resulting integral (complex)
    complex(8), optional, intent( in)  :: matc( ld, nb, nk) ! the matrix elements (complex)

    integer :: nn, it, ib, ie, idx(4)
    real(8) :: et(nb,4), eti(4), wti(4)
    real(8) :: mtir(ld,4)
    complex(8) :: mtic(ld,4)
    logical :: r
    real(8), allocatable    :: mtr(:,:,:), mttr(:,:,:)
    complex(8), allocatable :: mtc(:,:,:), mttc(:,:,:)

    r = .false.
    if( present( resr)) r = .true.
    if( (.not. r) .and. (.not. present( resc))) then
      write(*,*)
      write(*,'("Error (opt_tetra_int_delta): Either a real or a complex result array must be provided.")')
      stop
    end if

    nn = 20
    if( self%ttype .eq. 1) nn = 4
    if( r) then
      resr = 0.d0
    else
      resc = cmplx( 0.d0, 0.d0, 8)
    end if
#ifdef USEOMP
!$omp parallel default( shared) private( it, ib, ie, et, mtr, mttr, mtc, mttc, mtir, mtic, idx, eti, wti)
#endif
    if( present( matr)) then
      allocate( mtr( ld, nb, 4))
      if( self%ttype .eq. 2) allocate( mttr( ld, nb, nn))
      mtr = 1.d0
    else if( present( matc)) then
      allocate( mtc( ld, nb, 4))
      if( self%ttype .eq. 2) allocate( mttc( ld, nb, nn))
      mtc = cmplx( 1.d0, 0.d0, 8)
    end if
#ifdef USEOMP
!$omp do
#endif
    do it = 1, self%ntetra
      if( self%ttype .eq. 1) then
        et = eb( :, self%tetra(1:4,it))
        if(         r .and. present( matr)) mtr = matr(:,:,self%tetra(1:4,it))
        if( (.not. r) .and. present( matc)) mtc = matc(:,:,self%tetra(1:4,it))
      else
        et = matmul( eb( :, self%tetra(:,it)), transpose( self%wlsm))
        if(         r .and. present( matr)) then
          mttr = matr(:,:,self%tetra(:,it))
          call dgemm( 'n', 't', ld*nb, 4, nn, 1.d0, &
                 mttr, ld*nb, &
                 self%wlsm, 4, 0.d0, &
                 mtr, ld*nb)
        end if
        if( (.not. r) .and. present( matc)) then
          mttc = matc(:,:,self%tetra(:,it))
          call zgemm( 'n', 't', ld*nb, 4, nn, cmplx( 1.d0, 0.d0, 8), &
                 mttc, ld*nb, &
                 cmplx( self%wlsm, 0.d0, 8), 4, cmplx( 0.d0, 0.d0, 8), &
                 mtc, ld*nb)
        end if
      end if
      do ib = 1, nb
        eti = et( ib, :)
        call opt_tetra_sort4( eti, idx)
        eti = eti( idx)
        if( r) then
          mtir = mtr( :, ib, idx)
        else
          mtic = mtc( :, ib, idx)
        end if
        do ie = 1, ne
          call opt_tetra_getwgt_delta( eti-e( ie), wti)
          if( maxval( abs( wti)) .gt. eps) then
            if( r) then
#ifdef USEOMP
!$omp critical
#endif
              call dgemv( 'n', ld, 4, self%tetwgt( it), &
                     mtir, ld, &
                     wti, 1, 1.d0, &
                     resr( ie, :, ib), 1)
#ifdef USEOMP
!$omp end critical
#endif
            else
#ifdef USEOMP
!$omp critical
#endif
              call zgemv( 'n', ld, 4, cmplx( self%tetwgt( it), 0.d0, 8), &
                     mtic, ld, &
                     cmplx( wti, 0.d0, 8), 1, cmplx( 1.d0, 0.d0, 8), &
                     resc( ie, :, ib), 1)
#ifdef USEOMP
!$omp end critical
#endif
            end if
          end if
        end do
      end do
    end do
#ifdef USEOMP
!$omp end do
#endif
    if( present( matr)) then
      deallocate( mtr)
      if( self%ttype .eq. 2) deallocate( mttr)
    else if( present( matc)) then
      deallocate( mtc)
      if( self%ttype .eq. 2) deallocate( mttc)
    end if
#ifdef USEOMP
!$omp end parallel
#endif

    return
  end subroutine opt_tetra_int_delta

  subroutine opt_tetra_getwgt_delta( e, wt)
    real(8), intent( in)  :: e(4)
    real(8), intent( out) :: wt(4)
    
    real(8) :: c, w(3,4), a(4,4)

    wt = 0.d0 
    if( (0.d0 .le. e(1)) .or. (0.d0 .ge. e(4))) return
    a = coeff(4,e)
    if( (e(1) .lt. 0.d0) .and. (0.d0 .lt. e(2))) then
      call trian_a1( e, a, c, w)
      wt = c*sum( w, 1)
      return
    else if( (e(2) .le. 0.d0) .and. (0.d0 .lt. e(3))) then
      call trian_b1( e, a, c, w)
      wt = c*sum( w, 1)
      call trian_b2( e, a, c, w)
      wt = wt + c*sum( w, 1)
      return
    else if( (e(3) .le. 0.d0) .and. (0.d0 .lt. e(4))) then
      call trian_c1( e, a, c, w)
      wt = c*sum( w, 1)
      return
    end if
  end subroutine opt_tetra_getwgt_delta

!**************************************************************
!* Dirac delta difference integral
!**************************************************************
  subroutine opt_tetra_wgt_deltadiff( self, nk, nb1, eb1, nb2, eb2, ne, e, wgt)
    type( t_set), intent( in) :: self                   ! tetrahedron set
    integer, intent( in)      :: nk                     ! number of k-points
    integer, intent( in)      :: nb1, nb2               ! number of bands
    real(8), intent( in)      :: eb1( nb1, nk)          ! first set of band energies
    real(8), intent( in)      :: eb2( nb2, nk)          ! second set of band energies
    integer, intent( in)      :: ne                     ! number of energies/frequencies
    real(8), intent( in)      :: e(*)                   ! energy/frequency to evaluate at
    real(8), intent( out)     :: wgt( nb1, nb2, nk, *)  ! integration weights

    integer :: nn, i, it, ib1, ib2, ie, idx(4)
    real(8) :: et1(nb1,4), et2(nb2,4), eti(4), wti(4), add(20)

    nn = 20
    if( self%ttype .eq. 1) nn = 4
    wgt(:,:,:,1:ne) = 0.d0
    do it = 1, self%ntetra
      if( self%ttype .eq. 1) then
        et1 = eb1( :, self%tetra(1:4,it))
        et2 = eb2( :, self%tetra(1:4,it))
      else
        et1 = matmul( eb1( :, self%tetra(:,it)), transpose( self%wlsm))
        et2 = matmul( eb2( :, self%tetra(:,it)), transpose( self%wlsm))
      end if
#ifdef USEOMP
!$omp parallel default( shared) private( i, ib1, ib2, ie, eti, idx, wti, add)
!$omp do collapse(2)
#endif
      do ib2 = 1, nb2
        do ib1 = 1, nb1
          eti = et1( ib1, :)-et2( ib2, :)
          call opt_tetra_sort4( eti, idx)
          eti = eti( idx)
          do ie = 1, ne
            call opt_tetra_getwgt_delta( eti-e( ie), wti)
            wti( idx) = self%tetwgt( it)*wti
            if( maxval( abs( wti)) .gt. eps) then
              if( self%ttype .eq. 1) then
                add(1:4) = wti
              else
                add = matmul( wti, self%wlsm)
              end if
              do i = 1, nn
                wgt( ib1, ib2, self%tetra( i, it), ie) = wgt( ib1, ib2, self%tetra( i, it), ie) + add(i)
              end do
            end if
          end do
        end do
      end do
#ifdef USEOMP
!$omp end do
!$omp end parallel
#endif
    end do
    return
  end subroutine opt_tetra_wgt_deltadiff

  subroutine opt_tetra_int_deltadiff( self, nk, nb1, eb1, nb2, eb2, ne, e, ld, resr, matr, resc, matc)
    type( t_set), intent( in)          :: self                         ! tetrahedron set
    integer, intent( in)               :: nk                           ! number of k-points
    integer, intent( in)               :: nb1, nb2                     ! number of bands
    real(8), intent( in)               :: eb1( nb1, nk)                ! first set of band energies
    real(8), intent( in)               :: eb2( nb2, nk)                ! second set of band energies
    integer, intent( in)               :: ne                           ! number of energies/frequencies
    real(8), intent( in)               :: e( ne)                       ! energies/frequencies to evaluate at
    integer, intent( in)               :: ld                           ! leading dimension of matrix elements
    real(8), optional, intent( out)    :: resr( ne, ld, nb1, nb2)      ! the resulting integral (real)
    real(8), optional, intent( in)     :: matr( ld, nb1, nb2, nk)      ! the matrix elements (real)
    complex(8), optional, intent( out) :: resc( ne, ld, nb1, nb2)      ! the resulting integral (complex)
    complex(8), optional, intent( in)  :: matc( ld, nb1, nb2, nk)      ! the matrix elements (complex)

    integer :: nn, it, ib1, ib2, ie, idx(4)
    real(8) :: et1(nb1,4), et2(nb2,4), eti(4), wti(4), mtdr(ld,4)
    complex(8) :: mtdc(ld,4)
    logical :: r
    real(8), allocatable    :: mtr(:,:,:,:), mttr(:,:,:,:)
    complex(8), allocatable :: mtc(:,:,:,:), mttc(:,:,:,:)

    r = .false.
    if( present( resr)) r = .true.
    if( (.not. r) .and. (.not. present( resc))) then
      write(*,*)
      write(*,'("Error (opt_tetra_int_deltadiff): Either a real or a complex result array must be provided.")')
      stop
    end if

    nn = 20
    if( self%ttype .eq. 1) nn = 4
    if( r) then
      resr = 0.d0
    else
      resc = cmplx( 0.d0, 0.d0, 8)
    end if
#ifdef USEOMP
!$omp parallel default( shared) private( it, ib1, ib2, ie, et1, et2, mtr, mttr, mtdr, mtc, mttc, mtdc, idx, eti, wti)
#endif
    if( present( matr)) then
      allocate( mtr( ld, nb1, nb2, 4))
      if( self%ttype .eq. 2) allocate( mttr( ld, nb1, nb2, nn))
      mtr = 1.d0
    else if( present( matc)) then
      allocate( mtc( ld, nb1, nb2, 4))
      if( self%ttype .eq. 2) allocate( mttc( ld, nb1, nb2, nn))
      mtc = cmplx( 1.d0, 0.d0, 8)
    end if
#ifdef USEOMP
!$omp do
#endif
    do it = 1, self%ntetra
      if( self%ttype .eq. 1) then
        et1 = eb1( :, self%tetra(1:4,it))
        et2 = eb2( :, self%tetra(1:4,it))
        if(         r .and. present( matr)) mtr = matr(:,:,:,self%tetra(1:4,it))
        if( (.not. r) .and. present( matc)) mtc = matc(:,:,:,self%tetra(1:4,it))
      else
        et1 = matmul( eb1( :, self%tetra(:,it)), transpose( self%wlsm))
        et2 = matmul( eb2( :, self%tetra(:,it)), transpose( self%wlsm))
        if(         r .and. present( matr)) then
          mttr = matr(:,:,:,self%tetra(:,it))
          call dgemm( 'n', 't', ld*nb1*nb2, 4, nn, 1.d0, &
                 mttr, ld*nb1*nb2, &
                 self%wlsm, 4, 0.d0, &
                 mtr, ld*nb1*nb2)
        end if
        if( (.not. r) .and. present( matc)) then
          mttc = matc(:,:,:,self%tetra(:,it))
          call zgemm( 'n', 't', ld*nb1*nb2, 4, nn, cmplx( 1.d0, 0.d0, 8), &
                 mttc, ld*nb1*nb2, &
                 cmplx( self%wlsm, 0.d0, 8), 4, cmplx( 0.d0, 0.d0, 8), &
                 mtc, ld*nb1*nb2)
        end if
      end if
      do ib2 = 1, nb2
        do ib1 = 1, nb1
          eti = et1( ib1, :)-et2( ib2, :)
          call opt_tetra_sort4( eti, idx)
          eti = eti( idx)
          if( r) then
            mtdr = mtr( :, ib1, ib2, idx)
          else
            mtdc = mtc( :, ib1, ib2, idx)
          end if
          do ie = 1, ne
            call opt_tetra_getwgt_delta( eti-e( ie), wti)
            if( maxval( abs( wti)) .gt. eps) then
              if( r) then
#ifdef USEOMP
!$omp critical
#endif
                call dgemv( 'n', ld, 4, self%tetwgt( it), &
                       mtdr, ld, &
                       wti, 1, 1.d0, &
                       resr( ie, :, ib1, ib2), 1)
#ifdef USEOMP
!$omp end critical
#endif
              else
#ifdef USEOMP
!$omp critical
#endif
                call zgemv( 'n', ld, 4, cmplx( self%tetwgt( it), 0.d0, 8), &
                       mtdc, ld, &
                       cmplx( wti, 0.d0, 8), 1, cmplx( 1.d0, 0.d0, 8), &
                       resc( ie, :, ib1, ib2), 1)
#ifdef USEOMP
!$omp end critical
#endif
              end if
            end if
          end do
        end do
      end do
    end do
#ifdef USEOMP
!$omp end do
#endif
    if( present( matr)) then
      deallocate( mtr)
      if( self%ttype .eq. 2) deallocate( mttr)
    else if( present( matc)) then
      deallocate( mtc)
      if( self%ttype .eq. 2) deallocate( mttc)
    end if
#ifdef USEOMP
!$omp end parallel
#endif
    return
  end subroutine opt_tetra_int_deltadiff

!**************************************************************
!* Double Dirac delta integral
!**************************************************************
  subroutine opt_tetra_wgt_dbldelta( self, nk, nb1, eb1, nb2, eb2, e1, e2, wgt)
    type( t_set), intent( in) :: self          ! tetrahedron set
    integer, intent( in)          :: nk                  ! number of k-points
    integer, intent( in)          :: nb1, nb2            ! number of bands
    real(8), intent( in)          :: eb1( nb1, nk)       ! first set of band energies
    real(8), intent( in)          :: eb2( nb2, nk)       ! first set of band energies
    real(8), intent( in)          :: e1, e2              ! energy/frequency to evaluate at
    real(8), intent( out)         :: wgt( nb1, nb2, nk)  ! integration weights

    integer :: nn, i, it, ib, idx(4)
    real(8) :: et1(nb1,4), et2(4), wt(nb1,4), add(nb1,20)

    nn = 20
    if( self%ttype .eq. 1) nn = 4
    wgt = 0.d0
#ifdef USEOMP
!$omp parallel default( shared) private( i, ib, it, et1, et2, idx, wt, add)
!$omp do
#endif
    do ib = 1, nb2
      do it = 1, self%ntetra
        if( self%ttype .eq. 1) then
          et1 = eb1( :, self%tetra(1:4,it))
          et2 = eb2( ib, self%tetra(1:4,it))
        else
          et1 = matmul( eb1( :, self%tetra(:,it)), transpose( self%wlsm))
          et2 = matmul( self%wlsm, eb2( ib, self%tetra(:,it)))
        end if
        call opt_tetra_sort4( et2, idx)
        et1 = et1(:,idx)
        et2 = et2( idx)
        call opt_tetra_getwgt_dbldelta( et2-e2, nb1, et1-e1, wt)
        wt(:,idx) = self%tetwgt( it)*wt
        if( maxval( abs( wt)) .gt. eps) then
          if( self%ttype .eq. 1) then
            add(:,1:4) = wt
          else
            add = matmul( wt, self%wlsm)
          end if
          do i = 1, nn
            wgt( :, ib, self%tetra( i, it)) = wgt( :, ib, self%tetra( i, it)) + add(:,i)
          end do
        end if
      end do
    end do
#ifdef USEOMP
!$omp end do
!$omp end parallel
#endif
    return
  end subroutine opt_tetra_wgt_dbldelta

  subroutine opt_tetra_int_dbldelta( self, nk, nb1, eb1, nb2, eb2, ne1, e1, ne2, e2, ld, res, mat)
    type( t_set), intent( in)         :: self                         ! tetrahedron set
    integer, intent( in)              :: nk                           ! number of k-points
    integer, intent( in)              :: nb1, nb2                     ! number of bands
    real(8), intent( in)              :: eb1( nb1, nk)                ! first set of band energies
    real(8), intent( in)              :: eb2( nb2, nk)                ! second set of band energies
    integer, intent( in)              :: ne1, ne2                     ! number of energies/frequencies
    real(8), intent( in)              :: e1( ne1)                     ! energies/frequencies to evaluate at
    real(8), intent( in)              :: e2( ne2)                     ! energies/frequencies to evaluate at
    integer, intent( in)              :: ld                           ! leading dimension of matrix elements
    complex(8), intent( out)          :: res( ne1, ne2, ld, nb1, nb2) ! the resulting integral
    complex(8), optional, intent( in) :: mat( ld, nb1, nb2, nk)       ! the matrix elements

    integer :: nn, it, ib, ie1, ie2, idx(4)
    real(8) :: et1(nb1,4), et2(nb2,4), ets1(nb1,4), ets2(4), wt(nb1,4)
    complex(8) :: mts(ld,nb1,4)
    complex(8), allocatable :: mt(:,:,:,:), mtt(:,:,:,:)

    nn = 20
    if( self%ttype .eq. 1) nn = 4
    res = cmplx( 0.d0, 0.d0, 8)
#ifdef USEOMP
!$omp parallel default( shared) private( ib, it, ie1, ie2, et1, et2, mt, mtt, idx, wt)
#endif
    allocate( mt( ld, nb1, nb2, 4))
    if( self%ttype .eq. 2) allocate( mtt( ld, nb1, nb2, nn))
    mt = cmplx( 1.d0, 0.d0, 8)
#ifdef USEOMP
!$omp do
#endif
    do it = 1, self%ntetra
      if( self%ttype .eq. 1) then
        et1 = eb1( :, self%tetra(1:4,it))
        et2 = eb2( :, self%tetra(1:4,it))
        if( present( mat)) mt = mat(:,:,:,self%tetra(1:4,it))
      else
        et1 = matmul( eb1( :, self%tetra(:,it)), transpose( self%wlsm))
        et2 = matmul( eb2( :, self%tetra(:,it)), transpose( self%wlsm))
        if( present( mat)) then
          mtt = mat(:,:,:,self%tetra(:,it))
          call zgemm( 'n', 't', ld*nb1*nb2, 4, nn, cmplx( 1.d0, 0.d0, 8), &
                 mtt, ld*nb1*nb2, &
                 cmplx( self%wlsm, 0.d0, 8), 4, cmplx( 0.d0, 0.d0, 8), &
                 mt, ld*nb1*nb2)
        end if
      end if
      do ib = 1, nb2
        call opt_tetra_sort4( et2(ib,:), idx)
        ets1 = et1( :, idx)
        ets2 = et2( ib, idx)
        mts = mt( :, :, ib, idx)
        do ie2 = 1, ne2
          do ie1 = 1, ne1
            call opt_tetra_getwgt_dbldelta( ets2-e2( ie2), nb1, ets1-e1( ie1), wt)
            if( maxval( abs( wt)) .gt. eps) then
              call zgemv( 'n', ld*nb1, 4, cmplx( self%tetwgt( it), 0.d0, 8), &
                     mts, ld*nb1, &
                     cmplx( wt, 0.d0, 8), 1, cmplx( 1.d0, 0.d0, 8), &
                     res( ie1, ie2, :, :, ib), 1)
            end if
          end do
        end do
      end do
    end do
#ifdef USEOMP
!$omp end do
#endif
    deallocate( mt)
    if( self%ttype .eq. 2) deallocate( mtt)
#ifdef USEOMP
!$omp end parallel
#endif
    return
  end subroutine opt_tetra_int_dbldelta

  subroutine opt_tetra_getwgt_dbldelta( e2, nb1, e1, wt)
    real(8), intent( in)  :: e2(4), e1( nb1, 4)
    integer, intent( in)  :: nb1
    real(8), intent( out) :: wt( nb1, 4)
    
    integer :: ib, idx(3)
    real(8) :: c, w(3,4), a(4,4), ec1( 3, nb1), wc( 3, nb1), wct(3)

    wt = 0.d0 
    if( (0.d0 .le. e2(1)) .or. (0.d0 .ge. e2(4))) return
    a = coeff(4,e2)
    if( (e2(1) .lt. 0.d0) .and. (0.d0 .lt. e2(2))) then
      call trian_a1( e2, a, c, w)
      ec1 = matmul( w, transpose( e1))
      do ib = 1, nb1
        call opt_tetra_sort3( ec1(:,ib), idx)
        ec1(:,ib) = ec1( idx, ib)
        call opt_tetra_getwgt_dbldelta2( ec1(:,ib), wct)
        wc(:,ib) = wct( idx)
      end do
      call dgemm( 't', 'n', nb1, 4, 3, c, &
             wc, 3, &
             w, 3, 0.d0, &
             wt, nb1)
      return
    else if( (e2(2) .le. 0.d0) .and. (0.d0 .lt. e2(3))) then
      call trian_b1( e2, a, c, w)
      ec1 = matmul( w, transpose( e1))
      do ib = 1, nb1
        call opt_tetra_sort3( ec1(:,ib), idx)
        ec1(:,ib) = ec1( idx, ib)
        call opt_tetra_getwgt_dbldelta2( ec1(:,ib), wct)
        wc(:,ib) = wct( idx)
      end do
      call dgemm( 't', 'n', nb1, 4, 3, c, &
             wc, 3, &
             w, 3, 0.d0, &
             wt, nb1)

      call trian_b2( e2, a, c, w)
      ec1 = matmul( w, transpose( e1))
      do ib = 1, nb1
        call opt_tetra_sort3( ec1(:,ib), idx)
        ec1(:,ib) = ec1( idx, ib)
        call opt_tetra_getwgt_dbldelta2( ec1(:,ib), wct)
        wc(:,ib) = wct( idx)
      end do
      call dgemm( 't', 'n', nb1, 4, 3, c, &
             wc, 3, &
             w, 3, 1.d0, &
             wt, nb1)
      return
    else if( (e2(3) .le. 0.d0) .and. (0.d0 .lt. e2(4))) then
      call trian_c1( e2, a, c, w)
      ec1 = matmul( w, transpose( e1))
      do ib = 1, nb1
        call opt_tetra_sort3( ec1(:,ib), idx)
        ec1(:,ib) = ec1( idx, ib)
        call opt_tetra_getwgt_dbldelta2( ec1(:,ib), wct)
        wc(:,ib) = wct( idx)
      end do
      call dgemm( 't', 'n', nb1, 4, 3, c, &
             wc, 3, &
             w, 3, 0.d0, &
             wt, nb1)
      return
    end if
  end subroutine opt_tetra_getwgt_dbldelta

  subroutine opt_tetra_getwgt_dbldelta2( e, wc)
    real(8), intent( in)  :: e(3)
    real(8), intent( out) :: wc(3)
    
    real(8) :: c, w(2,3), a(3,3)

    wc = 0.d0 
    if( (0.d0 .le. e(1)) .or. (0.d0 .ge. e(3))) return
    a = coeff(3,e)
    if( (e(1) .lt. 0.d0) .and. (0.d0 .lt. e(2))) then
      call line_a1( e, a, c, w)
      wc = c*sum( w, 1)
      return
    else if( (e(2) .le. 0.d0) .and. (0.d0 .lt. e(3))) then
      call line_b1( e, a, c, w)
      wc = c*sum( w, 1)
      return
    end if
  end subroutine opt_tetra_getwgt_dbldelta2

!##############################################################
!# ADVANCED TASKS
!##############################################################
  
!**************************************************************
!* Fermi energy
!**************************************************************
  subroutine opt_tetra_efermi( self, ne, nk, nb, eb, ef, occ, ef0, df0)
    type( t_set), intent( in)     :: self          ! tetrahedron set
    real(8), intent(in)           :: ne            ! number of electrons
    integer, intent(in)           :: nk            ! number of k-points in irr-BZ
    integer, intent(in)           :: nb            ! number of bands
    real(8), intent(in)           :: eb( nb, nk)   ! band energies
    real(8), intent(in), optional :: ef0, df0      ! initial guess and window
    real(8), intent(out)          :: ef            ! the fermi energy
    real(8), intent(out)          :: occ( nb, nk)  ! occupation numbers of each k
  
    integer :: iter, maxiter = 500
    real(8) :: elw, eup, sumk, his(2,2)
    complex(8) :: cocc( nb)
  
    ! find bounds for the fermi energy
    if( present( ef0) .and. present( df0)) then
      elw = ef0 - df0
      eup = ef0 + df0
    else
      elw = minval( eb)
      eup = maxval( eb)
    end if
  
    ! bisection method
    do iter = 1, maxiter
      
      if( (iter .gt. 2) .and. ( abs( his(2,2) - his(2,1)) .gt. 1d-8)) then
        ef = (ne - his(2,1))/(his(2,2) - his(2,1))*(his(1,2) - his(1,1)) + his(1,1)
        if( ef .lt. elw) elw = ef 
        if( ef .gt. eup) eup = ef
      else
        ef = (eup + elw)*0.5d0
      end if
      
      ! calc. # of electrons 
      call opt_tetra_wgt_theta( self, nk, nb, eb, 1, (/ef/), occ)
      !call opt_tetra_int_theta( self, nk, nb, eb, 1, (/ef/), 1, cocc)
      sumk = sum( dble( occ))
      !write(*,'(i,5f23.16)') iter, ef, elw, eup, sumk, ne
      his(1,1) = his(1,2)
      his(2,1) = his(2,2)
      his(1,2) = ef
      his(2,2) = sumk
      
      ! convergence check
      if( abs( sumk - ne) .lt. 1.d-8) then
        exit
      else if( sumk .lt. ne) then
        elw = ef
      else
        eup = ef
      end if
      
    end do
    if( iter .ge. maxiter) then
      write(*,*)
      write(*,*) "Warning (opt_tetra_efermi): Not converged", iter
    end if
  end subroutine opt_tetra_efermi

!##############################################################
!# HELPER ROUTINES
!##############################################################

!**************************************************************
!* Initialization
!**************************************************************
  subroutine opt_tetra_init( self, kset, tetra_type, reduce)
    use sorting, only: sort_index_1d, sort_index_2d
    !
    ! This routine set the corners and additional points for each tetrahedra
    !
    implicit none
    type( t_set), intent( out) :: self
    
    type( k_set), intent( in) :: kset
    integer, intent( in) :: tetra_type ! = 1 for Linear tetrahedron method
                                       ! = 2 for Optimized tetrahedron method
    logical, optional, intent( in) :: reduce
  
    ! local variables
    integer :: i1, i2, i3, itet, i, j, k, touch
    integer :: jk, nk1, nk2, nk3, ikv(3), nn, idx(4)
    integer :: ivvec( 3, 20, 6), divvec(4,4), ivvec0(4), tet(20)
    real(8) :: l(4), bvec2(3,3), bvec3(3,4)
    logical :: reduce_, exitloop
  
    integer, allocatable :: tetra(:,:), tetwgt(:), srt(:,:)
  
    reduce_ = .false.
    if( present( reduce)) reduce_ = reduce
    if( .not. kset%isreduced) reduce_ = .false.
  
    nk1 = kset%ngridk(1)
    nk2 = kset%ngridk(2)
    nk3 = kset%ngridk(3)
    
    ! Take the shortest diagonal line as the "shaft" of tetrahedral devision
    bvec2(:,1) = kset%bvec(:,1)/dble( nk1)
    bvec2(:,2) = kset%bvec(:,2)/dble( nk2)
    bvec2(:,3) = kset%bvec(:,3)/dble( nk3)
    
    bvec3(:,1) = -bvec2(:,1) + bvec2(:,2) + bvec2(:,3)
    bvec3(:,2) =  bvec2(:,1) - bvec2(:,2) + bvec2(:,3)
    bvec3(:,3) =  bvec2(:,1) + bvec2(:,2) - bvec2(:,3)
    bvec3(:,4) =  bvec2(:,1) + bvec2(:,2) + bvec2(:,3)
    
    do i = 1, 4
      l(i) = dot_product( bvec3( :, i), bvec3( :, i))
    end do
    
    i = minloc( l, 1)
    
    ivvec0 = (/ 0, 0, 0, 0 /)
    
    divvec(:,1) = (/ 1, 0, 0, 0 /)
    divvec(:,2) = (/ 0, 1, 0, 0 /)
    divvec(:,3) = (/ 0, 0, 1, 0 /)
    divvec(:,4) = (/ 0, 0, 0, 1 /)
    
    ivvec0(i) = 1
    divvec(i,i) = -1
    
    ! Divide a subcell into 6 tetrahedra
    itet = 0
    do i1 = 1, 3
      do i2 = 1, 3
        if(i2 == i1) cycle
        do i3 = 1, 3
          if(i3 == i1 .or. i3 == i2) cycle
          itet = itet + 1
          ivvec( :, 1, itet) = ivvec0( 1:3)
          ivvec( :, 2, itet) = ivvec( :, 1, itet) + divvec( 1:3, i1)
          ivvec( :, 3, itet) = ivvec( :, 2, itet) + divvec( 1:3, i2)
          ivvec( :, 4, itet) = ivvec( :, 3, itet) + divvec( 1:3, i3)
          !
        end do
      end do
    end do
    
    ! set weights for optimized tetrahedron method
    self%wlsm(1, 1: 4) = dble((/1440,    0,   30,    0/))
    self%wlsm(2, 1: 4) = dble((/   0, 1440,    0,   30/))
    self%wlsm(3, 1: 4) = dble((/  30,    0, 1440,    0/))
    self%wlsm(4, 1: 4) = dble((/   0,   30,    0, 1440/))
    !                                              
    self%wlsm(1, 5: 8) = dble((/ -38,    7,   17,  -28/))
    self%wlsm(2, 5: 8) = dble((/ -28,  -38,    7,   17/))
    self%wlsm(3, 5: 8) = dble((/  17,  -28,  -38,    7/))
    self%wlsm(4, 5: 8) = dble((/   7,   17,  -28,  -38/))
    !                                              
    self%wlsm(1, 9:12) = dble((/ -56,    9,  -46,    9/))
    self%wlsm(2, 9:12) = dble((/   9,  -56,    9,  -46/))
    self%wlsm(3, 9:12) = dble((/ -46,    9,  -56,    9/))
    self%wlsm(4, 9:12) = dble((/   9,  -46,    9,  -56/))
    !                                              
    self%wlsm(1,13:16) = dble((/ -38,  -28,   17,    7/))
    self%wlsm(2,13:16) = dble((/   7,  -38,  -28,   17/))
    self%wlsm(3,13:16) = dble((/  17,    7,  -38,  -28/))
    self%wlsm(4,13:16) = dble((/ -28,   17,    7,  -38/))
    !                                              
    self%wlsm(1,17:20) = dble((/ -18,  -18,   12,  -18/))
    self%wlsm(2,17:20) = dble((/ -18,  -18,  -18,   12/))
    self%wlsm(3,17:20) = dble((/  12,  -18,  -18,  -18/))
    self%wlsm(4,17:20) = dble((/ -18,   12,  -18,  -18/))
    !
    self%wlsm = self%wlsm/dble( 1260)
  
    if (tetra_type .eq. 1) then
      self%ttype = 1
      nn = 4
    else if(tetra_type .eq. 2) then
      self%ttype = 2
      nn = 20
    else
      self%ttype = 2
      write(*,*)
      write( *, '("Warning (opt_tetra_init): invalid tetra_type (",i3,")! I use type 2 (improved tetrahedron).")'), tetra_type
    end if
    
    ! construct tetrahedra
    self%ntetra = 6*nk1*nk2*nk3
    allocate( tetra( nn, self%ntetra), tetwgt( self%ntetra))
  
    self%ntetra = 0
    tetra = 0
    tetwgt = 1
  
    do i3 = 0, nk3-1
      do i2 = 0, nk2-1
        do i1 = 0, nk1-1
  
          do itet = 1, 6
  
            self%ntetra = self%ntetra + 1
            do i = 1, 4
              ikv = (/i1, i2, i3/) + ivvec( :, i, itet)
              ikv = modulo( ikv, kset%ngridk)
              tetra( i, self%ntetra) = kset%ikmap( ikv(1), ikv(2), ikv(3))
            end do
            if( self%ttype .eq. 2) then
              call opt_tetra_sort4( dble( tetra( 1:4, self%ntetra)), idx)
              tetra( 1:4, self%ntetra) = tetra( idx, self%ntetra)
              call opt_tetra_addCorners( ivvec( :, :, itet), order=idx)
              do i = 5, nn
                ikv = (/i1, i2, i3/) + ivvec( :, i, itet)
                ikv = modulo( ikv, kset%ngridk)
                tetra( i, self%ntetra) = kset%ikmap( ikv(1), ikv(2), ikv(3))
              end do
            end if
  
          end do
  
        end do
      end do
    end do
  
    if( reduce_) then
      allocate( srt( self%ntetra, 2))
      srt(:,1) = sort_index_2d( 4, self%ntetra, tetra, nn)
      srt(:,2) = sort_index_1d( self%ntetra, srt(:,1))
      ! Note: tetra = tetra(:,srt) sometimes does not work for dense grids
      do i = 1, self%ntetra
        j = srt(i,1)
        tet(1:nn) = tetra(:,j)
        tetra(:,j) = tetra(:,i)
        tetra(:,i) = tet(1:nn)
        srt(i,1) = srt(srt(i,2),1)
        srt(srt(i,2),1) = j
        srt(j,2) = srt(i,2)
      end do
      k = self%ntetra
      self%ntetra = 1
      outer: do i = 2, k
        do j = 1, 4
          if( tetra( j, i) .ne. tetra( j, self%ntetra)) then
            self%ntetra = self%ntetra + 1
            tetra( :, self%ntetra) = tetra( :, i)
            cycle outer
          end if
        end do
        tetwgt( self%ntetra) = tetwgt( self%ntetra) + 1
      end do outer
      deallocate( srt)
    end if
  
    if (allocated( self%tetra)) deallocate( self%tetra)
    allocate( self%tetra( nn, self%ntetra))
    if (allocated( self%tetwgt)) deallocate( self%tetwgt)
    allocate( self%tetwgt( self%ntetra))
  
    self%tetra = tetra( :, 1:self%ntetra)
    self%tetwgt = dble( tetwgt( 1:self%ntetra))/dble( 6*nk1*nk2*nk3)
    self%initialized = .true.
  
    deallocate( tetra, tetwgt)
    return
  
    contains
      subroutine opt_tetra_addCorners( ivvec, order)
        integer, intent( inout) :: ivvec( 3, 20)
        integer, optional, intent( in) :: order(4)
  
        integer :: o(4)
        
        o = (/1, 2, 3, 4/)
        if( present( order)) o = order
  
        ivvec( :,  5) = 2*ivvec( :, o(1)) - ivvec( :, o(2))
        ivvec( :,  6) = 2*ivvec( :, o(2)) - ivvec( :, o(3))
        ivvec( :,  7) = 2*ivvec( :, o(3)) - ivvec( :, o(4))
        ivvec( :,  8) = 2*ivvec( :, o(4)) - ivvec( :, o(1))
        
        ivvec( :,  9) = 2*ivvec( :, o(1)) - ivvec( :, o(3))
        ivvec( :, 10) = 2*ivvec( :, o(2)) - ivvec( :, o(4))
        ivvec( :, 11) = 2*ivvec( :, o(3)) - ivvec( :, o(1))
        ivvec( :, 12) = 2*ivvec( :, o(4)) - ivvec( :, o(2))
        
        ivvec( :, 13) = 2*ivvec( :, o(1)) - ivvec( :, o(4))
        ivvec( :, 14) = 2*ivvec( :, o(2)) - ivvec( :, o(1))
        ivvec( :, 15) = 2*ivvec( :, o(3)) - ivvec( :, o(2))
        ivvec( :, 16) = 2*ivvec( :, o(4)) - ivvec( :, o(3))
        
        ivvec( :, 17) =   ivvec( :, o(4)) - ivvec( :, o(1)) + ivvec( :, o(2))
        ivvec( :, 18) =   ivvec( :, o(1)) - ivvec( :, o(2)) + ivvec( :, o(3))
        ivvec( :, 19) =   ivvec( :, o(2)) - ivvec( :, o(3)) + ivvec( :, o(4))
        ivvec( :, 20) =   ivvec( :, o(3)) - ivvec( :, o(4)) + ivvec( :, o(1))
  
        return
      end subroutine opt_tetra_addCorners
  end subroutine opt_tetra_init

  subroutine opt_tetra_destroy( self)
    type( t_set), intent( inout) :: self

    if( allocated( self%tetra))  deallocate( self%tetra)
    if( allocated( self%tetwgt)) deallocate( self%tetwgt)
    self%initialized = .false.
  end subroutine opt_tetra_destroy

!**************************************************************
!* Miscellaneous
!**************************************************************

  function coeff(n,e) result(a)
    integer, intent( in) :: n
    real(8), intent( in) :: e(n)
    integer :: i, j
    real(8) :: de, a(n,n)
    a = 0.d0
    do j = 1, n
      do i = j+1, n
        de = e(i) - e(j)
        if( de .lt. eps) cycle
        a(i,j) = -e(j)/de
        a(j,i) = 1.d0-a(i,j)
      end do
    end do
    return
  end function coeff

  subroutine tetra_a1( e, a, c, w)
    real(8), intent( in)  :: e(4), a(4,4)
    real(8), intent( out) :: c, w(4,4)
    c = a(2,1)*a(3,1)*a(4,1)
    w(1,:) = (/  1.d0,   0.d0,   0.d0,   0.d0/)
    w(2,:) = (/a(1,2), a(2,1),   0.d0,   0.d0/)
    w(3,:) = (/a(1,3),   0.d0, a(3,1),   0.d0/)
    w(4,:) = (/a(1,4),   0.d0,   0.d0, a(4,1)/)
    return
  end subroutine tetra_a1

  subroutine tetra_b1( e, a, c, w)
    real(8), intent( in)  :: e(4), a(4,4)
    real(8), intent( out) :: c, w(4,4)
    c = a(3,1)*a(4,1)*a(2,4)
    w(1,:) = (/  1.d0,   0.d0,   0.d0,   0.d0/)
    w(2,:) = (/a(1,3),   0.d0, a(3,1),   0.d0/)
    w(3,:) = (/a(1,4),   0.d0,   0.d0, a(4,1)/)
    w(4,:) = (/  0.d0, a(2,4),   0.d0, a(4,2)/)
    return
  end subroutine tetra_b1

  subroutine tetra_b2( e, a, c, w)
    real(8), intent( in)  :: e(4), a(4,4)
    real(8), intent( out) :: c, w(4,4)
    c = a(3,2)*a(4,2)
    w(1,:) = (/  1.d0,   0.d0,   0.d0,   0.d0/)
    w(2,:) = (/  0.d0,   1.d0,   0.d0,   0.d0/)
    w(3,:) = (/  0.d0, a(2,3), a(3,2),   0.d0/)
    w(4,:) = (/  0.d0, a(2,4),   0.d0, a(4,2)/)
    return
  end subroutine tetra_b2

  subroutine tetra_b3( e, a, c, w)
    real(8), intent( in)  :: e(4), a(4,4)
    real(8), intent( out) :: c, w(4,4)
    c = a(2,3)*a(3,1)*a(4,2)
    w(1,:) = (/  1.d0,   0.d0,   0.d0,   0.d0/)
    w(2,:) = (/a(1,3),   0.d0, a(3,1),   0.d0/)
    w(3,:) = (/  0.d0, a(2,3), a(3,2),   0.d0/)
    w(4,:) = (/  0.d0, a(2,4),   0.d0, a(4,2)/)
    return
  end subroutine tetra_b3

  subroutine tetra_c1( e, a, c, w)
    real(8), intent( in)  :: e(4), a(4,4)
    real(8), intent( out) :: c, w(4,4)
    c = a(4,3)
    w(1,:) = (/  1.d0,   0.d0,   0.d0,   0.d0/)
    w(2,:) = (/  0.d0,   1.d0,   0.d0,   0.d0/)
    w(3,:) = (/  0.d0,   0.d0,   1.d0,   0.d0/)
    w(4,:) = (/  0.d0,   0.d0, a(3,4), a(4,3)/)
    return
  end subroutine tetra_c1

  subroutine tetra_c2( e, a, c, w)
    real(8), intent( in)  :: e(4), a(4,4)
    real(8), intent( out) :: c, w(4,4)
    c = a(3,4)*a(4,2)
    w(1,:) = (/  1.d0,   0.d0,   0.d0,   0.d0/)
    w(2,:) = (/  0.d0,   1.d0,   0.d0,   0.d0/)
    w(3,:) = (/  0.d0, a(2,4),   0.d0, a(4,2)/)
    w(4,:) = (/  0.d0,   0.d0, a(3,4), a(4,3)/)
    return
  end subroutine tetra_c2

  subroutine tetra_c3( e, a, c, w)
    real(8), intent( in)  :: e(4), a(4,4)
    real(8), intent( out) :: c, w(4,4)
    c = a(3,4)*a(2,4)*a(4,1)
    w(1,:) = (/  1.d0,   0.d0,   0.d0,   0.d0/)
    w(2,:) = (/a(1,4),   0.d0,   0.d0, a(4,1)/)
    w(3,:) = (/  0.d0, a(2,4),   0.d0, a(4,2)/)
    w(4,:) = (/  0.d0,   0.d0, a(3,4), a(4,3)/)
    return
  end subroutine tetra_c3

  subroutine trian_a1( e, a, c, w)
    real(8), intent( in)  :: e(4), a(4,4)
    real(8), intent( out) :: c, w(3,4)
    c = a(2,1)*a(3,1)/(e(4)-e(1))
    w(1,:) = (/a(1,2), a(2,1),   0.d0,   0.d0/)
    w(2,:) = (/a(1,3),   0.d0, a(3,1),   0.d0/)
    w(3,:) = (/a(1,4),   0.d0,   0.d0, a(4,1)/)
    return
  end subroutine trian_a1

  subroutine trian_b1( e, a, c, w)
    real(8), intent( in)  :: e(4), a(4,4)
    real(8), intent( out) :: c, w(3,4)
    c = a(4,1)*a(2,4)/(e(3)-e(1))
    w(1,:) = (/a(1,3),   0.d0, a(3,1),   0.d0/)
    w(2,:) = (/a(1,4),   0.d0,   0.d0, a(4,1)/)
    w(3,:) = (/  0.d0, a(2,4),   0.d0, a(4,2)/)
    return
  end subroutine trian_b1

  subroutine trian_b2( e, a, c, w)
    real(8), intent( in)  :: e(4), a(4,4)
    real(8), intent( out) :: c, w(3,4)
    c = a(2,3)*a(4,2)/(e(3)-e(1))
    w(1,:) = (/a(1,3),   0.d0, a(3,1),   0.d0/)
    w(2,:) = (/  0.d0, a(2,3), a(3,2),   0.d0/)
    w(3,:) = (/  0.d0, a(2,4),   0.d0, a(4,2)/)
    return
  end subroutine trian_b2

  subroutine trian_c1( e, a, c, w)
    real(8), intent( in)  :: e(4), a(4,4)
    real(8), intent( out) :: c, w(3,4)
    c = a(1,4)*a(2,4)/(e(4)-e(3))
    w(1,:) = (/a(1,4),   0.d0,   0.d0, a(4,1)/)
    w(2,:) = (/  0.d0, a(2,4),   0.d0, a(4,2)/)
    w(3,:) = (/  0.d0,   0.d0, a(3,4), a(4,3)/)
    return
  end subroutine trian_c1

  subroutine line_a1( e, a, c, w)
    real(8), intent( in)  :: e(3), a(3,3)
    real(8), intent( out) :: c, w(2,3)
    c = 0.5d0*a(2,1)/(e(3)-e(1))
    w(1,:) = (/a(1,2), a(2,1),   0.d0/)
    w(2,:) = (/a(1,3),   0.d0, a(3,1)/)
  end subroutine line_a1

  subroutine line_b1( e, a, c, w)
    real(8), intent( in)  :: e(3), a(3,3)
    real(8), intent( out) :: c, w(2,3)
    c = 0.5d0*a(2,3)/(e(3)-e(1))
    w(1,:) = (/a(1,3),   0.d0, a(3,1)/)
    w(2,:) = (/  0.d0, a(2,3), a(3,2)/)
  end subroutine line_b1

  subroutine opt_tetra_sort3( a, i)
    real(8), intent( in)  :: a(3)
    integer, intent( out) :: i(3)
  
    integer :: t

    i = (/1, 2, 3/)
    if( a(i(2)) .lt. a(i(1))) then
      i(1) = 2
      i(2) = 1
    end if
    if( a(3) .lt. a(i(2))) then
      i(3) = i(2)
      i(2) = 3
    end if
    if( a(i(2)) .lt. a(i(1))) then
      t = i(1)
      i(1) = i(2)
      i(2) = t
    end if
  end subroutine opt_tetra_sort3

  subroutine opt_tetra_sort4( a, i)
    real(8), intent( in)  :: a(4)
    integer, intent( out) :: i(4)
  
    integer :: t
  
    i = (/1, 2, 3, 4/)
    if( a(i(2)) .lt. a(i(1))) then
      i(1) = 2
      i(2) = 1
    end if
    if( a(i(4)) .lt. a(i(3))) then
      i(3) = 4
      i(4) = 3
    end if
    if( a(i(3)) .lt. a(i(1))) then
      t = i(1)
      i(1) = i(3)
      i(3) = t
    end if
    if( a(i(4)) .lt. a(i(2))) then
      t = i(2)
      i(2) = i(4)
      i(4) = t
    end if
    if( a(i(3)) .lt. a(i(2))) then
      t = i(2)
      i(2) = i(3)
      i(3) = t
    end if
  end subroutine opt_tetra_sort4

end module mod_opt_tetra

