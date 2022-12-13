
! Copyright (C) 2002-2010 J. K. Dewhurst, S. Sharma, C. Meisenbichler and
! C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!
!
!> spherical harmonic transform (SHT) matrices 
Module mod_SHT
      implicit none 
! real backward SHT matrix for lmaxapw
      Real (8), Allocatable :: rbshtapw (:, :)
! real forward SHT matrix for lmmaxapw
      Real (8), Allocatable :: rfshtapw (:, :)
! real backward SHT matrix for lmaxvr
      Real (8), Allocatable :: rbshtvr (:, :)
! real forward SHT matrix for lmaxvr
      Real (8), Allocatable :: rfshtvr (:, :)
! complex backward SHT matrix for lmaxapw
      Complex (8), Allocatable :: zbshtapw (:, :)
! complex forward SHT matrix for lmaxapw
      Complex (8), Allocatable :: zfshtapw (:, :)
! complex backward SHT matrix for lmaxvr
      Complex (8), Allocatable :: zbshtvr (:, :)
! complex forward SHT matrix for lmaxvr
      Complex (8), Allocatable :: zfshtvr (:, :)

  contains

    !> Generate forward and backward matrices for a real spherical harmonics
    !> transform.
    !> For a given \(l_{\rm max}\) a spherical covering with \(lm' = 1,\dots,(l_{\rm max}+1)^2\)
    !> pairs \((\theta_{lm'},\varphi_{lm'})\) is created. The backward transformation
    !> matrix is given by \(\texttt{bsht}_{lm',lm} = R_{lm}(\theta_{lm'},\varphi_{lm'})\).
    !> The forward matrix is given by the inverse of the backward matrix, i.e., 
    !> \(\sum_{lm''} \texttt{fsht}_{lm,lm''} \, \texttt{bsht}_{lm'',lm'} = \delta_{lm,lm'}\).
    subroutine gen_rshtmat( lmax, fsht, bsht)
      use precision, only: dp
      use modmpi, only: terminate_if_false
      !> maximum angular momentum \(l_{\rm max}\)
      integer, intent(in) :: lmax
      !> real forward SHT matrix
      real(dp), allocatable, intent(out) :: fsht(:,:)
      !> real backward SHT matrix
      real(dp), allocatable, intent(out) :: bsht(:,:)
    
      integer :: lmmax, itp, lwork, info
    
      integer, allocatable :: ipiv(:)
      real(dp), allocatable :: tp(:,:)
      real(dp), allocatable :: rlm(:), work(:)
    
      lmmax = (lmax + 1)**2
    
      if( allocated( fsht)) deallocate( fsht)
      if( allocated( bsht)) deallocate( bsht)
      allocate( fsht(lmmax,lmmax), bsht(lmmax,lmmax))
    
      lwork = 2*lmmax
      allocate( tp(2,lmmax), rlm(lmmax), ipiv(lmmax), work(lwork))
    
      ! generate spherical covering set for lmaxapw
      call sphcover( lmmax, tp)
      ! generate real spherical harmonics and set the backward SHT matrix
      do itp = 1, lmmax
        call genrlm( lmax, tp(:,itp), rlm)
        bsht(itp,:) = rlm
      end do
      ! find the forward SHT matrix
      fsht = bsht
      call dgetrf( lmmax, lmmax, fsht, lmmax, ipiv, info)
      if( info == 0) &
        call dgetri( lmmax, fsht, lmmax, ipiv, work, lwork, info)
      call terminate_if_false( info == 0, &
        message="(gen_rshtmat) Improper spherical covering. Unable to find &
                 inverse spherical harmonic transform.")
    
      deallocate( tp, rlm, ipiv, work)
    end subroutine gen_rshtmat
    
    !> Generate forward and backward matrices for a complex spherical harmonics
    !> transform.
    !> For a given \(l_{\rm max}\) a spherical covering with \(lm' = 1,\dots,(l_{\rm max}+1)^2\)
    !> pairs \((\theta_{lm'},\varphi_{lm'})\) is created. The backward transformation
    !> matrix is given by \(\texttt{bsht}_{lm',lm} = Y_{lm}(\theta_{lm'},\varphi_{lm'})\).
    !> The forward matrix is given by the inverse of the backward matrix, i.e., 
    !> \(\sum_{lm''} \texttt{fsht}_{lm,lm''} \, \texttt{bsht}_{lm'',lm'} = \delta_{lm,lm'}\).
    subroutine gen_zshtmat( lmax, fsht, bsht)
      use precision, only: dp
      use modmpi, only: terminate_if_false
      !> maximum angular momentum \(l_{\rm max}\)
      integer, intent(in) :: lmax
      !> real forward SHT matrix
      complex(dp), allocatable, intent(out) :: fsht(:,:)
      !> real backward SHT matrix
      complex(dp), allocatable, intent(out) :: bsht(:,:)
    
      integer :: lmmax, itp, lwork, info
    
      integer, allocatable :: ipiv(:)
      real(dp), allocatable :: tp(:,:)
      complex(dp), allocatable :: ylm(:), work(:)
    
      lmmax = (lmax + 1)**2
    
      if( allocated( fsht)) deallocate( fsht)
      if( allocated( bsht)) deallocate( bsht)
      allocate( fsht(lmmax,lmmax), bsht(lmmax,lmmax))
    
      lwork = 2*lmmax
      allocate( tp(2,lmmax), ylm(lmmax), ipiv(lmmax), work(lwork))
    
      ! generate spherical covering set for lmaxapw
      call sphcover( lmmax, tp)
      ! generate real spherical harmonics and set the backward SHT matrix
      do itp = 1, lmmax
        call genylm( lmax, tp(:,itp), ylm)
        bsht(itp,:) = ylm
      end do
      ! find the forward SHT matrix
      fsht = bsht
      call zgetrf( lmmax, lmmax, fsht, lmmax, ipiv, info)
      if( info == 0) &
        call zgetri( lmmax, fsht, lmmax, ipiv, work, lwork, info)
      call terminate_if_false( info == 0, &
        message="(gen_zshtmat) Improper spherical covering. Unable to find &
                 inverse spherical harmonic transform.")
    
      deallocate( tp, ylm, ipiv, work)
    end subroutine gen_zshtmat
End Module
!
