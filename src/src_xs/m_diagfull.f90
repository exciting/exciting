module m_diagfull
  use modmpi

  implicit none

  private

  public :: diagfull

  contains

    !BOP
    ! !ROUTINE: diagfull
    ! !INTERFACE:
    subroutine diagfull(n, ham, evalre, evalim, evecr, evecl, fbalance, frcond, fsort)
    ! !INPUT/OUTPUT PARAMETERS:
    ! IN:
    !   integer(4) :: n
    !   complex(8) :: ham(n,n)          ! Complex square matrix 
    !   logical, optional :: fbalance   ! Balance matrix
    !   logical, optional :: frcond     ! Calculate reciprocal conditioning numbers
    !   logical, optional :: fsort      ! Sort solutions by real part of eigenvalues
    ! OUT:
    !   real(8) :: evalre(n)                        ! Real part of eigenvalues
    !   real(8), optional :: evalim(n)              ! Imaginary part of eigenvalues
    !   complex(8), target, optional :: evecr(n, n) ! Right eigenvectors as columns
    !   complex(8), target, optional :: evecl(n, n) ! Left eigenvectors as columns
    !
    ! !DESCRIPTION:
    !   This is a wrapper routine around the lapack routine {\tt ZGEEVX} for the
    !   digitalization of a general complex double precision matrix.
    !   It returns real, and optionally imaginary, part of the eigenvalues and 
    !   optionally right and/or left eigenvectors. 
    !   There are flags for balancing the matrix and calculating the 
    !   reciprocal conditioning numbers of eigenvalues and eigenvectors.
    !   Note: If frcond = .ture. the left and right eigenvectors are 
    !         still computed internally.
    !
    ! !REVISION HISTORY:
    !   Created 2016 (Aurich)
    !
    !EOP
    !BOC

      implicit none

      ! I/O
      integer(4), intent(in) :: n
      complex(8), intent(in) :: ham(n, n)
      logical, intent(in), optional :: frcond, fbalance, fsort

      real(8), intent(out) :: evalre(n)
      real(8), intent(out), optional :: evalim(n)
      complex(8), intent(out), target, optional :: evecr(n, n)
      complex(8), intent(out), target, optional :: evecl(n, n)

      ! Local variables
      real(8) :: absnorm
      integer(4) :: lwork, info
      integer(4) :: ilo, ihi
      logical :: frc, fba, fso
      character(1) :: bchar, vrchar, vlchar, schar
      
      ! Arrays
      real(8) :: scalefac(n), rconde(n), rcondv(n)
      real(8) :: rwork(2*n)
      complex(8) :: evals(n)

      ! Allocatable arrays
      integer(4), allocatable :: idxsort(:)
      complex(8), allocatable :: work(:)
      complex(8), allocatable, target :: evecright(:,:)
      complex(8), allocatable, target :: evecleft(:,:)

      complex(8), pointer :: pevecr(:, :) => null()
      complex(8), pointer :: pevecl(:, :) => null()

      ! Check flags
      if(present(fbalance)) then
        fba = fbalance
      else
        fba = .false.
      end if

      if(present(frcond)) then
        frc = frcond
      else
        frc = .false.
      end if

      if(present(fsort)) then
        fso = fsort
      else
        fso = .false.
      end if

      ! Set chars
      if(fba) then
        bchar = 'B'
      else
        bchar = 'N'
      end if

      !! Check what is to be done
      ! Only eigenvalues
      if(.not. present(evecl) .and. .not. present(evecr)) then
        schar = 'N'
        vrchar = 'N'
        vlchar = 'N'
        ! For calculation of reciprocal conditioning numbers
        ! also the left eigenvectors are needed.
        if(frc) then
          schar = 'B'
          allocate(evecright(n,n))
          allocate(evecleft(n,n))
          vrchar = 'V'
          vlchar = 'V'
        end if
        pevecr => evecright
        pevecl => evecleft
      ! Only eigenvalues and right eigenvectors
      else if(.not. present(evecl)) then
        schar = 'N'
        vrchar = 'V'
        vlchar = 'N'
        ! For calculation of reciprocal conditioning numbers
        ! also the left eigenvectors are needed.
        if(frc) then
          schar = 'B'
          allocate(evecleft(n,n))
          vlchar = 'V'
        end if
        pevecr => evecr
        pevecl => evecleft
      ! Only eigenvalues and left eigenvectors
      else if(.not. present(evecr)) then
        schar = 'N'
        vrchar = 'N'
        vlchar = 'V'
        ! For calculation of reciprocal conditioning numbers
        ! also the right eigenvectors are needed.
        if(frc) then
          schar = 'B'
          allocate(evecright(n,n))
          vrchar = 'V'
        end if
        pevecr => evecright
        pevecl => evecl
      ! Eigenvalues and both right and left eigenvectors
      else
        schar = 'N'
        vrchar = 'V'
        vlchar = 'V'
        if(frc) then
          schar = 'B'
        end if
        pevecr => evecr
        pevecl => evecl
      end if

      ! Workspace query 
      allocate(work(1))
      call zgeevx(bchar, vlchar, vrchar, schar, n, ham, n,&
        & evals, pevecl, n, pevecr, n, ilo, ihi, scalefac,&
        & absnorm, rconde, rcondv, work, -1, rwork, info)
      lwork = int(work(1),4)
      deallocate(work)
      allocate(work(lwork))

      ! Diagonalize
      call zgeevx(bchar, vlchar, vrchar, schar, n, ham, n,&
        & evals, pevecl, n, pevecr, n, ilo, ihi, scalefac,&
        & absnorm, rconde, rcondv, work, lwork, rwork, info)
      deallocate(work)

      if(info .ne. 0) then
       write(*,*)
       write(*, '("Error(bsesoldiagfull): zgeevx returned non-zero info:", i6)') info
       write(*,*)
       call terminate
      end if

      ! Get real and imaginary parts
      evalre = dble(evals)
      if(present(evalim)) then
        evalim = aimag(evals)
      end if

      ! Sort according to real part
      if(fso) then 
        allocate(idxsort(n))
        call sortidx(n, evalre, idxsort)
        evalre = evalre(idxsort)
        if(present(evalim)) then
          evalim = evalim(idxsort)
        end if
        if(present(evecr)) then 
          pevecr(:,:) = pevecr(:, idxsort)
        end if
        if(present(evecl)) then 
          pevecl(:,:) = pevecl(:, idxsort)
        end if
      end if

      ! Deallocate eigenvector arrays that were only needed for
      ! the calculation of condition numbers.
      if(allocated(evecright)) deallocate(evecright)
      if(allocated(evecleft)) deallocate(evecleft)
      nullify(pevecr, pevecl)

      ! Write out coniditioning information
      if(frc) then
        call writecondition(n, absnorm, rconde, evalre, rcondv)
      end if
      
    end subroutine diagfull
    !EOC

    subroutine writecondition(n, abnrm, rconde, eigvalre, rcondv)
      use m_getunit

      integer(4), intent(in) :: n
      real(8), intent(in) :: abnrm, rconde(n), rcondv(n), eigvalre(n)

      integer(4) :: un, i
      real(8) :: eps, eerrbd(n), verrbd(n)

      ! External functions
      real(8), external :: dlamch

      ! Get machine precision
      eps = dlamch('E')

      ! Get simple error bounds
      do i = 1, n
        eerrbd(i) = eps*abnrm/rconde(i)
        verrbd(i) = eps*abnrm/rcondv(i)
      end do

      call getunit(un)

      open(un, file='COND.OUT', status='replace', action='write', form='formatted')
      write(un, '("One-norm of the matrix: ",E23.16)') abnrm
      write(un, '("Errorbound of eignevalues&
        & (modulus of difference between real an numeric)")') 
      write(un, '("eval, err, err/eval")') 
      do i=1,n
        write(un, '(E23.16,1x,E23.16,1x,E23.16)') eigvalre(i), eerrbd(i),&
          & eerrbd(i)/eigvalre(i)
      end do
      write(un, '("Errorbound of eigenvectors&
        & (angle between real and numeric solution")') 
      write(un, '("errorangle")') 
      do i=1,n
        write(un, '(E23.16)') verrbd(i)
      end do
      close(un)

    end subroutine writecondition

end module m_diagfull
