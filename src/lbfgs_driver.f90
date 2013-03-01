!                                                                                      
!  L-BFGS-B is released under the "New BSD License" (aka "Modified BSD License"
!  or "3-clause license")
!  Please read attached file License.txt
!                                        
!    --------------------------------------------------------------
!             CUSTOMIZED DRIVER FOR L-BFGS-B
!    --------------------------------------------------------------
!
!       L-BFGS-B is a code for solving large nonlinear optimization
!            problems with simple bounds on the variables.
!
!       The code can also be used for unconstrained problems and is
!       as efficient for these problems as the earlier limited memory
!                         code L-BFGS.
!
!       This driver illustrates how to control the termination of the
!       run and how to design customized output.
!
!    References:
!
!       [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
!       memory algorithm for bound constrained optimization'',
!       SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.
!
!       [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: FORTRAN
!       Subroutines for Large Scale Bound Constrained Optimization''
!       Tech. Report, NAM-11, EECS Department, Northwestern University,
!       1994.
!
!
!         (Postscript files of these papers are available via anonymous
!          ftp to eecs.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.)
!
!                             *  *  *
!
!         February 2011   (latest revision)
!         Optimization Center at Northwestern University
!         Instituto Tecnologico Autonomo de Mexico
!
!         Jorge Nocedal and Jose Luis Morales
!
!    **************

subroutine lbfgs_driver
!
!     Exciting interface created by DIN (February 2013)
!
      use modinput
      use modmain      
      Use modmpi

      implicit none
 
!     Declare variables and parameters needed by the code.
!
!     We suppress both code-supplied stopping tests because the
!     user is providing his/her own stopping criteria.
 
      integer,  parameter    :: dp = kind(1.0d0)
      real(dp), parameter    :: factr  = 0.0d0, pgtol  = 0.0d0
      
      integer                :: n
      integer                :: m, iprint
      integer                :: is, ia, ias, ik, ispn
      
      character(len=60)      :: ctask, csave
      logical                :: lsave(4)
      integer                :: isave(44)
      real(dp)               :: f
      real(dp)               :: dsave(29)
      integer,  allocatable  :: nbd(:), iwa(:)
      real(dp), allocatable  :: x(:), g(:)
      real(dp), allocatable  :: l(:), u(:), wa(:)
!
      real(dp)               :: v(3)
      integer                :: i, j
      character(1024)        :: message
      integer, allocatable   :: amap(:,:)
      logical                :: force_conv

      m = input%structureoptimization%lbfgsnumcor
      iprint = input%structureoptimization%lbfgsverbosity
      force_conv = .false.

!     Total number of variables taking into account constraints
      n = 0
      do is = 1, nspecies
        do ia = 1, natoms(is)
          do i = 1, 3
            if (.not.input%structure%speciesarray(is)%species%atomarray(ia)%atom%lock(i)) then
              n = n+1
            end if
          end do
        end do
      end do
      if (n==0) then
        call warning('WARNING(lbfgs_driver):')
        write(message,'(" No active degrees of freedom = Nothing to relax! Check lock options in your input file")')
        call warning(message)
        return
      end if
      
      allocate( nbd(n), x(n), l(n), u(n), g(n) )
      allocate( iwa(3*n) )
      allocate( wa(2*m*n + 5*n + 11*m*m + 8*m) )
      allocate( amap(3,n) )

      j = 0
      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas(ia,is)
          do i = 1, 3
            if (.not.input%structure%speciesarray(is)%species%atomarray(ia)%atom%lock(i)) then
              j = j+1
              amap(1,j)=i; amap(2,j)=ia; amap(3,j)=is
              nbd(j) = 2 ! constraint optimization (see src/Lbfgsb.3.0/README)
              x(j) = atposc(i,ia,is)
              ! l and u boundaries are used only when nbd > 0
              l(j) = x(j)-input%structureoptimization%tau0atm
              u(j) = x(j)+input%structureoptimization%tau0atm
            end if
          end do
        end do
      end do

!     BEGIN the loop
      ctask = 'START'
      do while (ctask(1:2).eq.'FG'    .or. &
                ctask(1:5).eq.'NEW_X' .or. &
                ctask(1:5).eq.'START')
      
!     This is the call to the L-BFGS-B code

        call setulb(n,m,x,l,u,nbd,f,g,factr,pgtol,wa,iwa,ctask,iprint, &
        &   csave,lsave,isave,dsave)

        if (ctask(1:2) .eq. 'FG') then

!         the minimization routine has returned to request the
!         function f and gradient g values at the current x.
          call calcEnergyForces
          
        else 
          
          if (ctask(1:5) .eq. 'NEW_X') then   
            
            call updatepositions
            
            If (rank .Eq. 0) Then
              Write (60,*)
              Write (60, '("+--------------------------+")')
              Write (60, '("| Updated atomic positions |")')
              Write (60, '("+--------------------------+")')

              do is = 1, nspecies
                write (60,*)
                write (60, '("Species : ", I4, " (", A, ")")') &
                &   is, trim (input%structure%speciesarray(is)%species%chemicalSymbol)
                write (60, '(" atomic positions (lattice) :")')
                do ia = 1, natoms(is)
                  write (60, '(I4, " : ", 3F14.8)') ia, &
                  &   input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:)
                end do ! ia
              end do ! is

              ! write lattice vectors and optimised atomic positions to file
              Call writehistory
              Call writegeometryxml (.True.)
              ! write the optimized interatomic distances to file
              Call writeiad (.True.)
            End If
            
            ! check force convergence
            If (forcemax .Le. input%structureoptimization%epsforce) Then
              If (rank .Eq. 0) Then
                Write (60,*)
                Write (60, '("Force convergence target achieved")')
              End If
              ctask = 'STOP'
              force_conv = .true.
            End If

          end if ! 'NEW_X'
          
        end if

      end do
!     END the loop

      if (.not.force_conv) then
        if (rank .Eq. 0) Then
          write(60,*)
          write(60,'("ATTENTION(L-BFGS): Required force convergence has not been reached!")')
          write(60,'(" forcemax=",f12.8," > epsforce=",f12.8)') forcemax, input%structureoptimization%epsforce
          write(60,*)
        end if
      end if
      
      deallocate( nbd, x, l, u, g )
      deallocate( iwa )
      deallocate( wa )

contains

    subroutine calcEnergyForces
     
        implicit none
        integer :: is, ia, i, j

! update atomic positions
        call updatepositions

!-----------------------!
!   Reinitialization
!-----------------------!
! check for overlapping muffin-tins
        call checkmt
! generate structure factors for G-vectors
        Call gensfacgp (ngvec, vgc, ngvec, sfacg)
! generate the characteristic function
        Call gencfun
! generate structure factors for G+k-vectors
        Do ik = 1, nkpt
            Do ispn = 1, nspnfv
                Call gensfacgp (ngk(ispn, ik), vgkc(:, :, ispn, ik), &
               &  ngkmax, sfacgk(:, :, ispn, ik))
            End Do
        End Do
! determine the new nuclear-nuclear energy
        Call energynn

!-----------------------!
!   SCF cycle
!-----------------------!
        call scf_cycle

!-----------------------!
!   Output
!-----------------------!
        ! total energy
        f = engytot
        
        ! energy gradients = total forces
        do j = 1, n
          i = amap(1,j); ia = amap(2,j); is = amap(3,j)
          ias = idxas(ia,is)
          g(j) = -forcetot(i,ias)
        end do
      
        return
    end subroutine calcEnergyForces
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
    subroutine updatepositions
        
        implicit none
      
        do j = 1, n
            i = amap(1,j); ia = amap(2,j); is = amap(3,j)
            atposc(i,ia,is) = x(j)
        end do
        do is = 1, nspecies
            do ia = 1, natoms(is)
! compute the lattice coordinates of the atomic positions
                call r3mv (ainv, atposc(:, ia, is), &
                &   input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:))
            end do
        end do

! find the crystal symmetries and shift atomic positions if required
        call findsymcrys
    
    end subroutine updatepositions
      
end subroutine
