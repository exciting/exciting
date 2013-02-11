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
 
      integer,  parameter    :: iprint = -1
      integer,  parameter    :: dp = kind(1.0d0)
      real(dp), parameter    :: factr  = 0.0d0, pgtol  = 0.0d0
      
      integer                :: n, m
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
      real(dp)               :: xc(3), xl(3), fc(3)
      real(dp)               :: t1, t2
      real(dp)               :: alpha
      integer                :: i, j
      Logical                :: force_converged

!    Total number of variables
      n = 3*natmtot
      m = 4 ! <-- this should go into input parameters
      
      alpha=10.d0

      allocate( nbd(n), x(n), l(n), u(n), g(n) )
      allocate( iwa(3*n) )
      allocate( wa(2*m*n + 5*n + 11*m*m + 8*m) )

      j = 0
      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas(ia,is)
          do i = 1, 3
!           apply the constrain
            if (input%structure%speciesarray(is)%species%atomarray(ia)%atom%lock(i)) cycle
            j = j+1
            nbd(j) = 0
            x(j) = atposc(i,ia,is)
            l(j) = atposc(i,ia,is)-input%structureoptimization%tau0atm
            u(j) = atposc(i,ia,is)+input%structureoptimization%tau0atm
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
       &  csave,lsave,isave,dsave)

        if (ctask(1:2) .eq. 'FG') then

!         the minimization routine has returned to request the
!         function f and gradient g values at the current x.
          call calcEnergyForces(n,x,f,g)
          
        else 
          
          if (ctask(1:5) .eq. 'NEW_X') then   
            
            ! write lattice vectors and optimised atomic positions to file
            Call writehistory
            Call writegeometryxml (.True.)
            ! write the optimised interatomic distances to file
            Call writeiad (.True.)

            Write (60,*)
            Write (60, '("+--------------------------+")')
            Write (60, '("| Updated atomic positions |")')
            Write (60, '("+--------------------------+")')
            
            j = 0
            do is = 1, nspecies
              Write (60,*)
              Write (60, '("Species : ", I4, " (", A, ")")') &
             &  is, trim (input%structure%speciesarray(is)%species%chemicalSymbol)
              Write (60, '(" atomic positions (lattice) :")')
              do ia = 1, natoms(is)
                ias = idxas(ia,is)
                do i = 1, 3
                  j = j+1
                  xc(i) = x(j)
                end do
                ! compute the lattice coordinates of the atomic positions
                Call r3mv(ainv, xc, xl)
                Write (60, '(I4, " : ", 3F14.8)') ia, xl
              end do
            end do
            
            ! check force convergence
            If (forcemax .Le. input%structureoptimization%epsforce) Then
               If (rank .Eq. 0) Then
                  Write (60,*)
                  Write (60, '("Force convergence target achieved")')
               End If
               ctask = 'STOP'
            End If
            
          end if ! 'NEW_X'
          
        end if

      end do
!     END the loop

      j = 0
      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas(ia,is)
          do i = 1, 3
            j = j+1
            atposc(i,ia,is) = x(j)
          end do
          ! compute the lattice coordinates of the atomic positions
          Call r3mv(ainv, atposc(:, ia, is), &
         &  input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:))
        end do
      end do

      deallocate( nbd, x, l, u, g )
      deallocate( iwa )
      deallocate( wa )

contains

    subroutine calcEnergyForces(ndim,x,f,g)
     
        implicit none
    
        integer :: ndim
        real(8) :: x(ndim), f, g(ndim)
        Integer :: is, ia, i, j
    
        j = 0
        do is = 1, nspecies
            do ia = 1, natoms(is)
                ias = idxas(ia,is)
                do i = 1, 3
                    j = j+1
                    atposc(i,ia,is) = x(j)
                end do
! compute the lattice coordinates of the atomic positions
                Call r3mv(ainv, atposc(:, ia, is), &
               &  input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:))
            end do
        end do

!-----------------------!
!       GROUNDSTATE
!-----------------------!

! check for overlapping muffin-tins
        Call checkmt
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
!       SCF cycle
!-----------------------!
        
        task = 1
        call scf_cycle

!-----------------------!
!       OUTPUT
!-----------------------!

        f = engytot
        j = 0
        do is = 1, nspecies
            do ia = 1, natoms(is)
                ias = idxas(ia,is)
                do i = 1, 3
                    j = j+1
                    g(j) = -forcetot(i,ias)
                end do
            end do
        end do
      
        return
    end subroutine calcEnergyForces
      
end subroutine
