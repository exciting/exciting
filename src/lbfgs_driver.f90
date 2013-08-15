!
!-----------------------------------------------------------------------------80
!                                        
! CUSTOMIZED DRIVER FOR L-BFGS-B
!
! L-BFGS-B is a code for solving large nonlinear optimization
! problems with simple bounds on the variables.
!
! The code can also be used for unconstrained problems and is as 
! efficient for these problems as the earlier limited memory code L-BFGS.
!
! This driver illustrates how to control the termination of the
! run and how to design customized output.
!
! References:
!
! [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
!     memory algorithm for bound constrained optimization'',
!     SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.
!
! [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: FORTRAN
!     Subroutines for Large Scale Bound Constrained Optimization''
!     Tech. Report, NAM-11, EECS Department, Northwestern University, 1994.
!
! Postscript files of these papers are available via anonymous
! ftp to eecs.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.)
!
! February 2011   (latest revision)
! Optimization Center at Northwestern University
! Instituto Tecnologico Autonomo de Mexico
!
! Jorge Nocedal and Jose Luis Morales
!
!-----------------------------------------------------------------------------80

subroutine lbfgs_driver
!
!     Exciting interface created by DIN (February 2013)
!
      use modinput
      use modmain      
      Use modmpi

      implicit none
      
      integer,  parameter    :: dp = kind(1.0d0)

      real(dp)               :: factr, pgtol
      integer                :: n
      integer                :: m, iprint
      integer                :: is, ia, ias, ik, ispn
      
      character(len=77)      :: string
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
      integer                :: nscf, nconf
      logical                :: lnconv

      if (istep+1>input%relax%maxsteps) return

      if (rank==0) then
          write(string,'("Optimization step ", I4,"    (method = bfgs)")') istep+1
          call printbox(60,"-",string)
      end if

!_______________________________________________________________________
! Initialize some L-BFGS-B library parameters (see src/Lbfgs.3.0/README)
  
      m = 5         ! number of corrections used in the limited memory matrix
      iprint = -1   ! controls the frequency and type of output generated
      factr = 0.d0  ! suppress termination test controlled by machine precision
      pgtol = 0.d0  ! suppress termination test controlled by the component 
                    ! of the projected gradient
      
!____________________________________________________________________
! Determine total number of variables taking into account constraints

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

!________________
! Allocate memory
      
      allocate( nbd(n), x(n), l(n), u(n), g(n) )
      allocate( iwa(3*n) )
      allocate( wa(2*m*n + 5*n + 11*m*m + 8*m) )
      allocate( amap(3,n) )

!_______________________
! Set search constraints

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
              l(j) = x(j)-input%relax%taubfgs
              u(j) = x(j)+input%relax%taubfgs
            end if
          end do
        end do
      end do

!!!!!!INITIALIZE the loop

      nscf = 0
      nconf = 0
      ctask = 'START'
      lnconv = .True.

!!!!!!BEGIN the loop

      do while ( (ctask(1:5).eq.'START'  .or. &
                  ctask(1:2).eq.'FG'     .or. &
                  ctask(1:5).eq.'NEW_X') .and. &
                  lnconv )

                                  !write(60,*) "loop starts ", ctask

!______________________________________     
! This is the call to the L-BFGS-B code

                                  !write(60,*) " before"
                                  !write(60,'(6f12.6)') x
      
        call setulb(n,m,x,l,u,nbd,f,g,factr,pgtol,wa,iwa,ctask,iprint, &
        &           csave,lsave,isave,dsave)

                                  !write(60,*) " after"
                                  !write(60,'(6f12.6)') x
                                  !write(60,*) ctask, nconf

!_____________________________________________________
! the minimization routine has returned to request the
! function f and gradient g values at the current x

        if (ctask(1:2) .eq. 'FG') then
          call calcEnergyForces
        else 
          
!_____________________________________________________
! the minimization has found a new configuration 
! to be used in the next optimization step

          if (ctask(1:5) .eq. 'NEW_X') then   
           
            call updatepositions
            istep = istep+1

!____________
! output info

            if (rank==0) then
                if (input%relax%outputlevelnumber>1)  write(60,*)
                write(60,'(" Number of investigated configurations  : ",I5)') nconf
                write(60,'(" Number of total scf iterations         : ",I5)') nscf
                write(60,'(" Maximum force magnitude       (target) : ",F14.8,"    (", F14.8, ")")') &
               &  forcemax, input%relax%epsforce
                write(60,'(" Total energy at this optimization step :",F19.9)') engytot
                if (input%relax%outputlevelnumber>0)  then 
                    call writepositions(60,input%relax%outputlevelnumber) 
                    call writeforce(60,input%relax%outputlevelnumber)                    
                end if
                call flushifc(60)

!_____________________________________________________________
! write lattice vectors and optimised atomic positions to file

!                Call writehistory(istep) 
                Call writehistory
                Call writegeometryxml(.True.)
!__________________________________________________
! write the optimized interatomic distances to file

                Call writeiad(.True.)
            end if
            
!______________________________________
! check if force convergence is reached

            if (forcemax <= input%relax%epsforce) ctask = 'STOP'

!_________________________________________________
! check if maximum number of iterations is reached

            if (istep    >= input%relax%maxsteps) ctask = 'STOP'
            
            nconf = 0
            nscf = 0

            if ((rank==0).and.(ctask(1:5).eq.'NEW_X')) then
                write(string,'("Optimization step ", I4,"    (method = bfgs)")') istep+1
                call printbox(60,"-",string)
            end if

            if (ctask(1:4).eq.'CONV') lnconv = .False.

          end if ! 'NEW_X'
          
                                  !write(60,*) "else  ", ctask 

        end if

                                  !write(60,*) "final ", ctask

      end do

!!!!!!END the loop

!________________________________________________________
! Use Newton or harmonic method if BFGS does not converge

      if (ctask(1:4).eq.'CONV') then
        istep = istep+1
        if (input%relax%endbfgs.eq.'harmonic') then
          string = "BFGS scheme not converged -> Switching to harmonic method"
        else 
          string = "BFGS scheme not converged -> Switching to newton method"
        end if
       if (rank .Eq. 0) then
          write(60,*)
          write(60,*) string
          write(60,*)
          call flushifc(60)
          lstep = .True.
        end if
        if (input%relax%endbfgs.eq.'harmonic') then
           call harmonic(input%relax%epsforce)
        else 
          call newton(input%relax%epsforce)
        end if
      end if
      
!__________________
! Deallocate memory

      deallocate( nbd, x, l, u, g )
      deallocate( iwa )
      deallocate( wa )

contains
!
!-----------------------------------------------------------------------------80
!
    subroutine calcEnergyForces
     
        implicit none
        integer :: is, ia, i, j

!________________________
! update atomic positions

        call updatepositions

!__________________________________
! check for overlapping muffin-tins

        call checkmt

!_________________________________________
! generate structure factors for G-vectors

        Call gensfacgp (ngvec, vgc, ngvec, sfacg)

!_____________________________________
! generate the characteristic function

        Call gencfun

!___________________________________________
! generate structure factors for G+k-vectors

        Do ik = 1, nkpt
            Do ispn = 1, nspnfv
                Call gensfacgp (ngk(ispn, ik), vgkc(:, :, ispn, ik), &
               &  ngkmax, sfacgk(:, :, ispn, ik))
            End Do
        End Do

!_________________________________________
! determine the new nuclear-nuclear energy

        Call energynn

!__________
! SCF cycle

        call scf_cycle(-1)
        nconf = nconf + 1
        if ((rank==0).and.(input%relax%outputlevelnumber>1))  then 
            write(60,'(" Investigating configuration            # ",I5,"    (# of SCF cicles =",I5,")")') &
           &      nconf, iscl
            call flushifc(60)
        end if
        nscf = nscf+iscl

!_______
! Output

        f = engytot                            ! total energy
        
        do j = 1, n
          i = amap(1,j); ia = amap(2,j); is = amap(3,j)
          ias = idxas(ia,is)
          g(j) = -forcetot(i,ias)              ! energy gradients = total forces
        end do
      
        return
    end subroutine calcEnergyForces
!
!-----------------------------------------------------------------------------80
!
    subroutine updatepositions
        
        implicit none
      
        do j = 1, n
            i = amap(1,j); ia = amap(2,j); is = amap(3,j)
            atposc(i,ia,is) = x(j)
        end do

!________________________________________________________
! compute the lattice coordinates of the atomic positions

        do is = 1, nspecies
            do ia = 1, natoms(is)
                call r3mv (ainv, atposc(:, ia, is), &
               &     input%structure%speciesarray(is)%species%atomarray(ia)%atom%coord(:))
            end do
        end do

!___________________________________________________________________
! find the crystal symmetries and shift atomic positions if required

        call findsymcrys
    
    end subroutine updatepositions
      
end subroutine
