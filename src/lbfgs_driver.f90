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
!-----------------------------------------------------------------------------80.
! This file is distributed under the terms of the GNU General Public License.
! Last modified on 14-06-2014 Pasquale Pavone (exciting team)
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
      real(dp)               :: v(3), forcesave
      integer                :: i, j
      character(1024)        :: message
      integer, allocatable   :: amap(:,:)
      integer                :: nscf, nconf, ncheckconv
      logical                :: lstart

!________________
! Starting checks

      if (istep+1>input%relax%maxsteps) return

!______________________________________________________________________________
! save initial configuration and forces in atposc_1 and forcetp, respecctively

      atposc_1(:,:,:) = atposc(:,:,:)
      forcetp(:,:) = forcetot(:,:)
      forcesave = forcemax

      if (rank==0) then
          write(string,'("Optimization step ", I4,"    (method = bfgs)")') istep+1
          call printbox(60,"-",string)
          call flushifc(60)
      end if

      call timesec(tsec1)

      if (input%groundstate%epsengy/input%relax%epsforce .gt. 0.020001) then
          input%groundstate%epsengy = max(input%relax%epsforce*0.02,1.d-9)
          if (rank==0) then
              write(60,'(" Convergence target for the total energy decreased to ",G13.6," Ha")') &
             & input%groundstate%epsengy
              write(60,*)
              call flushifc(60)
          end if
      end if

!________________________________________________________________________
! Initialize some L-BFGS-B library parameters (see src/Lbfgsb.3.0/README)
  
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
            if (.not.input%structure%speciesarray(is)%species%atomarray(ia)%atom%lockxyz(i)) n = n+1
          end do
        end do
      end do

      if (n==0) then
        call warning(' ')
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

!!!!!!INITIALIZE the searching loop

      ncheckconv = 0

99    continue

      nscf = 0
      nconf = 0
      ctask = 'START'
      lstart = .True.

!_______________________
! Set search constraints

     j = 0
     do is = 1, nspecies
        do ia = 1, natoms(is)
          do i = 1, 3
            if (.not.input%structure%speciesarray(is)%species%atomarray(ia)%atom%lockxyz(i)) then
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

!!!!!!BEGIN the searching loop

      do while ( (ctask(1:5).eq.'START'  .or. &
                  ctask(1:2).eq.'FG'     .or. &
                  ctask(1:5).eq.'NEW_X') )

!______________________________________     
! This is the call to the L-BFGS-B code
      
        call setulb(n,m,x,l,u,nbd,f,g,factr,pgtol,wa,iwa,ctask,iprint, &
        &           csave,lsave,isave,dsave)

        !write(60,'(" ctask = ",A)') trim(ctask)
        !call flushifc(60)

!_____________________________________________________
! ctask(1:2) .eq. 'FG'
! The minimization routine has returned to request the
! function f and gradient g values at the current x

        if (ctask(1:2) .eq. 'FG') then
          if (nconf .ge. input%relax%maxbfgs) then
            if (rank==0) then 
              call warning(' ')
              call warning('Warning(lbfgs_driver):')
              call warning(' Reached maximum number of investigated configurations')
              write(message,'(" for a single BFGS relaxation step :         ",I3,2X,A4)') nconf, ctask(1:4)
              call warning(message)
            end if 
            ctask(1:4) = 'NCFG'
          else
            call calcEnergyForces
          end if
        end if
         
!_______________________________________________
! ctask(1:5) .eq. 'NEW_X'
! The minimization has found a new configuration 
! to be used in the next optimization step

        if (ctask(1:5) .eq. 'NEW_X') then 

          istep = istep+1
          if (lstart) then
              nconf = nconf-1
              lstart = .False.
          end if

!______________________________________________________________________________
! save accepted configuration and forces in atposc_1 and forcetp, respecctively

          atposc_1(:,:,:) = atposc(:,:,:)
          forcetp(:,:) = forcetot(:,:)
          forcesave = forcemax

! output info

          if (rank==0) then

              if (input%relax%outputlevelnumber>1) write(60,*)
              write(60,'(" Number of investigated configurations",T45,": ",I9)') nconf
              write(60,'(" Number of total scf iterations",T45,": ",I9)') nscf
              write(60,'(" Maximum force magnitude",T36,"(target) : ",F18.8,"  (", F12.8, ")")') &
             &  forcemax, input%relax%epsforce
              write(60,'(" Total energy at this optimization step",T45,": ",F18.8)') engytot
              if (input%relax%outputlevelnumber>0) then 
                  call writepositions(60,input%relax%outputlevelnumber) 
                  call writeforce(60,input%relax%outputlevelnumber)                    
              end if
              if (input%relax%outputlevelnumber>1) call writechg (60,input%relax%outputlevelnumber)          
              if (input%relax%printtorque) call writetorque(60)          
              call flushifc(60)

!_____________________________________________________________
! write lattice vectors and optimised atomic positions to file

              Call writehistory
              Call writegeometryxml(.True.)

!__________________________________________________
! write the optimized interatomic distances to file

              Call writeiad(.True.)
          end if

!_______________________________________________
! write the time spent in this optimization step 

          call timesec(tsec2)
          if (rank==0) then
              write(60,*)
              write(60,'(" Time spent in this optimization step",T45,": ",F12.2," seconds")') tsec2-tsec1
              call flushifc(60)
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
              call flushifc(60)
          end if

          call timesec(tsec1)

        end if !!!!!!!! 'NEW_X'

!____________________________________________________________________________________
! ctask(1:4) .eq. 'CONV'
! The mininum configuration is likely outside the range of variability of parameters, 
! the searching loop is re-initialized (unless the number of re-initialization is 
! greater then 20

        if (ctask(1:4) .eq. 'CONV') then 
            ncheckconv = ncheckconv+1
            if (ncheckconv .le. 20) then

!__________________________________________________
! restart bfgs from the last accepted configuration

                atposc(:,:,:) = atposc_1(:,:,:)
                forcetot(:,:) = forcetp(:,:)
                forcemax = forcesave
                if (rank==0) then
                    call warning(' ')
                    call warning('Warning(lbfgs_driver):')
                    call warning(' Mininum configuration is likely to lie outside the searching region.')
                    write(message,'(" Restarting BFGS at step ", I3)') istep
                    call warning(message) 
                    write(message,'(" ctask = ",A)') trim(ctask)
                    call warning(message) 
                end if
                goto 99
            end if 
        end if

!_________________________________________________________________
! ctask(1:4) .eq. 'ABNO'
! The searching loop is broken, end with the current configuration

      end do

!!!!!!END the searching loop

!________________________________________________________
! Use Newton or harmonic method if BFGS does not converge

      if ((ctask(1:4).eq.'CONV') .or. (ctask(1:4).eq.'ABNO') .or. (ctask(1:4).eq.'NCFG')) then
        istep = istep+1

!___________________________________________________________________________________________
! either reset configuration to the last accepted one or 
! accept new configuration if the forcemax il lower than the value for the last accepted one

        if ( forcemax .gt. (forcesave-10*input%relax%epsforce) ) then
          atposc(:,:,:) = atposc_1(:,:,:)
          forcetot(:,:) = forcetp(:,:)
          forcemax = forcesave
        else
          if (rank==0) then
            if (input%relax%outputlevelnumber>1) write(60,*)
            write(60,'(" Warning: Accepted configuration has higher energy than previous one")')
            write(60,*)
            write(60,'(" Number of investigated configurations",T45,": ",I9)') nconf
            write(60,'(" Number of total scf iterations",T45,": ",I9)') nscf
            write(60,'(" Maximum force magnitude",T36,"(target) : ",F18.8,"  (", F12.8, ")")') &
           &  forcemax, input%relax%epsforce
            write(60,'(" Total energy at this optimization step",T45,": ",F18.8)') engytot
            if (input%relax%outputlevelnumber>0) then 
              call writepositions(60,input%relax%outputlevelnumber) 
              call writeforce(60,input%relax%outputlevelnumber)                    
            end if
            if (input%relax%outputlevelnumber>1) call writechg (60,input%relax%outputlevelnumber)          
            if (input%relax%printtorque) call writetorque(60)          
            call flushifc(60)
            Call writehistory
            Call writegeometryxml(.True.)
            Call writeiad(.True.)
            call timesec(tsec2)
            if (rank==0) then
              write(60,*)
              write(60,'(" Time spent in this optimization step",T45,": ",F12.2," seconds")') tsec2-tsec1
              call flushifc(60)
            end if
            if (forcemax <= input%relax%epsforce) return
            if (rank==0) then
              write(string,'("Optimization step ", I4,"    (method = bfgs)")') istep+1
              call printbox(60,"-",string)
              call flushifc(60)
            end if
            istep = istep+1
          end if
        end if

        if (input%relax%endbfgs.eq.'stop') then

          if (rank .Eq. 0) then
            write(60,'(" Number of investigated configurations",T45,": ",I9)') nconf
            write(60,*)
            write(60,'(A)') " BFGS scheme not converged -> Stopping BFGS"
            write(60,*)
            call flushifc(60)
            lstep = .True.
            call warning(' ')
            call warning('Warning(lbfgs_driver):')
            call warning(' BFGS scheme not converged')
            write(message,'(" ctask = ",A)') trim(ctask)
            call warning(message) 
            write(message,'(" -> Stopping BFGS at optimization step ",I3)') istep
            call warning(message)
          end if

        else

          if (rank .Eq. 0) then
            write(60,'(" Number of investigated configurations",T45,": ",I9)') nconf
            write(60,*)
            write(60,'(3A)') " BFGS scheme not converged -> Switching to ", trim(input%relax%endbfgs), " method"
            write(60,*)
            call flushifc(60)
            lstep = .True.
            call warning(' ')
            call warning('Warning(lbfgs_driver):')
            call warning(' BFGS scheme not converged')
            write(message,'(" ctask = ",A)') trim(ctask)
            call warning(message) 
            write(message,'(" -> Switching to ",A," method at optimization step ",I3)') &
           & trim(input%relax%endbfgs), istep
            call warning(message) 
          end if
          if (input%relax%endbfgs.eq.'harmonic') then
            call harmonic(input%relax%epsforce)
          else 
            call newton(input%relax%epsforce)
          end if

        end if
      end if

!____________________________________________________________________________________________
! Stop if ctask(1:5). eq. 'ERROR' : the routine has detected an error in the input parameters

      if (ctask(1:5).eq.'ERROR') then
        if (rank .Eq. 0) Then
          write(6,*)
          write(60,*)
          call printline(60,"#")
          call printline(6,"#")
          write(string,'("ERROR(lbfgs_driver): Routine has detected an error in the input parameters")') 
          call printtext(60,"#",string)
          call printtext(6,"#",string)  
          call printline(60,"#")
          call printline(6,"#")
          write(6,*)
          call flushifc(60)
          call flushifc(6)
          stop
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

!_______________________
! restart initialization

        Call init_relax

!__________
! SCF cycle

        call scf_cycle(-1)
        nconf = nconf + 1
        if ((rank==0).and.(input%relax%outputlevelnumber>1))  then 
            write(60,'(" Investigating configuration",T45,"# ",I9,"  (# of SCF cicles =",I4,")")') &
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

                  !call writepositions(60,input%relax%outputlevelnumber) 
                  !call writeforce(60,input%relax%outputlevelnumber)        

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

    end subroutine updatepositions
      
end subroutine
