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
!     Exciting interface created by DIN (31.01.2013)
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
            j = j+1
            nbd(j) = 0
            x(j) = atposc(i,ia,is)
            l(j) = atposc(i,ia,is)-input%structureoptimization%tau0atm
            u(j) = atposc(i,ia,is)+input%structureoptimization%tau0atm
          end do
        end do
      end do

!        ------- the beginning of the loop ----------
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
            Call writegeom (.True.)
            ! write the optimised interatomic distances to file
            Call writeiad (.True.)

            Write (60,*)
            Write (60, '("+--------------------------+")')
            Write (60, '("| Updated atomic positions |")')
            Write (60, '("+--------------------------+")')
            
            j = 0
            do is = 1, nspecies
              Write (60,*)
              Write (60, '("Species : ", I4, " (", A, ")")') is, trim (input%structure%speciesarray(is)%species%chemicalSymbol)
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
               Return
            End If
            
          end if ! 'NEW_X'
          
        end if

      end do
!           ---------- the end of the loop -------------

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
    
    ! local variables
    Logical :: exist, tibs
    Integer :: ik, is, ia, idm, i, j
    Integer :: n, nwork
    Real(8) :: et, fm, timetot
    Real(8) :: deltae, dforcemax
    
    ! allocatable arrays
    Real(8),    Allocatable :: v (:)
    Real(8),    Allocatable :: evalfv (:, :)
    Complex(8), Allocatable :: evecfv (:, :, :)
    Complex(8), Allocatable :: evecsv (:, :)
    
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
 !     GROUNDSTATE       !
 !-----------------------!
    task = 1

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

! size of mixing vector
    n = lmmaxvr*nrmtmax*natmtot+ngrtot
    If (associated(input%groundstate%spin)) n = n*(1+ndmag)
    If (ldapu .Ne. 0) n = n+2*lmmaxlu*lmmaxlu*nspinor*nspinor*natmtot
! allocate mixing arrays
    if(allocated(v))deallocate(v)
    Allocate (v(n))

! set stop flag
    tstop = .False.
    
    If (rank .Eq. 0) Then
       Write (60,*)
       Write (60, '("+-----------------------+")')
       Write (60, '("| L-BFGS-B Optimization |")')
       Write (60, '("+-----------------------+")')
    End If
     
! call mixing array allocation functions by setting
    nwork = - 1
! and call interface
    If (rank .Eq. 0) Call mixerifc (input%groundstate%mixernumber, n, v, currentconvergence, nwork)
    et = 0.d0
    fm = 0.d0

! set last iteration flag
    tlast = .False.

! delete any existing eigenvector files
    If ((rank .Eq. 0) .And. ((task .Eq. 0) .Or. (task .Eq. 2))) Call delevec

!-----------------------------------
! begin the self-consistent loop
!-----------------------------------
    If (rank .Eq. 0) Then
       Write (60,*)
       Write (60, '("+------------------------------+")')
       Write (60, '("| Self-consistent loop started |")')
       Write (60, '("+------------------------------+")')
    End If
    Do iscl = 1, input%groundstate%maxscl
       If (rank .Eq. 0) Then
          Write (60,*)
          Write (60, '("+-------------------------+")')
          Write (60, '("| Iteration number : ", I4, " |")') iscl
          Write (60, '("+-------------------------+")')
       End If
       If (iscl .Ge. input%groundstate%maxscl) Then
          If (rank .Eq. 0) Then
             Write (60,*)
             Write (60, '("Reached self-consistent loops maximum")')
             Write (100,*)
             Write (100, '("Warning(lbfgs_driver): Reached self-consistent loops maximum")')
          End If
          tlast = .True.
       End If
       If (rank .Eq. 0) Call flushifc (60)
       ! generate the core wavefunctions and densities
       Call gencore
       ! find the new linearisation energies
       Call linengy
       ! write out the linearisation energies
       if (rank .eq. 0) Call writelinen
       ! generate the APW radial functions
       Call genapwfr
       ! generate the local-orbital radial functions
       Call genlofr
       ! compute the overlap radial integrals
       Call olprad
       ! compute the Hamiltonian radial integrals
       Call hmlrad
#ifdef MPI
       Call MPI_barrier (MPI_COMM_WORLD, ierr)
       If (rank .Eq. 0) Call delevec ()
#endif

!------------------------------------
! begin parallel loop over k-points
!------------------------------------

#ifdef MPISEC
       splittfile = .True.
       Do ik = firstk (rank), lastk (rank)
!
#endif
#ifndef MPISEC
       splittfile = .False.
       Do ik = 1, nkpt
#endif
         ! every thread should allocate its own arrays
         Allocate (evalfv(nstfv, nspnfv))
         Allocate (evecfv(nmatmax, nstfv, nspnfv))
         Allocate (evecsv(nstsv, nstsv))
         ! solve the first- and second-variational secular equations
         Call seceqn (ik, evalfv, evecfv, evecsv)
         ! write the eigenvalues/vectors to file
         Call putevalfv (ik, evalfv)
         Call putevalsv (ik, evalsv(:, ik))
         Call putevecfv (ik, evecfv)
         Call putevecsv (ik, evecsv)
         Deallocate (evalfv, evecfv, evecsv)
       End Do
#ifdef MPISEC
       call mpi_allgatherv_ifc(nkpt,nstsv,rbuf=evalsv)
       Call MPI_barrier (MPI_COMM_WORLD, ierr)
#endif
       ! find the occupation numbers and Fermi energy
       Call occupy
       If (rank .Eq. 0) Then
         ! write out the eigenvalues and occupation numbers
         Call writeeval
         ! write the Fermi energy to file
         Call writefermi
       End If
       ! set the charge density and magnetisation to zero
       rhomt (:, :, :) = 0.d0
       rhoir (:) = 0.d0
       If (associated(input%groundstate%spin)) Then
          magmt (:, :, :, :) = 0.d0
          magir (:, :) = 0.d0
       End If
#ifdef MPIRHO
       Do ik = firstk (rank), lastk (rank)
         !write the occupancies to file
          Call putoccsv (ik, occsv(:, ik))
       End Do
       Do ik = firstk (rank), lastk (rank)
#endif
#ifndef MPIRHO
       If (rank .Eq. 0) Then
         Do ik = 1, nkpt
         !write the occupancies to file
           Call putoccsv (ik, occsv(:, ik))
         End Do
       End If
       Do ik = 1, nkpt
#endif
         Allocate (evecfv(nmatmax, nstfv, nspnfv))
         Allocate (evecsv(nstsv, nstsv))
         ! get the eigenvectors from file
         Call getevecfv (vkl(:, ik), vgkl(:, :, :, ik), evecfv)
         Call getevecsv (vkl(:, ik), evecsv)
         ! add to the density and magnetisation
         Call rhovalk (ik, evecfv, evecsv)
         Deallocate (evecfv, evecsv)
       End Do
#ifdef MPIRHO
       Call mpisumrhoandmag
#endif
       ! symmetrise the density
       Call symrf (input%groundstate%lradstep, rhomt, rhoir)
       ! symmetrise the magnetisation
       If (associated(input%groundstate%spin)) Call symrvf (input%groundstate%lradstep, magmt, magir)
       ! convert the density from a coarse to a fine radial mesh
       Call rfmtctof (rhomt)
       ! convert the magnetisation from a coarse to a fine radial mesh
       Do idm = 1, ndmag
         Call rfmtctof (magmt(:, :, :, idm))
       End Do
       ! add the core density to the total density
       Call addrhocr
       ! calculate the charges
       Call charge
       ! calculate the moments
       If (associated(input%groundstate%spin)) Call moment
       ! normalise the density
       Call rhonorm
       ! LDA+U
       If (ldapu .Ne. 0) Then
         ! generate the LDA+U density matrix
         Call gendmatlu
         ! generate the LDA+U potential matrix
         Call genvmatlu
         ! write the LDA+U matrices to file
         if (rank .eq. 0) Call writeldapu
       End If
       ! generate charge distance
       call chgdist
       ! store density to reference
       rhoirref(:)=rhoir(:)
       rhomtref(:,:,:)=rhomt(:,:,:)
       ! compute the effective potential
       Call poteff
       ! pack interstitial and muffin-tin effective potential and field into one array
       Call packeff (.True., n, v)
       ! mix in the old potential and field with the new
       If (rank .Eq. 0) Call mixerifc(input%groundstate%mixernumber, n, v, currentconvergence, nwork)
#ifdef MPI
       Call MPI_bcast (v(1), n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
#endif
 ! unpack potential and field
       Call packeff (.False., n, v)
       ! add the fixed spin moment effect field
       If (getfixspinnumber() .Ne. 0) Call fsmfield
       ! Fourier transform effective potential to G-space
       Call genveffig
       ! reduce the external magnetic fields if required
       If (associated(input%groundstate%spin)) Then
          If (input%groundstate%spin%reducebf .Lt. 1.d0) Then
              input%groundstate%spin%bfieldc(:) = &
             &  input%groundstate%spin%bfieldc(:) * input%groundstate%spin%reducebf
             Do is = 1, nspecies
                Do ia = 1, natoms (is)
                   input%structure%speciesarray(is)%species%atomarray(ia)%atom%bfcmt(:) =   &
                  &  input%structure%speciesarray(is)%species%atomarray(ia)%atom%bfcmt(:) * &
                  &  input%groundstate%spin%reducebf
                End Do
             End Do
          End If
       End If
       ! compute the energy components
       Call energy
       tibs=input%groundstate%tfibs
       input%groundstate%tfibs=.false.
       call force
       input%groundstate%tfibs=tibs
       If (rank .Eq. 0) Then
         ! output forces to INFO.OUT
         Call writeforce (60)
         ! output energy components
         Call writeengy (60)
         Write (60,*)
         Write (60, '("Density of states at Fermi energy : ", G18.10)') fermidos
         Write (60, '(" (states/Hartree/unit cell)")')
         ! output charges and moments
         Call writechg (60)
         ! output effective fields for fixed spin moment calculations
         If (getfixspinnumber() .Ne. 0) Call writefsm (60)
         ! check for WRITE file
         Inquire (File='WRITE', Exist=Exist)
         If (exist) Then
            Write (60,*)
            Write (60, '("WRITE file exists - writing STATE.OUT")')
            Call writestate
            Open (50, File='WRITE')
            Close (50, Status='DELETE')
         End If
         ! write STATE.OUT file if required
         If (input%groundstate%nwrite .Ge. 1) Then
            If (Mod(iscl, input%groundstate%nwrite) .Eq. 0) Then
               Call writestate
               Write (60,*)
               Write (60, '("Wrote STATE.OUT")')
            End If
         End If
         ! update convergence criteria
         deltae=abs(et-engytot)
         dforcemax=abs(fm-forcemax)
       End If
       ! exit self-consistent loop if last iteration is complete
       If (tlast) Go To 20
       If (rank .Eq. 0) Then
         ! check for convergence
         If (iscl .Ge. 2) Then
            Write (60,*)
            Write (60, '("RMS change in effective potential (target) : ", G18.10, " (", G18.10, ")")') &
           &  currentconvergence, input%groundstate%epspot
            write(60,'("Absolute change in total energy (target)   : ",G18.10," (",G18.10,")")') &
           &  deltae, input%groundstate%epsengy
            write(60,'("Absolute change in |max. force| (target)   : ",G18.10," (",G18.10,")")') &
           &  dforcemax, input%groundstate%epsforce
            write(60,'("Charge distance (target)                   : ",G18.10," (",G18.10,")")') &
           &  chgdst, input%groundstate%epschg
            If ((currentconvergence .Lt. input%groundstate%epspot).and. &
           &    (deltae .lt. input%groundstate%epsengy).and. &
           &    (chgdst .lt. input%groundstate%epschg).and. &
           &    (dforcemax .lt. input%groundstate%epsforce)) Then
               Write (60,*)
               Write (60, '("Convergence targets achieved")')
               tlast = .True.
            End If
         End If
         et = engytot
         fm = forcemax
         ! check for STOP file
         Inquire (File='STOP', Exist=Exist)
         If (exist) Then
            Write (60,*)
            Write (60, '("STOP file exists - stopping self-consistent loop")')
            tstop = .True.
            tlast = .True.
            Open (50, File='STOP')
            Close (50, Status='DELETE')
         End If
       End If
#ifdef MPI
       Call MPI_bcast (tstop, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
       Call MPI_bcast (tlast, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
#endif
     End Do ! iscl

!-----------------------------------
! end the self-consistent loop
!-----------------------------------

20   Continue
     If (rank .Eq. 0) Then
        Write (60,*)
        Write (60, '("+------------------------------+")')
        Write (60, '("| Self-consistent loop stopped |")')
        Write (60, '("+------------------------------+")')
 ! write density and potentials to file only if maxscl > 1
        If (input%groundstate%maxscl .Gt. 1) Then
           Call writestate
           Write (60,*)
           Write (60, '("Wrote STATE.OUT")')
        End If
     end if

 !-----------------------!
 !     compute forces    !
 !-----------------------!

     If ( .Not. tstop) Then
        Call force
        If (rank .Eq. 0) Then
          ! output forces to INFO.OUT
          Call writeforce (60)
        End If
     End If

     !set nwork to -2 to tell interface to call the deallocation functions
     If (rank .Eq. 0) Call mixerifc (input%groundstate%mixernumber, n, v, currentconvergence, -2)
     Deallocate (v)
     
     Call mpiresumeevecfiles ()

 !-----------------------!
 !        OUTPUT         !
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
