!BOP
!
!!ROUTINE: calcevalqp
!
!!INTERFACE:
!
subroutine calcevalqp
!
!!DESCRIPTION:
! 
! Given the matrix elements $\langle
! \Psi_{n\vec{k}}|\Sigma(\vec{k},\omega)|\Psi_{n\vec{k}}\rangle$, $\langle
! \Psi_{n\vec{k}}|V^{xc}|\Psi_{n\vec{k}}\rangle$ and
! $\varepsilon^{DFT}_{n\vec{k}}$. this subroutine calculates the
! quasi-particle energies $\varepsilon^{qp}_{n\vec{k}}$
!
!!USES:
    use modinput
    use modmain, only : evalsv, efermi, zzero
    use modgw,   only : ibgw, nbgw, kset, evalqp, eferqp, &
    &                   sigc, znorm, selfex, selfec, vxcnn, &
    &                   nbandsgw, nvelgw, &
    &                   sigsx, sigch, fgw
       
    implicit none
    integer :: nb, ie, ik
    real(8) :: egap, df
 
!!REVISION HISTORY:
!
! Created: 16.08.05 by RGA
! Revisited June 2011 by DIN
!
!EOP
!BOC
    select case (input%gw%taskname)
      
      case('g0w0','gw0','acon')
        call calcevalqp_gw
         
      case('g0w0_x')
      
        do ik = 1, kset%nkpt
          do ie = ibgw, nbgw
            evalqp(ie,ik) = evalsv(ie,ik)+ &
            &  dble(selfex(ie,ik))-dble(vxcnn(ie,ik))
          end do ! ie
        end do ! ik
          
      case('cohsex')
      
        do ik = 1, kset%nkpt
          do ie = ibgw, nbgw
            sigsx(ie,ik) = sigsx(ie,ik)+selfex(ie,ik)
            !evalqp(ie,ik) = evalsv(ie,ik)+ &
            !&  dble(selfex(ie,ik))+dble(selfec(ie,1,ik))- &
            !&  dble(vxcnn(ie,ik))
            evalqp(ie,ik) = evalsv(ie,ik)+ &
            &  dble(sigsx(ie,ik))+dble(sigch(ie,ik))- &
            &  dble(vxcnn(ie,ik))
          end do ! ie
        end do ! ik
                  
    end select

    ! Calculate Fermi energy
    call fermi_exciting(input%groundstate%tevecsv, &
    &                   nvelgw, &
    &                   nbandsgw,kset%nkpt,evalqp(ibgw:nbgw,:), &
    &                   kset%ntet,kset%tnodes,kset%wtet,kset%tvol, &
    &                   eferqp,egap,df)  

contains

    subroutine calcevalqp_gw
    
    ! A critical parameters used in this subroutine is iopes, 
    ! which control how Fermi energy shift is treated 
    ! input%gw%selfenergy%iopes =    
    !          0 -- perturbative G0W0 without energy shift
    !          1 -- perturbative G0W0 with energy shift
    !          2 -- iterative G0W0 with energy shift
    !          3 -- iterative G0W0 without energy shift
    !         -1 -- selfconsitent GW0 without energy shift
    
        use modmain, only : zzero
        use modgw,   only : freq, evalks, evalqp, eferqp, nvelgw, &
        &                   nomax, numin, hev, selfex, selfec, &
        &                   iopac, vxcnn
       
        implicit none
        integer(4) :: ie 
        integer(4) :: ikp
        integer(4) :: npar ! Number of parameters of the analitic cont. of
                           ! the selfenergy
        integer(4) :: it, ierror, nb
        real(8)    :: enk, enk0, snk, vxcnk, znk
        real(8)    :: delta 
        real(8)    :: es, egap0, egapold
        complex(8) :: ein, dsig, sigma

        ! parameters of the functions fitting the selfenergy
        complex(8), allocatable :: a(:), sc(:), poles(:)
        complex(8), allocatable :: sacpar(:,:,:)

        real(8), parameter :: etol=1.0d-6
        integer, parameter :: nitmax=100
        
        write(fgw,*)
        select case(input%gw%selfenergy%iopes)
          case(0)
            write(fgw,*) 'Perform perturbative G0W0 without energy shift'
          case(1)
            write(fgw,*) 'Perform perturbative G0W0 with energy shift'
          case(2)
            write(fgw,*) 'Perform iterative G0W0 without energy shift'
          case(3)
            write(fgw,*) 'Perform iterative G0W0 with energy shift'
        end select
        write(fgw,*)
      
        ! initialize QP energies with KS values
        evalqp(ibgw:nbgw,1:kset%nkpt) = evalsv(ibgw:nbgw,1:kset%nkpt)
      
        ! Parameters of perform the analytic function of the correlation self-energy
        npar = 2*input%gw%selfenergy%npol
        allocate(a(npar),sc(freq%nomeg),poles(npar))
        allocate(sacpar(npar,ibgw:nbgw,kset%nkpt))

        !----------------------------------
        ! Start iterative procedure
        !----------------------------------      
        es = 0.0d0
        egap = 0.0d0
        ierror = 0
        do it = 0, nitmax
          do ikp = 1, kset%nkpt
            do ie = ibgw, nbgw
            
              ! energy from previous iteration 
              enk  = evalsv(ie,ikp)-efermi
              ! KS energy
              enk0 = evalks(ie,ikp)-efermi
              
              !---------------------------------------------
              ! Set up the analytic continuation parameters 
              !---------------------------------------------
              if (it==0) then
                sc(1:freq%nomeg) = selfec(ie,1:freq%nomeg,ikp)
                call setsac(iopac,freq%nomeg, &
                &           npar,enk,sc,freq%freqs,a,poles)
                sacpar(:,ie,ikp) = a
              else
                a = sacpar(:,ie,ikp)
              end if
              
              !--------------------------
              ! Apply the energy shift
              !--------------------------
              if (input%gw%selfenergy%iopes==2) then 
                ein = cmplx(evalqp(ie,ikp),0.0d0,8)
              else if (input%gw%selfenergy%iopes==3) then 
                ein = cmplx(evalqp(ie,ikp)-es,0.0d0,8)
              else 
                ein = cmplx(enk,0.0d0,8)
              endif
              
              !---------------------------------------
              ! Calculate the correlation self-energy   
              !---------------------------------------
              call getsac(iopac,freq%nomeg, &
              &           npar,enk,ein,freq%freqs,a,sigma,dsig)
              
              ! Set the correlation self-energy
              sigc(ie,ikp) = sigma

              ! Set the renormalization factor              
              znk = 1.0d0/(1.0d0-real(dsig))
              if ((znk>1.d0) .or. (znk<0.5d0)) then
                write(fgw,*) 'WARNING(calcevalqp):'
                write(fgw,100) ikp, ie, enk, znk, sigma, dsig
                100 format(' Suspicious Znk',' irk=',i4,' ie=',i4, &
                &          ' enk=',f8.3,' eV',' Znk=',f6.2," ReSc=",f8.3,  &
                &          ' ImSc=',f8.3," ReSc'=",f8.3," ImSc'=",f8.3)
                write(fgw,*)
              endif
              znorm(ie,ikp) = znk
              
              !------------------------------------------            
              ! Calculate the new quasi-particle energy 
              !------------------------------------------
              snk = dble(selfex(ie,ikp))+dble(sigc(ie,ikp))
              vxcnk = dble(vxcnn(ie,ikp))
              
              select case(input%gw%selfenergy%iopes)
              
                case(-2) ! no renormalization
                  delta = snk-vxcnk
              
                case(-1) ! self-consistent GW0
                  delta = znk*(snk-vxcnk+enk0-enk)

                case(0)
                  delta = znk*(snk-vxcnk)

                case(1) 
                  delta = znk*(snk-vxcnk)+(1.d0-znk)*es

                case(2) 
                  delta = snk-vxcnk

                case(3) 
                  if (it==0) then
                    delta = znk*(snk-vxcnk)
                  else
                    delta = snk-vxcnk
                  endif
                  
              end select 
              evalqp(ie,ikp) = enk+delta
            
          enddo ! ie
        enddo ! ikp
        
        ! if non-iterative scheme is applied
        if (input%gw%selfenergy%iopes<2) exit 
        
        ! Calculate Fermi energy
        call fermi_exciting(input%groundstate%tevecsv, &
        &                   nvelgw, &
        &                   nbandsgw,kset%nkpt,evalqp(ibgw:nbgw,:), &
        &                   kset%ntet,kset%tnodes,kset%wtet,kset%tvol, &
        &                   eferqp,egap,df)

        if (it==0) then 
          egap0 = egap
          write(fgw,*) ' # iter    Ef_QP    es    Eg/eV    Eg(k=0)/eV)' 
        endif
       
        if (egap>0.d0) then                                                
          write(fgw,'(i4,4f12.6)') it, eferqp, es, egap, &
          &  evalqp(numin,1)-evalqp(nomax,1)
        else                                                                 
          write(fgw,*) ' WARNING(calcevalqp): Metallic behaviour!'
          write(fgw,*) ' DOS at Fermi level', df                
        endif 

        if (it.ne.0) then 
          if (abs(egap-egap0)>0.2d0) then
            write(fgw,*) 'WARNING(calcevalqp): Band gap deviates from initial &
            &GW gap more than 0.2 Ha'
            ierror = 1
            exit
          end if
        end if
        
        ! convergence check
        if ((abs(es-eferqp)<etol).and.(abs(egap-egapold)<=etol)) exit

        ! Perform next iteration
        es = eferqp
        egapold = egap
         
      enddo ! it
      
      if (it>nitmax) ierror = 1
      if (ierror.ne.0) write(fgw,*) 'WARNING(calcevalqp): Failed to converge!'
      
      if (allocated(a)) deallocate(a)
      if (allocated(sc)) deallocate(sc)
      if (allocated(poles)) deallocate(poles)
      if (allocated(sacpar)) deallocate(sacpar)

      ! *** New: Return absolute values of the QP energies
      ! evalqp(:,:) = evalqp(:,:)+efermi
      
      return
    end subroutine

end subroutine
!EOC        
          
