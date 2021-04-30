  

!> Calculate the electronic transport coefficients
!>
!> Calculate the electronic transport coefficients, specifically:
!>  * the electronic conductivity
!>  * the Seebeck coefficient
!>  * the electronic part of the thermal conductivity.
!>
!> They are obtained by solving the linearized Boltzmann transport equation in the relaxation time approximation. 
!> TODO(Maria) Issue 55 Document equations 
!> 
!> For more information, see:
!>  - B.R. Nag, "Electron Transport in Compound Semiconductors" (Springer Verlag, Berlin, 1980), Chap. 7
!>  - G.D. Mahan and J.O Sofo, "The best thermoelectric", PNAS 93, 7436 (1996)

!TODO(Maria) Issue 55 Pass inputs and outputs
subroutine boltzequ()
    ! TODO(Maria) Issue 55 Remove use of modmain
    use modmain
    use modinput
    use modmpi
    use modxs, only: symt2
    use precision, only: dp
    use unit_conversion, only: hartree_to_ev, hartree_to_j
    use physical_constants, only: kboltz, hbar_si, elec_charge, bohr_radius_si
    implicit none

    integer :: l, a, b, ncomp, etCoeffComponents(2,9)
    integer :: ik, jk, isym
    integer :: ist, jst, iw, itemp, imu
    integer :: recl, iostat
    integer :: ncol, it
    integer :: n_mu_points, n_temp_points, n_tdf_points
    logical :: exist, tefermi, tdoping, ttdf
    real(dp) :: swidth, t1, t2, t3
    real(dp) :: v1 (3), v2 (3), sc(3,3)
    real(dp) :: chg, efermi2, e0, e1, x, nvm, charge, doping, e_fermi_shift, level
    real(dp) :: tdf_spacing, tdf_init
    real(dp) :: temp, temp_spacing, temp_init
    real(dp) :: mu, nsum, mu_spacing, mu_init
    complex(dp) :: zt1, zt2, zt3
    character(256) :: fname
    real(dp), allocatable :: wtdf(:)
    real(dp), allocatable :: temparray(:), muarray(:), mu_grid(:)
    real(dp), allocatable :: evalsvt(:)
    complex(dp), allocatable :: pmat(:,:,:) 

    !> conductivity 
    complex(dp), allocatable :: sigma(:,:)
    !> Seebeck coefficient
    complex(dp), allocatable :: seebeck(:,:)
    !> thermal conductivity
    complex(dp), allocatable :: thermalcond(:,:)
    !> ZT Figure of merit
    complex(dp), allocatable :: zt(:,:)
    !> Transport distribution function
    complex(dp), allocatable :: td(:)
    !> chemical potential adjusted to the doping level
    real(dp), allocatable :: chelevel(:)
    
    integer :: kfirst, klast    
    real(dp) sdelta, stheta

    !> Maximum number of iterations to find the position of the chemical potential
    integer, parameter :: maxit = 2000
    !> Boltzmann distribution exponent 
    real(dp), parameter :: critical_factor = 40.0

    
    ! initialise universal variables
    call init0
    call init1

    ! read Fermi energy from file
    call readfermi
    
    ! booleans from the input
    tefermi = .False.
    if (input%properties%boltzequ%energyReference == "efermi") then
       tefermi = .True.
    end if
    tdoping = input%properties%boltzequ%useDopingConcentration
    ttdf = input%properties%boltzequ%useTransportDf

    e_fermi_shift = 0
    if (tefermi) then
       e_fermi_shift = efermi
    end if

    ! TODO(Maria) Issue 55 Replace with calls to linear grid routine
    ! Mu and temperature grids 
    if (tdoping) then
       allocate(muarray(1))
       muarray(1) = input%properties%boltzequ%chemicalPotentialRange(1)
       
    else
       mu_spacing = input%properties%boltzequ%chemicalPotentialSpacing
       mu_init=input%properties%boltzequ%chemicalPotentialRange(1)
       if (mu_init == input%properties%boltzequ%chemicalPotentialRange(2) ) then
          n_mu_points = 1
       else
          n_mu_points = ceiling((input%properties%boltzequ%chemicalPotentialRange(2)-mu_init)/(1.0*mu_spacing))+1
       end if
       
       allocate(muarray(n_mu_points))
       Do imu = 1, size(muarray)
          muarray(imu) = mu_init + mu_spacing*(imu-1)
       End Do
       
       !mu_grid = linspace(input%properties%boltzequ%windmu(1), input%properties%boltzequ%windmu(2), n_mu_points)
       !if (mu_init == muf) then
       !   if (size(muarray) > 1) then
       !      deallocate(muarray)
       !      allocate(muarray(1))
       !   end if
       !   muarray(1) = mui
       !else
       !   mus=(muf-mui)/(size(muarray)-1)
       !   Do imu = 1, size(muarray)
       !      muarray(imu) = mui + mus*(imu-1)
       !   End Do
       !end if
    end if

    temp_spacing = input%properties%boltzequ%temperatureSpacing
    temp_init=input%properties%boltzequ%temperatureRange(1)
    if (temp_init == input%properties%boltzequ%temperatureRange(2) ) then
       n_temp_points = 1
    else
       n_temp_points = ceiling((input%properties%boltzequ%temperatureRange(2)-temp_init)/(1.0*temp_spacing))+1
    end if
    
    allocate(temparray(n_temp_points))
    Do itemp = 1, size(temparray)
       temparray(itemp) = temp_init + temp_spacing*(itemp-1)
    End Do

    ncol = size(muarray)*size(temparray)
  
    if (ttdf) then 
       tdf_spacing = input%properties%boltzequ%transportDfSpacing
       tdf_init=input%properties%boltzequ%transportDfRange(1)
       if (tdf_init == input%properties%boltzequ%transportDfRange(2) ) then
          n_tdf_points = 1
       else
          n_tdf_points = ceiling((input%properties%boltzequ%transportDfRange(2)-tdf_init)/(1.0*tdf_spacing))+1
       end if
       allocate(wtdf(n_tdf_points))       
       !wtdfi = input%properties%boltzequ%windtdf(1)
       !wtdff = input%properties%boltzequ%windtdf(2)
       !t1 = (wtdff-wtdfi)/dble(size(wtdf)-1)       
       Do iw = 1, size(wtdf)
          wtdf(iw) = tdf_init + tdf_spacing*(iw-1)
       End Do

    end if

    ! End of linear grid constructors 
    !------------------------------------
    
    allocate(td(size(wtdf)))
    allocate(sigma(size(muarray),size(temparray)))
    allocate(seebeck(size(muarray),size(temparray)))
    allocate(thermalcond(size(muarray),size(temparray)))
    allocate(zt(size(muarray),size(temparray)))
    ! KS eigenvalues and occupancies    
    allocate(evalsvt(nstsv))    
    call getevalsv(vkl(:,:),evalsv(:,:))
      
    if (tdoping) then
       !search for Fermi level at doping
       nvm  = nint(chgval/occmax)

       allocate(chelevel(size(temparray)))
       chelevel(:) = zzero
       charge=(chgval+((input%properties%boltzequ%dopingConcentration)*omega*(bohr_radius_si*100)**3))
    
       Do itemp = 1, size(temparray)
          temp = temparray(itemp)

          call getevalsv(vkl(:,1),evalsvt)
          e0 = evalsvt(1)
          e1 = e0
          Do ik = 1, nkpt
             call getevalsv(vkl(:,ik),evalsvt)
             Do ist = 1, nstsv
                e0 = Min (e0, evalsvt(ist))
                e1 = Max (e1, evalsvt(ist))
             End Do
          End Do

          t1 = 1.0_dp / (kboltz*temp)    
          Do it = 1, maxit
             efermi2 = 0.5_dp * (e0+e1)
             chg = 0._dp
             Do ik = 1, nkpt
                call getevalsv(vkl(:,ik),evalsvt)
                Do ist = 1, nstsv
                   x = (efermi2-evalsvt(ist)) * t1
                   occsv (ist, ik) = occmax * stheta(3, x)
                   chg = chg + wkpt (ik) * occsv (ist, ik)
                End Do
             End Do
             If (chg .Lt. charge) Then
                e0 = efermi2
             Else
                e1 = efermi2
             End If
             If ((e1-e0) .Lt. input%groundstate%epsocc) Exit
          End Do
          chelevel(itemp) = efermi2
       End Do
    end if

    ! smearing factor
    swidth = input%properties%boltzequ%transportDfBroadening

    ! the momentum matrix elements
    allocate(pmat(3,nstsv,nstsv))
    inquire(iolength=recl) pmat
    open(50,File='PMAT.OUT',Action='READ',Form='UNFORMATTED', &
    &    Access='DIRECT',Recl=recl,IOstat=iostat)

#ifdef MPI
    kfirst = firstofset(rank,nkptnr)
    klast = lastofset(rank,nkptnr)
#else
    kfirst = 1
    klast = nkptnr
#endif
 
    !-------------------------------------------------------------------------------
    !   Loop over optical components
    !-------------------------------------------------------------------------------
    ncomp = size(input%properties%boltzequ%etCoeffComponents,2)
    if (ncomp>0) then
        do l = 1, ncomp
            etCoeffComponents(1,l) = input%properties%boltzequ%etCoeffComponents(1,l)
            etCoeffComponents(2,l) = input%properties%boltzequ%etCoeffComponents(2,l)
        end do
    else

        ncomp = 0
        do a = 1, 3
           do b = 1, 3
              !TODO(Maria) Issue 55 Remove magic number 
                if (sum(abs(symt2(a,b,:,:)))>1.0d-8) then
                    ncomp = ncomp+1
                    etCoeffComponents(1,ncomp) = a
                    etCoeffComponents(2,ncomp) = b
                end if
            end do
        end do
    end if
    
    do l = 1, ncomp
        nsum = 0.0_dp

        ! sum over non-reduced k-points
        td = zzero
        sigma = zzero
        seebeck = zzero
        thermalcond = zzero
        zt = zzero
        !if (ttdf) then 
        !   td(:) = zzero
        !end if

        do ik = kfirst, klast


           ! equivalent reduced k-point
           call findkpt(vklnr(:,ik),isym,jk)
           
           ! read momentum matrix elements from direct-access file
           read(50,Rec=jk) pmat
           
           ! rotate the matrix elements from the reduced to non-reduced k-point
           if (isym > 1) then
               do ist = 1, nstsv
                  do jst = 1, nstsv
                     sc(:,:)=symlatc(:,:,lsplsymc(isym))
                     call r3mv(sc,dble(pmat(:,ist,jst)),v1)
                     call r3mv(sc,aimag(pmat(:,ist,jst)),v2)
                     pmat(:,ist,jst) = cmplx(v1(:),v2(:),8)
                  end do
               end do
             end if
            
             ! read eigenvalues and occupancies from files
             call getevalsv(vkl(:,jk),evalsvt)

             !--------------------------
             ! Transport distribution
             !--------------------------
            nsum=nsum+wkptnr(jk)
            a = etCoeffComponents(1,l)
            b = etCoeffComponents(2,l)

            If (ttdf) then
               do ist = 1, nstsv
                  zt1 = occmax*pmat(a,ist,ist)*conjg(pmat(b,ist,ist))
                  do iw = 1, size(wtdf)
                     t1 = (evalsvt(ist)-wtdf(iw)-e_fermi_shift)/swidth
                     !td(iw) = td(iw)+ wkptnr(jk)*zt1*sdelta(input%groundstate%stypenumber,t1)/swidth
                     td(iw) = td(iw)+ zt1*sdelta(input%groundstate%stypenumber,t1)/swidth
                  end do ! iw
               end do ! ist
               
            Else

               do imu = 1, size(muarray)
                  mu = muarray(imu)
                  
                  do itemp = 1, size(temparray)
                     temp = temparray(itemp)

                     if (tdoping) then
                        level = chelevel(itemp)
                     else
                        level = mu+e_fermi_shift
                     end if
                     
                     do ist = 1, nstsv
                        t1 = -(level-evalsvt(ist))/(temp*kboltz)

                        if (abs(t1) > critical_factor) then
                           t3=0
                        else
                           t2 = exp(t1)
                           t3 = t2/((t2+1)**2 * temp*kboltz)
                        end if
                        zt1=occmax*pmat(a,ist,ist)*conjg(pmat(b,ist,ist))
                        sigma(imu, itemp) = sigma(imu, itemp) + zt1*t3
                        seebeck(imu, itemp) = seebeck(imu, itemp) + zt1*t3*t1
                        thermalcond(imu, itemp) = thermalcond(imu, itemp) + zt1*t3*(t1**2)

                     end do ! ist
                  end do ! mu
               end do ! temp
               
            End If
         end do
                    
#ifdef MPI
        call MPI_AllReduce(MPI_IN_PLACE, &
        &                  td, &
        &                  size(wtdf), &  
        &                  MPI_DOUBLE_COMPLEX, &
        &                  MPI_SUM, &
        &                  MPI_COMM_WORLD, &
        &                  ierr)
       call MPI_AllReduce(MPI_IN_PLACE, &
       &                  sigma, &
       &                  ncol, &  
       &                  MPI_DOUBLE_COMPLEX, &
       &                  MPI_SUM, &
       &                  MPI_COMM_WORLD, &
       &                  ierr)
       call MPI_AllReduce(MPI_IN_PLACE, &
        &                  seebeck, &
        &                  ncol, &  
        &                  MPI_DOUBLE_COMPLEX, &
        &                  MPI_SUM, &
        &                  MPI_COMM_WORLD, &
        &                  ierr)
        call MPI_AllReduce(MPI_IN_PLACE, &
       &                  thermalcond, &
        &                  ncol, &  
        &                  MPI_DOUBLE_COMPLEX, &
        &                  MPI_SUM, &
        &                  MPI_COMM_WORLD, &
        &                  ierr)
        call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
        
        if (rank==0) then

           !n=0
           do imu = 1, size(muarray)
              mu = muarray(imu)
              
              do itemp = 1, size(temparray)
                 temp = temparray(itemp)

                 !n=n+1
                 if (ttdf) then
                    do iw = 1, size(wtdf)
                       if (tdoping) then
                          t1 = -(chelevel(itemp)-e_fermi_shift-wtdf(iw))/(temp*kboltz)
                       else
                          t1 = -(mu-wtdf(iw))/(temp*kboltz)
                       end if
                       t2 = exp(t1)
                       if (abs(t1) > critical_factor) then
                          t3=0
                       else
                          t2 = exp(t1)
                          t3 = t2/((t2+1)**2 * temp*kboltz)
                       end if
                       sigma(imu, itemp) = sigma(imu, itemp) + td(iw)*t3
                       seebeck(imu, itemp) = seebeck(imu, itemp) + td(iw)*t3*t1
                       thermalcond(imu, itemp) = thermalcond(imu, itemp) + td(iw)*t3*(t1**2)
                    end do ! iw
                    t1 = tdf_spacing !(wtdff-wtdfi)/dble(size(wtdf)-1)
                 else
                    t1 = 1.0
                 end if
                 !zt1 = 1.0/(1.0*omega)
                 zt1 = 1/(omega*dble(nkptnr))
                 sigma(imu, itemp) = t1*zt1*sigma(imu, itemp) 
                 seebeck(imu, itemp) = (-1)*t1*zt1*kboltz*seebeck(imu, itemp)/(sigma(imu, itemp))
                 thermalcond(imu, itemp) = t1*zt1*(kboltz**2)*temp*thermalcond(imu, itemp)
                 zt(imu, itemp) = (seebeck(imu, itemp)**2*sigma(imu, itemp))/(thermalcond(imu, itemp))

              end do ! temp
           end do ! mu

           !conductivity(:) = zt1*conductivity(:)

     
        end if


!-----------------------------------------------
! Print the output spectra
!-----------------------------------------------
      if (rank==0) then
! output energy units
            t1 = 1.0d0
            if (input%properties%boltzequ%evOutputEnergies) t1 = hartree_to_ev

            if (ttdf) then 
! write the transport distribution function from Boltzmann
            !zt1 = 1/(omega*dble(nkptnr))
               zt1 = 1.0/(omega*1.0*dble(nkptnr))
               write(fname, '("TD_", 2I1, ".OUT")') a, b
               write(*, '("  Transport distribution of component ", 2I1, " written to ", a)') a, b, trim(adjustl(fname))
               open(60, file=trim(fname), action='WRITE', form='FORMATTED')
               t3 = 1.0d0
               if (input%properties%boltzequ%siOutputUnits) t3 = hartree_to_j/(hbar_si**2*bohr_radius_si)
               do iw = 1, size(wtdf)
                  write(60, '(3G18.10)') t1*wtdf(iw), td(iw)*zt1*t3
               end do
               close(60)
            end if

! write the optical conductivity from Boltzmann
            write(fname, '("ELECTCOND_", 2I1, ".OUT")') a, b
            write(*, '("  Electrical conductivity over tau (sigma/tau) of component ", 2I1, " written to ", a)') a, b, trim(adjustl(fname))
            open(60, file=trim(fname), action='WRITE', form='FORMATTED')
            t3 = 1.0d0
            if (input%properties%boltzequ%siOutputUnits) t3 = hartree_to_j*(elec_charge)**2/((hbar_si)**2*bohr_radius_si)

            do imu = 1, size(muarray)
               do itemp = 1, size(temparray)
                  if (tdoping) then
                     write(60, '(4G18.10)') temparray(itemp), (chelevel(itemp)-e_fermi_shift)*t1, sigma(imu, itemp)*t3
                  else
                     write(60, '(4G18.10)') temparray(itemp), muarray(imu)*t1, sigma(imu, itemp)*t3
                  end if
               end do
            end do
           close(60)

! write the Seebeck coefficient from Boltzmann
            write(fname, '("SEEBECK_", 2I1, ".OUT")') a, b
            write(*, '("  Seebeck coefficient (S) of component ", 2I1, " written to ", a)') a, b, trim(adjustl(fname))
            open(60, file=trim(fname), action='WRITE', form='FORMATTED')
            t3 = 1.0d0
            if (input%properties%boltzequ%siOutputUnits) t3 = hartree_to_ev !/e!hartree/(elec_charge)

            do imu = 1, size(muarray)
               do itemp = 1, size(temparray)
                  if (tdoping) then
                     write(60, '(4G18.10)') temparray(itemp), (chelevel(itemp)-e_fermi_shift)*t1, seebeck(imu, itemp)*t3
                  else
                     write(60, '(4G18.10)') temparray(itemp), muarray(imu)*t1, seebeck(imu, itemp)*t3
                  end if
               end do
            end do
           close(60)
! write the thermal conductivity from Boltzmann
            write(fname, '("THERMALCOND_", 2I1, ".OUT")') a, b
            write(*, '("  Thermal conductivity over tau (kappa/tau) of component ", 2I1, " written to ", a)') a, b, trim(adjustl(fname))
            open(60, file=trim(fname), action='WRITE', form='FORMATTED')
            t3 = 1.0d0
            if (input%properties%boltzequ%siOutputUnits) t3 = (hartree_to_j)**3/(bohr_radius_si*(hbar_si)**2)

            do imu = 1, size(muarray)
               do itemp = 1, size(temparray)
                  if (tdoping) then
                     write(60, '(4G18.10)') temparray(itemp), (chelevel(itemp)-e_fermi_shift)*t1, thermalcond(imu, itemp)*t3
                  else
                     write(60, '(4G18.10)') temparray(itemp), muarray(imu)*t1, thermalcond(imu, itemp)*t3
                  end if
               end do
            end do
            close(60)            

! write the figure of merit from Boltzmann
            write(fname, '("Z_", 2I1, ".OUT")') a, b
            write(*, '("  Figure of merit Z of component ", 2I1, " written to ", a)') a, b, trim(adjustl(fname))
            open(60, file=trim(fname), action='WRITE', form='FORMATTED')

            do imu = 1, size(muarray)
               do itemp = 1, size(temparray)
                  if (tdoping) then
                     write(60, '(4G18.10)') temparray(itemp), (chelevel(itemp)-e_fermi_shift)*t1, zt(imu, itemp)
                  else
                     write(60, '(4G18.10)') temparray(itemp), muarray(imu)*t1, zt(imu, itemp)
                  end if
               end do
            end do
            close(60)
                      
        end if
        write(*,*)
#ifdef MPI        
        call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
     end do ! l, optical components
     
     if (input%properties%boltzequ%siOutputUnits) then 
        write(*, '("  Output of energy in eV")')
        write(*, '("  Output of transport coefficients in SI units: &
              & S in V/K , sigma/tau in S/(m s) , kappa/tau in W/(mK s) ")')
     endif     
     write(*,*)
        
end subroutine
