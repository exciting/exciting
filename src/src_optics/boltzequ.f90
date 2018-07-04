!
subroutine boltzequ
    use modmain
    use modinput
    use modmpi
    use modxs, only: symt2
    implicit none
! local variables
    integer :: l, a, b, ncomp, condcomp(2,9)
    integer :: ik, jk, isym
    integer :: ist, jst, iw
    integer :: recl, iostat
    integer :: mugrid, tempgrid, nwtdf
    integer :: n
    logical :: exist
    real(8) :: swidth, t1, t2, t3
    real(8) :: v1 (3), v2 (3), sc(3,3)
    real(8) :: kb, ab, hbar, ech, hartree
    real(8) :: wtdfi, wtdff
    real(8) :: tempi,tempf,temps,temp
    real(8) :: mui, muf, mus, mu
    complex(8) :: zt1, zt2, zt3, ncharge
    character(256) :: fname
! allocatable arrays
    real(8), allocatable :: wtdf(:)
    real(8), allocatable :: evalsvt(:), occsvt(:)
    complex (8), allocatable :: pmat(:,:,:) 
    complex (8), allocatable :: sigmab(:) ! conductivity
    complex (8), allocatable :: seebeck(:) ! Seebeck coefficient
    complex (8), allocatable :: thermalcond(:) ! thermal conductivity
    complex (8), allocatable :: conductivity(:) ! Conductivity avoiding the 
    complex (8), allocatable :: zt(:)          ! ZT
    complex (8), allocatable :: echarge(:) ! charge density
    complex (8), allocatable :: td(:) ! Transport distribution function
    complex (8), allocatable :: ndos(:) ! density of states
    integer :: kfirst, klast
    
! external functions
    real(8) sdelta

    kb =  3.1668114d-6 ! Boltzmann constant Ha/K
    
    ab = 5.2917721092d-11
    hbar = 1.054571726d-34
    ech = 1.602176565d-19
    hartree = 4.35974417d-18
    
! initialise universal variables
    call init0
    call init1

! read Fermi energy from file
    call readfermi

! working arrays
    mugrid = input%properties%boltzequ%mugrid
    tempgrid = input%properties%boltzequ%tgrid
    nwtdf = input%properties%boltzequ%nwtdf
    
    allocate(wtdf(nwtdf))
    allocate(td(nwtdf))
    allocate(ndos(nwtdf))
    allocate(sigmab(tempgrid*mugrid))
    allocate(seebeck(tempgrid*mugrid))
    allocate(thermalcond(tempgrid*mugrid))
    allocate(conductivity(tempgrid*mugrid))
    allocate(zt(tempgrid*mugrid))
    allocate(echarge(tempgrid*mugrid))

    wtdfi = input%properties%boltzequ%windtdf(1)
    wtdff = input%properties%boltzequ%windtdf(2)
    t1 = (wtdff-wtdfi)/dble(nwtdf-1)       
    Do iw = 1, nwtdf
       wtdf(iw) = wtdfi + t1*(iw-1)
    End Do

! smearing factor
    swidth = input%properties%boltzequ%swidth

! the momentum matrix elements
    allocate(pmat(3,nstsv,nstsv))
    inquire(iolength=recl) pmat
    open(50,File='PMAT.OUT',Action='READ',Form='UNFORMATTED', &
    &    Access='DIRECT',Recl=recl,IOstat=iostat)

! KS eigenvalues and occupancies    
    allocate(evalsvt(nstsv),occsvt(nstsv))

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

    ncomp = size(input%properties%boltzequ%condcomp,2)
    if (ncomp>0) then
        do l = 1, ncomp
            condcomp(1,l) = input%properties%boltzequ%condcomp(1,l)
            condcomp(2,l) = input%properties%boltzequ%condcomp(2,l)
        end do
    else

        ncomp = 0
        do a = 1, 3
            do b = 1, 3
                if (sum(abs(symt2(a,b,:,:)))>1.0d-8) then
                    ncomp = ncomp+1
                    condcomp(1,ncomp) = a
                    condcomp(2,ncomp) = b
                end if
            end do
        end do
    end if

! initialize 

    tempi=input%properties%boltzequ%windtemp(1)
    tempf=input%properties%boltzequ%windtemp(2)
    if (tempi == tempf) then
       tempf=tempi+1
       temps=temps+2
    else if (input%properties%boltzequ%tgrid == 1) then
       temps=tempf+1
    else
       temps=(tempf-tempi)/(input%properties%boltzequ%tgrid-1)
    end if

    mui=input%properties%boltzequ%windmu(1)
    muf=input%properties%boltzequ%windmu(2)
    if (mui == muf) then
       muf=mui+1
       mus=mus+2
    else if (input%properties%boltzequ%mugrid == 1) then
       mus=muf+1
    else
       mus=(muf-mui)/(input%properties%boltzequ%mugrid-1)
    end if

    
    do l = 1, ncomp
        
        a = condcomp(1,l)
        b = condcomp(2,l)

! sum over non-reduced k-points
        conductivity(:) = zzero
        sigmab(:) = zzero
        seebeck(:) = zzero
        thermalcond(:) = zzero
        zt(:) = zzero
        echarge(:) = zzero
        td(:) = zzero
        ndos(:) = zzero

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
            call getoccsv(vkl(:,jk),occsvt)

!--------------------------
! Transport distribution
!--------------------------       
            do ist = 1, nstsv
               zt1 = occmax*pmat(a,ist,ist)*conjg(pmat(b,ist,ist))
               do iw = 1, nwtdf
                  t1 = (evalsvt(ist)-wtdf(iw))/swidth
                  td(iw) = td(iw)+ zt1*sdelta(input%groundstate%stypenumber,t1)/swidth
                  ndos(iw) = ndos(iw) + occmax*sdelta(input%groundstate%stypenumber,t1)/swidth
               end do ! iw
            end do ! ist


            n=0
            do mu = mui, muf, mus
           
               do temp = tempi, tempf, temps

                  n=n+1
                  do ist = 1, nstsv
                  !do iw = 1, nwtdf
                     t1 = -(mu-evalsvt(ist))/(temp*kb)
                     t2 = exp(t1)
                     if (t1 > 40) then
                        t3=0
                     else if (t1 < -40 ) then
                        t3=0
                     else
                        t2 = exp(t1)
                        t3 = t2/((t2+1)**2 * temp*kb)
                     end if
                     zt1=occmax*pmat(a,ist,ist)*conjg(pmat(b,ist,ist))
                     conductivity(n) = conductivity(n) + zt1*t3
                  end do ! iw
               end do ! mu
            end do ! temp

         end do
                    
#ifdef MPI
        call MPI_AllReduce(MPI_IN_PLACE, &
        &                  td, &
        &                  nwtdf, &  
        &                  MPI_DOUBLE_COMPLEX, &
        &                  MPI_SUM, &
        &                  MPI_COMM_WORLD, &
        &                  ierr)
        call MPI_AllReduce(MPI_IN_PLACE, &
        &                  conductivity, &
        &                  nwtdf, &  
        &                  MPI_DOUBLE_COMPLEX, &
        &                  MPI_SUM, &
        &                  MPI_COMM_WORLD, &
        &                  ierr)
        call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
        
        if (rank==0) then

           ncharge=0.0d0
           do iw = 1, nwtdf
              if ( wtdf(iw) < efermi ) then
                 ncharge=ncharge+ndos(iw)
              else if ( wtdf(iw) == efermi) then
                 ncharge=ncharge+0.5*ndos(iw)
              end if
           end do ! iw

           n=0
           do mu = mui, muf, mus
           
              do temp = tempi, tempf, temps

                 n=n+1
                 do iw = 1, nwtdf
                    t1 = -(mu-wtdf(iw))/(temp*kb)
                    t2 = exp(t1)
                    if (t1 > 40) then
                       t3=0
                    else if (t1 < -40 ) then
                       t3=0
                    else
                       t2 = exp(t1)
                       t3 = t2/((t2+1)**2 * temp*kb)
                    end if
                    sigmab(n) = sigmab(n) + td(iw)*t3
                    seebeck(n) = seebeck(n) + td(iw)*t3*t1
                    thermalcond(n) = thermalcond(n) + td(iw)*t3*(t1**2)
                    echarge(n)=echarge(n) + ndos(iw)/(exp(t1)+1)
              
                 end do ! iw
                 zt1 = 1/(omega*dble(nkptnr))
                 t1 = (wtdff-wtdfi)/dble(nwtdf-1)
                 write(*,*) t1
                 sigmab(n) = t1*zt1*sigmab(n) 
                 seebeck(n) = (-1)*t1*zt1*kb*seebeck(n)/(sigmab(n))
                 thermalcond(n) = t1*zt1*(kb**2)*temp*thermalcond(n)
                 zt(n) = (seebeck(n)**2*sigmab(n))/(thermalcond(n))*temp
                 echarge(n) = zt1*t1*(ncharge-echarge(n))

                 conductivity(n) = zt1*conductivity(n)

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
            if (input%properties%boltzequ%tevout) t1 = h2ev

! write the transport distribution function from Boltzmann
            zt1 = 1/(omega*dble(nkptnr))
            write(fname, '("TD_", 2I1, ".OUT")') a, b
            write(*, '("  Transport distribution  of component ", 2I1, " written to ", a)') a, b, trim(adjustl(fname))
            open(60, file=trim(fname), action='WRITE', form='FORMATTED')
            t3 = 1.0d0
            if (input%properties%boltzequ%tsiout) t3 = 1/(hbar*ab)
            do iw = 1, nwtdf
                write(60, '(6G18.10)') t1*wtdf(iw), td(iw)*zt1*t3, ndos(iw)*zt1
            end do
           close(60)

! write the optical conductivity from Boltzmann
            write(fname, '("ELECTCOND_", 2I1, ".OUT")') a, b
            write(*, '("  Electrical conductivity of component ", 2I1, " written to ", a)') a, b, trim(adjustl(fname))
            open(60, file=trim(fname), action='WRITE', form='FORMATTED')
            t3 = 1.0d0
            if (input%properties%boltzequ%tsiout) t3 = hartree*(ech)**2/((hbar)**2*ab)

            n=0
            do mu = mui, muf, mus
               do temp = tempi, tempf, temps
                  n=n+1
                  write(60, '(4G18.10)') temp, mu*t1, sigmab(n)*t3
               end do
               !write(60, ' ')
            end do
           close(60)

! write the Seebeck coefficient from Boltzmann
            write(fname, '("SEEBECK_", 2I1, ".OUT")') a, b
            write(*, '("  Seebeck coefficient of component ", 2I1, " written to ", a)') a, b, trim(adjustl(fname))
            open(60, file=trim(fname), action='WRITE', form='FORMATTED')
            t3 = 1.0d0
            if (input%properties%boltzequ%tsiout) t3 = h2ev !/e!hartree/(ech)

            n=0
            do mu = mui, muf, mus
               do temp = tempi, tempf, temps
                  n=n+1
                  write(60, '(4G18.10)') temp, mu*t1, seebeck(n)*t3
               end do
            end do
           close(60)
! write the thermal conductivity from Boltzmann
            write(fname, '("THERMALCOND_", 2I1, ".OUT")') a, b
            write(*, '("  Thermal conductivity of component ", 2I1, " written to ", a)') a, b, trim(adjustl(fname))
            open(60, file=trim(fname), action='WRITE', form='FORMATTED')
            t3 = 1.0d0
            if (input%properties%boltzequ%tsiout) t3 = (hartree)**3/(ab*(hbar)**2)

            n=0
            do mu = mui, muf, mus
               do temp = tempi, tempf, temps
                  n=n+1
                  write(60, '(4G18.10)') temp, mu*t1, thermalcond(n)*t3
               end do
            end do
            close(60)            

! write the figure of merit from Boltzmann
            write(fname, '("ZT_", 2I1, ".OUT")') a, b
            write(*, '("  Figure of merit zT of component ", 2I1, " written to ", a)') a, b, trim(adjustl(fname))
            open(60, file=trim(fname), action='WRITE', form='FORMATTED')
            n=0
            do mu = mui, muf, mus
               do temp = tempi, tempf, temps
                  n=n+1
                  write(60, '(4G18.10)') temp, mu*t1, zt(n)
               end do
            end do
            close(60)

! write the electrical conductivity from Boltzmann
            !write(fname, '("SIGMA_B_", 2I1, ".OUT")') a, b
            !write(*, '("  Alternative calculation of the electrical conductivity for component ", 2I1, "written to ", a)') a, b, trim(adjustl(fname))
            !open(60, file=trim(fname), action='WRITE', form='FORMATTED')
            !t3 = 1.0d0
            !if (input%properties%boltzequ%tsiout) t3 = hartree*(ech)**2/((hbar)**2*ab)

            !n=0
            !do mu = mui, muf, mus
            !   do temp = tempi, tempf, temps
            !      n=n+1
            !      write(60, '(4G18.10)') temp, mu*t1, conductivity(n)*t3
            !   end do
               !write(60, ' ')
            !end do
           !close(60)

! write the total charge for Boltzmann
            write(fname, '("CHARGE_", 2I1, ".OUT")') a, b
            write(*, '("  Total charge for a given chemical potential and temperature ", a)') trim(adjustl(fname))
            open(60, file=trim(fname), action='WRITE', form='FORMATTED')
            t3 = 1.0d0
            if (input%properties%boltzequ%tsiout) t3 = 1/((ab)**3)*10.0d-6

            n=0
            do mu = mui, muf, mus
               do temp = tempi, tempf, temps
                  n=n+1
                  write(60, '(6G18.10)') temp, mu*t1, echarge(n), echarge(n)*t3
               end do
            end do
           close(60)
                      
            if (input%properties%boltzequ%tevout) write(*, '("  Output of energy is in eV")')
            if (input%properties%boltzequ%tsiout) write(*, '("  Output of transport coefficients is in SI")')
            write(*,*)
        end if
#ifdef MPI        
        call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
    end do ! l, optical components
        
    deallocate(pmat)
    deallocate(evalsvt,occsvt)
    deallocate(sigmab,seebeck,thermalcond,zt,td,ndos,echarge,conductivity)
    
end subroutine
