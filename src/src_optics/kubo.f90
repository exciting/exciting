!
subroutine kubo
    use modmain
    use modinput
    use modmpi
    use modxs, only: symt2
    implicit none
! local variables
    integer :: l, a, b, ncomp, epscomp(2,9)
    integer :: ik, jk, isym
    integer :: ist, jst, iw
    integer :: recl, iostat
    integer :: mugrid, tempgrid
    integer :: n
    logical :: exist, intraband, drude
    real(8) :: wplas, swidth, eji, t1, t2, t3
    real(8) :: v1 (3), v2 (3), sc(3,3)
    real(8) :: kb
    real(8) :: deltaw, dshift
    real(8) :: tempi,tempf,temps,temp
    real(8) :: mui, muf, mus, mu, deltamu
    complex(8) :: sigmaboltz(9)
    complex(8) :: sigmak(9)
    complex(8) :: sigmaktwo(9)
    complex(8) :: zt1, zt2, zt3
    character(256) :: fname
! allocatable arrays
    real(8), allocatable :: w(:)
    real(8), allocatable :: evalsvt(:), occsvt(:)
    complex (8), allocatable :: pmat(:,:,:)
    complex (8), allocatable :: pmatd(:,:,:)
    complex (8), allocatable :: sigma(:)
    complex (8), allocatable :: sigmab(:)
    complex (8), allocatable :: seebeck(:)
    complex (8), allocatable :: thermalcond(:)
    complex (8), allocatable :: zt(:)
    complex (8), allocatable :: echarge(:)
!    complex (8), allocatable :: nf(:)
    complex (8), allocatable :: td(:) ! Transport distribution function
    complex (8), allocatable :: ndos(:) ! Density of states
    complex (8), allocatable :: eta(:)
    integer :: wgrid
    integer :: kfirst, klast
    real(8) :: wmax, scissor
! external functions
    real(8) sdelta

    kb =  0.30036326d-5 ! Boltzmann constant Ha/K

    
    if (rank==0) then
      write(*,*)
      write(*,'("Info(kubo):")')
    end if
    
! initialise universal variables
    call init0
    call init1

! read Fermi energy from file
    call readfermi

! working arrays
    wgrid = input%properties%dielmat%wgrid
    mugrid = input%properties%dielmat%mugrid
    tempgrid = input%properties%dielmat%tempgrid
    dshift = input%properties%dielmat%dshift
!    write(*,*) wgrid
!    write(*,*) mugrid
!    write(*,*) tempgrid
    
    allocate(w(wgrid))
    allocate(td(wgrid))
    allocate(ndos(wgrid))
    allocate(sigma(wgrid))
    allocate(sigmab(tempgrid*mugrid))
    allocate(seebeck(tempgrid*mugrid))
    allocate(thermalcond(tempgrid*mugrid))
    allocate(zt(tempgrid*mugrid))
    allocate(echarge(tempgrid*mugrid))
    !allocate(nf(tempgrid*mugrid))

! generate energy grid
    !wmax = input%properties%dielmat%wmax
    deltaw = input%properties%dielmat%deltaw
    
    t1 = deltaw/(dble(wgrid))
    Do iw = 1, wgrid
        w(iw) = t1*(iw-1)-deltaw/(1.0*2)+efermi
    End Do

! if element "drude" present, use the specified parameters for the Drude term 
! instead of calculating the plasma frequency
    drude = .false.
    if (input%properties%dielmat%intraband) then
        ! if omega_p is different from default (negative) value
        if (input%properties%dielmat%drude(1)>1.d-8) then
            drude = .true.
            if (rank==0) then
                write(*,*) 
                write(*,*) '  Drude model parameters:'
                write(*,*) '     plasma frequency   : ', input%properties%dielmat%drude(1)
                write(*,*) '     lifetime broadening: ', input%properties%dielmat%drude(2)
            end if
        end if
    end if

! input%properties%dielmat%scissor
    scissor = input%properties%dielmat%scissor

! smearing factor
    swidth = input%properties%dielmat%swidth

! finite lifetime broadening to mimic the experimental resolution
! chosen according to PRB86, 125139 (2012)
    allocate(eta(wgrid))
    do iw = 1, wgrid
        eta(iw) = cmplx(0.d0,swidth)
        !eta(iw) = cmplx(0.d0,swidth+0.1*swidth*w(iw))
    end do

! calculate momentum matrix elements
    !if (rank==0) then
    !    write(*,*)
    !    write(*,'("  Calculate the momentum matrix elements")')
    !end if
    !call writepmat

! the momentum matrix elements
    allocate(pmat(3,nstsv,nstsv))
    inquire(iolength=recl) pmat
    open(50,File='PMAT.OUT',Action='READ',Form='UNFORMATTED', &
    &    Access='DIRECT',Recl=recl,IOstat=iostat)

! the momentum matrix elements with delta for second derivative
    if (dshift > 0) then
       allocate(pmatd(3,nstsv,nstsv))
       inquire(iolength=recl) pmatd
       open(80,File='PMATD.OUT',Action='READ',Form='UNFORMATTED', &
            &    Access='DIRECT',Recl=recl,IOstat=iostat)
    end if

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

    ncomp = size(input%properties%dielmat%epscomp,2)
    if (ncomp>0) then
        do l = 1, ncomp
            epscomp(1,l) = input%properties%dielmat%epscomp(1,l)
            epscomp(2,l) = input%properties%dielmat%epscomp(2,l)
        end do
    else
!
! calculate \eps_{ij}} for all symmetry inequivalent components,
! by analysing the symmetrization matrices
! 
        ncomp = 0
        do a = 1, 3
            do b = 1, 3
                if (sum(abs(symt2(a,b,:,:)))>1.0d-8) then
                    ncomp = ncomp+1
                    epscomp(1,ncomp) = a
                    epscomp(2,ncomp) = b
                end if
            end do
        end do
    end if

! initialize sigmaboltz to zeros
    sigmaboltz(:) = zzero
    sigmak(:) = zzero
    sigmaktwo(:) = zzero

    tempi=input%properties%dielmat%temprange(1)
    tempf=input%properties%dielmat%temprange(2)
    if (tempi == tempf) then
       tempf=tempi+1
       temps=temps+2
    else if (input%properties%dielmat%tempgrid == 1) then
       temps=tempf+1
    else
       temps=(tempf-tempi)/(input%properties%dielmat%tempgrid-1)
    end if
    
    deltamu=input%properties%dielmat%deltamu
    if (input%properties%dielmat%mugrid == 1) then
       mus=1
    else
       mus=(deltamu)/(input%properties%dielmat%mugrid-1)
    end if

    if (mus == 1) then
       mui=efermi
       muf=efermi+0.5
    else
       mui=efermi-deltamu/2.0
       muf=efermi+deltamu/2.0
    end if

    do l = 1, ncomp
        
        a = epscomp(1,l)
        b = epscomp(2,l)

! flag for calculating the intraband term
        intraband = input%properties%dielmat%intraband .and. (a==b)

! sum over non-reduced k-points
        sigma(:) = zzero
        sigmab(:) = zzero
        seebeck(:) = zzero
        thermalcond(:) = zzero
        zt(:) = zzero
        echarge(:) = zzero
        !nf(:) = zzero
        td(:) = zzero
        ndos(:) = zzero
        wplas = 0.0d0
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

            if ( dshift>0 ) then
! read momentum matrix elements of shifted k grid by delta from direct-access file
               read(80,Rec=jk) pmatd

! rotate the matrix elements from the reduced to non-reduced k-point
               if (isym > 1) then
                  do ist = 1, nstsv
                     do jst = 1, nstsv
                        sc(:,:)=symlatc(:,:,lsplsymc(isym))
                        call r3mv(sc,dble(pmatd(:,ist,jst)),v1)
                        call r3mv(sc,aimag(pmatd(:,ist,jst)),v2)
                        pmatd(:,ist,jst) = cmplx(v1(:),v2(:),8)
                     end do
                  end do
               end if
            end if

            
! read eigenvalues and occupancies from files
            call getevalsv(vkl(:,jk),evalsvt)
            call getoccsv(vkl(:,jk),occsvt)

!--------------------------
! PLASMA FREQUENCY
!--------------------------
            if (intraband.and.(.not.drude)) then
                do ist = 1, nstsv
                    zt1 = occsvt(ist)*pmat(a,ist,ist)*conjg(pmat(b,ist,ist))
                    t1 = (evalsvt(ist)-efermi)/input%groundstate%swidth
                    wplas = wplas+dble(zt1)*sdelta(input%groundstate%stypenumber,t1)/input%groundstate%swidth
                end do
            end if

!--------------------------
! Transport distribution
!--------------------------       
!--------------------------       
            do ist = 1, nstsv
               zt1 = pmat(a,ist,ist)*conjg(pmat(b,ist,ist))
               do iw = 1, wgrid
                  t1 = (evalsvt(ist)-w(iw))/input%groundstate%swidth
                  !td(iw) = td(iw)+ zt1/input%properties%dielmat%drude(2)*sdelta(input%groundstate%stypenumber,t1)/input%groundstate%swidth
                  !write(*,*) 1/input%properties%dielmat%drude(2)
                  td(iw) = td(iw)+ zt1*sdelta(input%groundstate%stypenumber,t1)/input%groundstate%swidth
                  ndos(iw) = ndos(iw) + sdelta(input%groundstate%stypenumber,t1)/input%groundstate%swidth
               end do ! iw
            end do ! ist

!--------------------------
! sigmak from Kubo
!--------------------------       
!--------------------------       
! formula with first derivative: produkt of pmat elements:
!-------------------------       
            do ist = 1, nstsv
               !zt1 = pmat(a,ist,ist)*conjg(pmat(b,ist,ist))
               !t1 = (evalsvt(ist)-efermi)/(tempi*kb)
               !sigmak(l) = sigmak(l)+ zt1*sdelta(3,t1)/(tempi*kb)
               !if (occsvt(ist) > 0 ) then
               zt1 = pmat(a,ist,ist)*conjg(pmat(b,ist,ist))
               t1 = (evalsvt(ist)-efermi)/(tempi*kb)
               if (t1 > 20) then
                  t3=0
               else if (t1 < -20 ) then
                  t3=0
               else
                  t2 = exp(t1)
                  t3 = t2/((t2+1)**2 * tempi*kb)
               end if
               !write(*,*) occsvt(ist), evalsvt(ist), efermi
               !write(*,*) zt1, t3
               sigmak(l) = sigmak(l)+ zt1*t3
               !end if
            end do ! ist

!--------------------------       
! formula with second derivative: determined by pmat and pmatd
!--------------------------
            if ( dshift > 0) then
               do ist = 1, nstsv
                  !zt1 = pmat(a,ist,ist)*conjg(pmat(b,ist,ist))
                  !t1 = (evalsvt(ist)-efermi)/(tempi*kb)
                  !sigmak(l) = sigmak(l)+ zt1*sdelta(3,t1)/(tempi*kb)
                  !if (occsvt(ist) > 0 ) then
                  zt3 = (pmatd(a,ist,ist)-pmat(b,ist,ist))/dshift
                  t1 = (evalsvt(ist)-efermi)/(tempi*kb)
                  !if (t1 > 20) then
                  !   t3=0
                  !else if (t1 < -20 ) then
                  !   t3=0
                  !else
                  t2 = exp(t1)
                  t3 = 1/(t2+1)
                  !end if
                  !write(*,*) occsvt(ist), evalsvt(ist), efermi
                  !write(*,*) zt1, t3
                  sigmaktwo(l) = sigmaktwo(l)+ zt1*t3
                  !end if
               end do ! ist
            end if

            
!--------------------------
! INTERBAND CONTRIBUTION
!--------------------------       
            do ist = 1, nstsv
            if (evalsvt(ist)<=efermi) then
            
                do jst = 1, nstsv
                if (evalsvt(jst)>efermi) then
            
                    zt1 = pmat(a,ist,jst)*conjg(pmat(b,ist,jst))
                    eji = evalsvt(jst)-evalsvt(ist)
                
                    ! scissor operator
                    if (dabs(scissor)>1.0d-8) then
                        if (dabs(eji)>1.0d-8) then
                            t1 = (eji+scissor)/eji
                            zt1 = zt1*t1*t1
                            eji = eji+scissor
                        end if
                    end if
                
                    if (dabs(eji)>1.0d-8) then
                        t1=occsvt(ist)*(1.0d0-occsvt(jst)/occmax)/eji
                        do iw = 1, wgrid
                            sigma(iw) = sigma(iw)+ &
                           &  t1*(zt1/(w(iw)-eji+eta(iw)) + &
                           &  conjg(zt1)/(w(iw)+eji+eta(iw)))
                        end do ! iw
                    end if

                end if
                end do ! jst
        
            end if
            end do ! ist
            
        end do ! ik


        !do iw = 1, wgrid
        !   t1 = (efermi-w(iw))/(input%properties%dielmat%t*kb)
        !   sigmab(iw) = sigmab(iw)+ td(iw)*sdelta(3,t1)/(input%properties%dielmat%t*kb)
           !sigmaboltz(l) = sigmaboltz(l) + td(iw)*sdelta(3,t1)/(input%properties%dielmat%t*kb)
           !write (*,*) td(iw), sdelta(3,t1), t1
        !end do ! iw
        !zt1 = twopi**3/(omega*dble(nkptnr))
        !t1 = deltaw/dble(wgrid)
        !sigmab(:) = t1*zt1*sigmab(:)
               
        


        !sigmab(:) = zt1*sigmab(:)
        !sigmaboltz(l) = zt1*sigmaboltz(l)
        
#ifdef MPI
        call MPI_AllReduce(MPI_IN_PLACE, &
        &                  sigma, &
        &                  wgrid, &
        &                  MPI_DOUBLE_COMPLEX, &
        &                  MPI_SUM, &
        &                  MPI_COMM_WORLD, &
        &                  ierr)
        call MPI_AllReduce(MPI_IN_PLACE, &
        &                  td, &
        &                  wgrid, &
        &                  MPI_DOUBLE_COMPLEX, &
        &                  MPI_SUM, &
        &                  MPI_COMM_WORLD, &
        &                  ierr)
        if (intraband.and.(.not.drude)) then
          call MPI_AllReduce(MPI_IN_PLACE, &
          &                  wplas, &
          &                  1, &
          &                  MPI_DOUBLE_PRECISION, &
          &                  MPI_SUM, &
          &                  MPI_COMM_WORLD, &
          &                  ierr)
        end if
        call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
        if (rank==0) then
        n=0
        do mu = mui, muf, mus
           
           do temp = tempi, tempf, temps

        !do temp = tempi, tempf, temps
           !do mu = mui, muf, mus
              n=n+1
              !write(*,*) n
              do iw = 1, wgrid
                 t1 = -(mu-w(iw))/(temp*kb)
                 t2 = exp(t1)
                 t3 = t2/((t2+1)**2 * temp*kb)
                 sigmab(n) = sigmab(n) + td(iw)*t3
                 seebeck(n) = seebeck(n) + td(iw)*t3*t1
                 thermalcond(n) = thermalcond(n) + td(iw)*t3*(t1**2)
                 echarge(n)=echarge(n) + ndos(iw)/(exp(t1)+1)
                 !nf(n)=nf(n)+ndos(iw)
                 !sigmab(n) = sigmab(n) + td(iw)*sdelta(3,t1)/(temp*kb)
                 !seebeck(n) = seebeck(n) + td(iw)*sdelta(3,t1)/(temp*kb)*(-1)*t1
                 !thermalcond(n) = thermalcond(n) + td(iw)*sdelta(3,t1)/(temp*kb)*(t1**2)
              !write (*,*) td(iw), sdelta(3,t1), t1, sigmaboltz(l), l
              end do ! iw
              zt1 = twopi**3/(omega*dble(nkptnr))
              t1 = deltaw/dble(wgrid)
              sigmab(n) = t1*zt1*sigmab(n) !/input%properties%dielmat%drude(2)
              seebeck(n) = t1*zt1*kb*temp*seebeck(n)/sigmab(n)
              thermalcond(n) = t1*zt1*(kb**2)*temp*temp*thermalcond(n)
              zt(n) = sigmab(n)*seebeck(n)**2/(thermalcond(n))
              echarge(n) = t1*zt1*echarge(n)
              !nf(n)=t1*zt1*nf(n)*0.5
              !write(*,*) nf(n)
              !zt(n) = t1*zt1*(kb**2)*temp*temp*thermalcond(n)
              !write(*,*) temp, mu, sigmab(n)

           end do ! mu
        end do ! temp

        zt1 = zi/(omega*dble(nkptnr))
        sigma(:) = zt1*sigma(:)
        write (*,*) tempi, kb, efermi
        do iw = 1, wgrid
           !t1 = (efermi-w(iw))/(tempi*kb)
           !sigmaboltz(l) = sigmaboltz(l) + td(iw)*sdelta(3,t1)/(tempi*kb)
           t1 = -(efermi-w(iw))/(tempi*kb)
           t2 = exp(t1)
           t3 = t2/((t2+1)**2 * tempi*kb)
           sigmaboltz(l) = sigmaboltz(l) + td(iw)*t3
           !write (*,*) td(iw), sdelta(3,t1), t1, sigmaboltz(l), l
        end do ! iw
        
        zt1 = twopi**3/(omega*dble(nkptnr))
        t1 = deltaw/dble(wgrid)
        sigmaboltz(l) = t1*zt1*sigmaboltz(l)
        write (*,*) sigmaboltz(l)

        sigmak(l) = zt1*sigmak(l)
        write (*,*) sigmak(l)

        if (dshift > 0 ) then
           sigmaktwo(l) = zt1*sigmaktwo(l)
           write (*,*) sigmaktwo(l)
        end if

        
! add the intraband contribution
        if (intraband) then
            if (drude) then
! use specified parameters for Drude term
                zt1 = zi/fourpi

                do iw = 1, wgrid
                    t1 = input%properties%dielmat%drude(1)
                    t2 = input%properties%dielmat%drude(2)
                    sigma(iw) = sigma(iw)+zt1*t1*t1/(w(iw)+zi*t2)
                end do
            else
! use the calculated plasma frequency            
                zt1 = fourpi/(omega*dble(nkptnr))
                wplas = zt1*dabs(wplas)
                if (rank==0) then
                    write(fname, '("PLASMA_", 2I1, ".OUT")') a, b
                    write(*, '("  Plasma frequency written to ", a)') trim(adjustl(fname))
                    open(60, File=trim(fname), Action='WRITE', Form='FORMATTED')
                    write(60, '(G18.10, " : plasma frequency")') sqrt(wplas)
                    close(60)
                end if
                zt1 = zi/fourpi
                do iw = 1, wgrid
                    sigma(iw) = sigma(iw)+zt1*wplas/(w(iw)+eta(iw))
                end do
            end if
         end if
         end if

        !-----------------------------------------------
! Print the output spectra
!-----------------------------------------------
        if (rank==0) then
! output energy units
            t1 = 1.0d0
            if (input%properties%dielmat%tevout) t1 = h2ev
! write the optical conductivity
            write(fname, '("SIGMA_", 2I1, ".OUT")') a, b
            write(*, '("  Optical conductivity tensor written to ", a)') trim(adjustl(fname))
            write(*,*) trim(fname)
            
            open(61, file=trim(fname), action='WRITE', form='FORMATTED')
            do iw = 1, wgrid
                write(61, '(3G18.10)') t1*w(iw), sigma(iw)
            end do
            close(61)
! write the transport distribution function from Boltzmann
            write(fname, '("TD_", 2I1, ".OUT")') a, b
            write(*, '("  Transport distribution written to ", a)') trim(adjustl(fname))
            open(60, file=trim(fname), action='WRITE', form='FORMATTED')
            do iw = 1, wgrid
                !write(60, '(4G18.10)') t1*w(iw), td(iw)
            end do
           close(60)
! write the optical conductivity from Boltzmann
            write(fname, '("SIGMAB_", 2I1, ".OUT")') a, b
            write(*, '("  Optical conductivity tensor written to ", a)') trim(adjustl(fname))
            open(60, file=trim(fname), action='WRITE', form='FORMATTED')
            n=0
            do mu = mui, muf, mus
               do temp = tempi, tempf, temps
               !do mu = mui, muf, mus
                  n=n+1
                  write(60, '(4G18.10)') temp, mu, sigmab(n)
               end do
               !write(60, ' ')
            end do
           close(60)
! write the Seebeck coefficient from Boltzmann
            write(fname, '("SEEBECK_", 2I1, ".OUT")') a, b
            write(*, '("  Seebeck coefficient tensor written to ", a)') trim(adjustl(fname))
            open(60, file=trim(fname), action='WRITE', form='FORMATTED')
            n=0
            do mu = mui, muf, mus
               do temp = tempi, tempf, temps
               !do mu = mui, muf, mus
                  n=n+1
                  write(60, '(4G18.10)') temp, mu, seebeck(n)
               end do
               !write(60, ' ')
            end do
           close(60)
! write the thermal conductivity from Boltzmann
            write(fname, '("THERMALCOND_", 2I1, ".OUT")') a, b
            write(*, '("  Thermal conductivity tensor written to ", a)') trim(adjustl(fname))
            open(60, file=trim(fname), action='WRITE', form='FORMATTED')
            n=0
            do mu = mui, muf, mus
               do temp = tempi, tempf, temps
               !do mu = mui, muf, mus
                  n=n+1
                  write(60, '(4G18.10)') temp, mu, thermalcond(n)
               end do
               !write(60, ' ')
            end do
           close(60)

! write the figure of merit from Boltzmann
            write(fname, '("ZT_", 2I1, ".OUT")') a, b
            write(*, '("  Figure of merit zT written to ", a)') trim(adjustl(fname))
            open(60, file=trim(fname), action='WRITE', form='FORMATTED')
            n=0
            do mu = mui, muf, mus
               do temp = tempi, tempf, temps
               !do mu = mui, muf, mus
                  n=n+1
                  write(60, '(4G18.10)') temp, mu, zt(n)
               end do
            end do
           close(60)

! write the total charge for Boltzmann
            write(fname, '("CHARGE_", 2I1, ".OUT")') a, b
            write(*, '("  Total charge for a given chemical potential and temperature ", a)') trim(adjustl(fname))
            open(60, file=trim(fname), action='WRITE', form='FORMATTED')
            n=0
            do mu = mui, muf, mus
               do temp = tempi, tempf, temps
               !do mu = mui, muf, mus
                  n=n+1
                  write(60, '(4G18.10)') temp, mu, echarge(n)
               end do
            end do
           close(60)
           
           
! write the dielectric function to file
            write(fname, '("EPSILON_", 2I1, ".OUT")') a, b
            write(*, '("  Dielectric tensor written to ", a)') trim(adjustl(fname))
            open(60, file=trim(fname), action='WRITE', form='FORMATTED')
            zt1 = zi*fourpi
            t2 = 0.0d0; if (a==b) t2 = 1.0d0
            do iw = 1, wgrid
                zt2 = t2+zt1*sigma(iw)/(w(iw)+eta(iw))
                write(60, '(3G18.10)') t1*w(iw), zt2
            end do
           close(60)
! write the EELS spectra
            write(fname, '("LOSS_", 2I1, ".OUT")') a, b
            write(*, '("  Electron energy loss spectrum written to ", a)') trim(adjustl(fname))
            open(60, file=trim(fname), action='WRITE', form='FORMATTED')
            zt1 = zi*fourpi
            t2 = 0.0d0; if (a==b) t2 = 1.0d0
            do iw = 1, wgrid
                zt2 = t2+zt1*sigma(iw)/(w(iw)+eta(iw))
                write(60, '(2G18.10)') t1*w(iw), -aimag(1.0d0/zt2)
            end do
           close(60)
            if (input%properties%dielmat%tevout) write(*, '("  Output energy is in eV")')
            write(*,*)
        end if
#ifdef MPI        
        call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
    
    end do ! l, optical components


    if (rank==0) then
       
    ! write the optical components of conductivity from Boltzmann     
       write(fname, '("SIGMABOLTZ.OUT")')
       write(*, '("  Optical conductivity components written to ", a)') trim(adjustl(fname))
       open(60, file=trim(fname), action='WRITE', form='FORMATTED')
       !write(*,*) sigmaboltz
       do l = 1, ncomp
          write(60, '(2I1,2G18.10)') epscomp(1,l), epscomp(2,l), sigmaboltz(l)
       end do
      close(60)
    ! write the optical components of conductivity from Kubo
       write(fname, '("SIGMAKUBO.OUT")')
       write(*, '("  Optical conductivity components written to ", a)') trim(adjustl(fname))
       open(60, file=trim(fname), action='WRITE', form='FORMATTED')
       !write(*,*) sigmaboltz
       do l = 1, ncomp
          write(60, '(2I1,2G18.10)') epscomp(1,l), epscomp(2,l), sigmak(l)
       end do
      close(60)
    ! write the optical components of conductivity from Kubo with second derivative
       if (dshift>0) then
          write(fname, '("SIGMAKUBOTWO.OUT")')
          write(*, '("  Optical conductivity components written to ", a)') trim(adjustl(fname))
          open(60, file=trim(fname), action='WRITE', form='FORMATTED')
          !write(*,*) sigmaboltz
          do l = 1, ncomp
             write(60, '(2I1,2G18.10)') epscomp(1,l), epscomp(2,l), sigmaktwo(l)
          end do
         close(60)
       end if
    end if
    
    close(50)
    if (dshift>0) then
       close(80)
       deallocate(pmatd)
    end if
    
    deallocate(pmat)
    deallocate(evalsvt,occsvt)
    deallocate(w,sigma,sigmab,seebeck,thermalcond,zt,eta,td,ndos,echarge)
    !deallocate(nf)
end subroutine
