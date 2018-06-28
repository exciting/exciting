
subroutine shg_new(a,b,c)

    use modinput
    use modmain
    use modmpi
    implicit none
    
    integer, intent(IN) :: a, b, c

    integer :: ik, jk
    integer :: m, n, l
    integer :: iw
    
    integer :: recl, iostat
    integer :: isym
    real(8) :: sc(3,3), v1(3), v2(3), v3(3)
    integer :: kfirst, klast
    integer :: wgrid
    
    character(256) :: fname
    logical :: exist

    real(8), allocatable :: w(:)
    real(8), allocatable :: evalsvt(:), occsvt(:), evalsv0(:)
    complex(8), allocatable :: pmat(:,:,:)
    complex(8), allocatable :: chi2b(:,:)
    complex(8), allocatable :: chi3b(:,:)
    
    real(8) :: t1, t2
    real(8) :: wmn, wml, wln, wdiff
    real(8) :: fnm, fln, fml
    real(8) :: const1, etol
    complex(8) :: pnm(3), dmn(3), pml(3), pln(3)
    complex(8) :: const, prefac, eta, zt1, zt2
    complex(8), allocatable :: eps(:,:)

   
    if (rank==0) then
      write(*,*)
      write(*,'("Info(shg_new):")')
    end if
    
    ! initialization
    call init0
    call init1
    call readfermi

    !----------------------
    ! generate energy grid
    !----------------------
    wgrid = input%properties%shg%wgrid
    allocate(w(wgrid))
    t1 = input%properties%shg%wmax/dble(wgrid)
    do iw = 1, wgrid
      w(iw) = t1*dble(iw-1)
    end do

    !--------------------------
    ! Momentum matrix elements
    !--------------------------
    allocate(pmat(3,nstsv,nstsv))
    inquire(iolength=recl) pmat
    deallocate(pmat)
    open(50,File='PMAT.OUT',Action='READ',Form='UNFORMATTED',Access='DIRECT', &
    &    Recl=recl,IOstat=iostat)
    if (iostat.ne.0) then
      write(*,*)
      write(*,'("Error(shg): error opening PMAT.OUT")')
      write(*,*)
      stop
    end if

    ! complex smearing factor
    eta  = cmplx(0.d0,input%properties%shg%swidth,8)
    etol = input%properties%shg%etol
    
    const = -zi/omega/dble(nkptnr)/occmax
    const1 = 4.0*pi*occmax/omega/dble(nkptnr)
    
    ! allocate response function arrays
    allocate(chi2b(wgrid,2))
    allocate(chi3b(wgrid,2))
    chi2b(:,:) = 0.d0
    chi3b(:,:) = 0.d0
    
    allocate(eps(wgrid,2))
    eps(:,:) = 0.d0

#ifdef MPI
    kfirst = firstofset(rank,nkptnr)
    klast = lastofset(rank,nkptnr)
#else
    kfirst = 1
    klast = nkptnr
#endif    
       
    ! parallel loop over non-reduced k-points
    do ik = kfirst, klast
     
      ! find the k-point number
      call findkpt(vklnr(:,ik),isym,jk)

      !---------------------------------------------
      ! read eigenvalues and occupancies from files
      !---------------------------------------------
      allocate(evalsvt(nstsv))
      call getevalsv(vkl(:,jk),evalsvt)
      
      allocate(occsvt(nstsv))
      call getoccsv(vkl(:,jk),occsvt)

      !-----------------------------------------
      ! read momentum matrix elements from file
      !-----------------------------------------
      allocate(pmat(3,nstsv,nstsv))
      read(50,rec=jk) pmat
      
      ! rotate the matrix elements from the reduced to non-reduced k-point
      if (isym > 1) then
        sc(:,:) = symlatc(:,:,lsplsymc(isym))
        do n = 1, nstsv
          do m = 1, nstsv
            v1(:) = dble(pmat(:,n,m))
            call r3mv(sc,v1,v2)
            v1(:) = aimag(pmat(:,n,m))
            call r3mv(sc,v1,v3)
            pmat(:,n,m) = cmplx(v2(:),v3(:),8)
          end do
        end do
      end if
      
      if (input%properties%shg%scissor > 1.d-4) then

        allocate(evalsv0(nstsv))
        evalsv0(:) = evalsvt(:)
        
        ! scissor correction applied to KS eigenvalues
        do m = 1, nstsv
          if (evalsvt(m) > efermi) then
            evalsvt(m) = evalsvt(m)+input%properties%shg%scissor
          end if
        end do
 
        ! renormalize momentum matrix elements
        do n = 1, nstsv
        if (evalsv0(n) <= efermi) then
          do m = 1, nstsv
          if (evalsv0(m) > efermi) then
            wmn = evalsv0(m)-evalsv0(n)
            if (dabs(wmn) > 1.d-8) then
              t1 = (evalsvt(m)-evalsvt(n))/(evalsv0(m)-evalsv0(n))
              pmat(:,m,n) = t1*pmat(:,m,n)
              pmat(:,n,m) = conjg(pmat(:,m,n))
            end if
          end if
          end do
        end if
        end do
        
        deallocate(evalsv0)
        
      end if
      
      !------------------------
      ! Dielectric function
      !------------------------
      do n = 1, nstsv
      if (evalsvt(n) <= efermi) then
        do m = 1, nstsv
        if (evalsvt(m) > efermi) then
          wmn = evalsvt(m)-evalsvt(n)
          t1 = wmn**2
          if (dabs(t1) > 1.d-8) then
            pnm(:) = pmat(:,n,m)
            zt1 = pnm(a)*conjg(pnm(a))
            prefac = const1*zt1/t1
            do iw = 1, wgrid
              ! resonance
              eps(iw,1) = eps(iw,1) + prefac/( wmn-w(iw)+eta)
              ! anti-resonance
              eps(iw,2) = eps(iw,2) - prefac/(-wmn-w(iw)+eta)
            end do
          end if
        end if
        end do ! m
      end if
      end do ! n
      
      !-------------------------------
      ! Two-bands contribution Eq.(31)
      !-------------------------------
      do m = 1, nstsv
      do n = 1, nstsv
        fnm = occsvt(n)-occsvt(m)
        if (dabs(fnm) > 1.d-4) then
          pnm(:) = pmat(:,n,m)
          dmn(:) = pmat(:,m,m)-pmat(:,n,n)
          zt1    = 0.5*pnm(a)*(dmn(b)*conjg(pnm(c))+dmn(c)*conjg(pnm(b)))
          wmn    = evalsvt(m)-evalsvt(n)
          prefac = const * fnm * zt1 / wmn**4
          do iw = 1, wgrid
            ! 2w-term
            chi2b(iw,2) = chi2b(iw,2) + 8.0*prefac/(wmn-2.0*w(iw)+eta)
            ! w-term
            chi2b(iw,1) = chi2b(iw,1) - prefac/(wmn-w(iw)+eta)
          end do
        end if
      end do
      end do
      
      !----------------------------------
      ! Three-bands contribution Eq. (32)
      !----------------------------------
      do m = 1, nstsv
      do n = 1, nstsv

        if (n==m) cycle
        
        fnm    = occsvt(n)-occsvt(m)
        wmn    = evalsvt(m)-evalsvt(n)
        pnm(:) = pmat(:,n,m)
          
        do l = 1, nstsv
      
          if ((l==n).or.(l==m)) cycle
         
          fln    = occsvt(l)-occsvt(n)
          wln    = evalsvt(l)-evalsvt(n)
          pln(:) = pmat(:,l,n)
        
          fml    = occsvt(m)-occsvt(l)
          wml    = evalsvt(m)-evalsvt(l)
          pml(:) = pmat(:,m,l)

          wdiff  = wln-wml
          
          zt1 = 0.5*pnm(a)*(pml(b)*pln(c)+pml(c)*pln(b))
          prefac = const * zt1
          
          ! 2w term
          if (abs(fnm) > 1.d-4) then
            t1 = wdiff*wmn**3
            if (abs(t1) > etol) then
              zt1 = prefac * fnm / t1
              do iw = 1, wgrid
                chi3b(iw,2) = chi3b(iw,2) + 8.0*zt1/(wmn-2.0*w(iw)+eta)
              enddo
            endif
          endif
        
          ! ml w-term
          if (abs(fml) > 1.d-4) then
            t1 = wdiff*wml**3
            if (abs(t1) > etol) then
              zt1 = prefac * fml / t1
              do iw = 1, wgrid
                chi3b(iw,1) = chi3b(iw,1) + zt1/(wml-w(iw)+eta)
              enddo
            endif
          endif
            
          ! ln w-term
          if (dabs(fln) > 1.d-4) then
            t1 = wdiff*wln**3
            if (abs(t1) > etol) then
              zt1 = prefac * fln / t1
              do iw = 1, wgrid
                chi3b(iw,1) = chi3b(iw,1) + zt1/(wln-w(iw)+eta)
              enddo
            endif
          endif
      
        end do
        
      end do
      end do
      
      deallocate(pmat)
      deallocate(evalsvt)
      deallocate(occsvt)

    ! end loop over k-points
    end do

#ifdef MPI    
    call MPI_ALLREDUCE(MPI_IN_PLACE, eps, 2*wgrid, &
    &                  MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, chi2b, 2*wgrid, &
    &                  MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, chi3b, 2*wgrid, &  
    &                  MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

    call barrier
    if (rank==0) then
        ! output energy units
        t1 = 1.0d0
        if (input%properties%shg%tevout) t1 = h2ev    
         
        ! Epsilon
        write(fname,'("EPSILON_",3I1,".OUT")') a,b,c
        open(50,file=trim(fname),action='WRITE',form='FORMATTED')
        eps(:,1) = 1.d0+eps(:,1)
        do iw = 1, wgrid
          write(50,'(7G18.10)') t1*w(iw), eps(iw,1), eps(iw,2), eps(iw,1)+eps(iw,2)
        end do
        close(50)
        
        ! chi2
        write(fname,'("CHI_2B_",3I1,".OUT")') a,b,c
        open(51,file=trim(fname),action='WRITE',form='FORMATTED')
        write(fname,'("CHI_3B_",3I1,".OUT")') a,b,c
        open(52,file=trim(fname),action='WRITE',form='FORMATTED')
        write(fname,'("CHI_",3I1,".OUT")') a,b,c
        open(53,file=trim(fname),action='WRITE',form='FORMATTED')
        do iw = 1, wgrid
            write(51,'(3G18.10)') t1*w(iw), dble(chi2b(iw,1)+chi2b(iw,2))
            write(52,'(3G18.10)') t1*w(iw), dble(chi3b(iw,1)+chi3b(iw,2))
            t2 = dble(chi2b(iw,1)+chi2b(iw,2)+chi3b(iw,1)+chi3b(iw,2))
            write(53,'(3G18.10)') t1*w(iw), t2
        end do
        write(51,'("     ")')
        write(52,'("     ")')
        write(53,'("     ")')
        do iw = 1, wgrid
            write(51,'(3G18.10)') t1*w(iw), aimag(chi2b(iw,1)+chi2b(iw,2))
            write(52,'(3G18.10)') t1*w(iw), aimag(chi3b(iw,1)+chi3b(iw,2))
            t2 = aimag(chi2b(iw,1)+chi2b(iw,2)+chi3b(iw,1)+chi3b(iw,2))
            write(53,'(3G18.10)') t1*w(iw), t2
        end do
        close(51); close(52)
        write(53,'("     ")')
        do iw = 1, wgrid
          t2 = abs(chi2b(iw,1)+chi2b(iw,2)+chi3b(iw,1)+chi3b(iw,2))
          write(53,'(3G18.10)') t1*w(iw), t2
        end do
        close(53)
        write(*,*)
        write(*,'("  2-bands contributions written to CHI_2B_abc.OUT")')
        write(*,'("  3-bands contributions written to CHI_3B_abc.OUT")')
        write(*,'("  Susceptibility (real,imaginary,module) tensor written to CHI_abc.OUT")')
        write(*,'("  for components")')
        write(*,'("  a = ",I1,", b = ",I1,", c = ",I1)') a, b, c
        if (input%properties%shg%tevout) write(*,'("  Output energy is in eV")')
        write(*,*)
    end if

    deallocate(w,eps,chi2b,chi3b)

    return
end subroutine
!EOC
