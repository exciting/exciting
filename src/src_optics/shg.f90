!BOP
!!ROUTINE: shg
!!INTERFACE:
!
subroutine shg(a,b,c)
!
!!USES:
    use modinput
    use modmain
    use modmpi
    
!!DESCRIPTION:
!   Calculates susceptibility tensor for non-linear optical second-harmonic
!   generation (SHG). The terms (ztm) are numbered according to Eqs. (49)-(51)
!   of the article {\it Physica Scripta} {\bf T109}, 128 (2004). Other good
!   references are {\it Phys. Rev. B} {\bf 48}, 11705 (1993) and
!   {\it Phys. Rev. B} {\bf 53}, 10751 (1996).

!!REVISION HISTORY:
!   Rewrote earlier version, June 2010 (Sharma)
!   Modified, April 2013 (DIN)

!EOP
!BOC
    implicit none

! Optical components
    integer, intent(IN) :: a, b, c

! local variables
    integer :: ik,jk,ist,jst,kst
    integer :: iw,l
    integer :: recl, iostat
    integer :: isym
    real(8) :: sc(3,3), v1(3), v2(3), v3(3)
    integer :: kfirst, klast
    integer :: wgrid
    
! smallest eigenvalue difference allowed in denominator
    real(8)    :: eji,eki,ekj,t1,t2
    complex(8) :: pii(3),dji(3),vji(3),vik(3),vkj(3)
    complex(8) :: eta,ztm(3,3),zt1
    character(256) :: fname
    logical :: exist

! allocatable arrays
    real(8), allocatable :: w(:)
    real(8), allocatable :: evalsvt(:)
    complex(8), allocatable :: pmat(:,:,:)
    complex(8), allocatable :: chiw(:,:)
    complex(8), allocatable :: chi2w(:,:)
    
    if (rank==0) then
      write(*,*)
      write(*,'("Info(shg):")')
    end if
    
! initialise universal variables
    call init0
    call init1

! read Fermi energy from file
    call readfermi

! generate energy grid (starting from zero)
    wgrid = input%properties%shg%wgrid
    allocate(w(wgrid))
    t1 = input%properties%shg%wmax/dble(wgrid)
    do iw = 1, wgrid
      w(iw) = t1*dble(iw-1)
    end do

! calculate the momentum matrix elements
    !if (rank==0) then
    !    write(*,*)
    !    write(*,'("  Calculate the momentum matrix elements")')
    !end if
    !call writepmat

! find the record length for momentum matrix element file
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

! allocate response function arrays
    allocate(chiw(wgrid,3))
    allocate(chi2w(wgrid,2))
    chiw(:,:)=0.d0
    chi2w(:,:)=0.d0

! i divided by the complex relaxation time
    eta = cmplx(0.d0,input%properties%shg%swidth,8)

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

! read momentum matrix elements from file
      allocate(pmat(3,nstsv,nstsv))
      read(50,rec=jk) pmat

! rotate the matrix elements from the reduced to non-reduced k-point
      if (isym > 1) then
        sc(:,:) = symlatc(:,:,lsplsymc(isym))
        do ist = 1, nstsv
          do jst = 1, nstsv
            v1(:) = dble(pmat(:,ist,jst))
            call r3mv(sc,v1,v2)
            v1(:) = aimag(pmat(:,ist,jst))
            call r3mv(sc,v1,v3)
            pmat(:,ist,jst) = cmplx(v2(:),v3(:),8)
          end do
        end do
      end if

! get SV eigenvalues from file
      allocate(evalsvt(nstsv))
      call getevalsv(vkl(:,jk),evalsvt)

! scissor correction applied to the momentum matrix elements
      do ist = 1, nstsv
        if (evalsvt(ist) < efermi) then
          do jst = 1, nstsv
            if (evalsvt(jst) > efermi) then
              eji = evalsvt(jst)-evalsvt(ist)
              t1 = (eji+input%properties%shg%scissor)/eji
              pmat(:,ist,jst) = t1*pmat(:,ist,jst)
            end if
          end do
        end if
      end do
      zt1 = zi*(occmax/omega)/dble(nkptnr)

! loop over valence states
      do ist = 1, nstsv
      
        if (evalsvt(ist) < efermi) then
          pii(:) = pmat(:,ist,ist)

! loop over conduction states
          do jst = 1, nstsv
            
            if (evalsvt(jst) > efermi) then
              eji = evalsvt(jst)-evalsvt(ist)+input%properties%shg%scissor
              vji(:) = pmat(:,jst,ist)
              dji(:) = pmat(:,jst,jst)-pii(:)

! loop over intermediate states
              do kst = 1, nstsv

                if ((kst.ne.ist).and.(kst.ne.jst)) then
                  eki = evalsvt(kst)-evalsvt(ist)
                  ekj = evalsvt(kst)-evalsvt(jst)

                  if (evalsvt(kst) > efermi) then
                    eki = eki+input%properties%shg%scissor
                  else
                    ekj = ekj-input%properties%shg%scissor
                  end if
                  vkj(:) = pmat(:,kst,jst)
                  vik(:) = pmat(:,ist,kst)
! interband terms
                  t1 = -eji*eki*(-ekj)*(eki+ekj)
                  if (abs(t1).gt.input%properties%shg%etol) then
                    t1 = 1.d0/t1
                  else
                    t1 = t1/input%properties%shg%etol**2
                  end if
                  ztm(1,1) = zt1*conjg(vji(a))*(conjg(vkj(b))*conjg(vik(c))+ &
                  &          conjg(vik(b))*conjg(vkj(c)))*t1
                  t1 = eji*(-eki)*ekj*(-eki-eji)
                  if (abs(t1).gt.input%properties%shg%etol) then
                    t1 = 1.d0/t1
                  else
                    t1 = t1/input%properties%shg%etol**2
                  end if
                  ztm(1,2) = 0.5d0*zt1*vkj(c)*(vji(a)*vik(b)+vik(a)*vji(b))*t1
                  t1 = eji*(-eki)*ekj*(ekj-eji)
                  if (abs(t1).gt.input%properties%shg%etol) then
                    t1 = 1.d0/t1
                  else
                    t1 = t1/input%properties%shg%etol**2
                  end if
                  ztm(1,3) = 0.5d0*zt1*vik(b)*(vkj(c)*vji(a)+vji(c)*vkj(a))*t1
! intraband terms
                  t1 = (-eki)*ekj*eji**3
                  if (abs(t1).gt.input%properties%shg%etol) then
                    t1 = 1.d0/t1
                  else
                    t1 = t1/input%properties%shg%etol**2
                  end if
                  ztm(2,1) = 0.5d0*zt1*(eki*vik(b)*(vkj(c)*vji(a)+vji(c)*vkj(a))+ &
                  &          ekj*vkj(c)*(vji(a)*vik(b)+vik(a)*vji(b)))*t1
                  t1 = (eki*(-ekj)*eji**3)/(ekj+eki)
                  if (abs(t1).gt.input%properties%shg%etol) then
                    t1 = 1.d0/t1
                  else
                    t1 = t1/input%properties%shg%etol**2
                  end if
                  ztm(2,3) = zt1*conjg(vji(a))*(conjg(vkj(b))*conjg(vik(c))+ &
                  &          conjg(vik(b))*conjg(vkj(c)))*t1
! modulation terms
                  t1 = ekj*(-eki)*eji**3
                  if (abs(t1).gt.input%properties%shg%etol) then
                    t1 = 1.d0/t1
                  else
                    t1 = t1/input%properties%shg%etol**2
                  end if
                  ztm(3,1) = -0.5d0*zt1*eki*vkj(a)*(vji(b)*vik(c)+vik(b)*vji(c))*t1
                  t1 = ekj*(-eki)*eji**3
                  if (abs(t1).gt.input%properties%shg%etol) then
                    t1 = 1.d0/t1
                  else
                    t1 = t1/input%properties%shg%etol**2
                  end if
                  ztm(3,2) = 0.5d0*zt1*ekj*vik(a)*(vkj(b)*vji(c)+vji(b)*vkj(c))*t1

                  do iw = 1, wgrid
! 2w interband
                    chi2w(iw,1) = chi2w(iw,1)+ztm(1,1)/(eji-2.d0*w(iw)+eta)
! 2w intraband
                    chi2w(iw,2) = chi2w(iw,2)+ztm(2,3)/(eji-2.d0*w(iw)+eta)
! w interband
                    chiw(iw,1) = chiw(iw,1)-(ztm(1,2)-ztm(1,3))/(eji-w(iw)+eta)
! w intraband
                    chiw(iw,2) = chiw(iw,2)+ztm(2,1)/(eji-w(iw)+eta)
! w modulation
                    chiw(iw,3) = chiw(iw,3)+0.5d0*(ztm(3,1)-ztm(3,2))/(eji-w(iw)+eta)
                  end do

                end if
! end loop over kst
              end do
              ztm(2,2) = 4.d0*zt1*conjg(vji(a))*(dji(b)*vji(c)+vji(b)*dji(c))/eji**4
              ztm(3,3) = 0.5d0*zt1*vji(a)*(vji(b)*dji(c)+dji(b)*vji(c))/eji**4

              do iw = 1, wgrid
! 2w intraband
                chi2w(iw,2)=chi2w(iw,2)+ztm(2,2)/(eji-2.d0*w(iw)+eta)
! w modulation
                chiw(iw,3)=chiw(iw,3)+0.5d0*ztm(3,3)/(eji-w(iw)+eta)
              end do

! end loop over jst
            end if
          end do
! end loop over ist
        end if
      end do
      
      deallocate(pmat)
      deallocate(evalsvt)
      
! end loop over k-points
    end do

    close(50)

#ifdef MPI
    call MPI_ALLREDUCE(MPI_IN_PLACE, chiw, 3*wgrid, &
    &                  MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, chi2w, 2*wgrid, &  
    &                  MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

    call barrier
    if (rank==0) then
! output energy units
        t1 = 1.0d0
        if (input%properties%shg%tevout) t1 = h2ev    
! write to files
        write(fname,'("CHI_INTER2w_",3I1,".OUT")') a,b,c
        open(51,file=trim(fname),action='WRITE',form='FORMATTED')
        write(fname,'("CHI_INTRA2w_",3I1,".OUT")') a,b,c
        open(52,file=trim(fname),action='WRITE',form='FORMATTED')
        write(fname,'("CHI_INTERw_",3I1,".OUT")') a,b,c
        open(53,file=trim(fname),action='WRITE',form='FORMATTED')
        write(fname,'("CHI_INTRAw_",3I1,".OUT")') a,b,c
        open(54,file=trim(fname),action='WRITE',form='FORMATTED')
        write(fname,'("CHI_MOD_",3I1,".OUT")') a,b,c
        open(55,file=trim(fname),action='WRITE',form='FORMATTED')
        write(fname,'("CHI_",3I1,".OUT")') a,b,c
        open(56,file=trim(fname),action='WRITE',form='FORMATTED')
        do iw = 1, wgrid
            write(51,'(3G18.10)') t1*w(iw), chi2w(iw,1)
            write(52,'(3G18.10)') t1*w(iw), chi2w(iw,2)
            write(53,'(3G18.10)') t1*w(iw), chiw(iw,1)
            write(54,'(3G18.10)') t1*w(iw), chiw(iw,2)
            write(55,'(3G18.10)') t1*w(iw), chiw(iw,3)
            zt1 = chi2w(iw,1)+chi2w(iw,2)+chiw(iw,1)+chiw(iw,2)+chiw(iw,3)
            write(56,'(4G18.10)') t1*w(iw), zt1, abs(zt1) 
        end do
        close(51); close(52); close(53); close(54); close(55); close(56)
        write(*,*)
        write(*,'("  Susceptibility (complex+module) tensor written to CHI_abc.OUT")')
        write(*,'("  Interband contributions written to CHI_INTERx_abc.OUT")')
        write(*,'("  Intraband contributions written to CHI_INTRAx_abc.OUT")')
        write(*,'("  Modulation contributions written to CHI_MOD_abc.OUT")')
        write(*,'("  for components")')
        write(*,'("  a = ",I1,", b = ",I1,", c = ",I1)') a, b, c
        if (input%properties%shg%tevout) write(*,'("  Output energy is in eV")')
        write(*,*)
    end if

    deallocate(w,chiw,chi2w)

    return
end subroutine
!EOC

