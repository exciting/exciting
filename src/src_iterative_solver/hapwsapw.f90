subroutine HapwSapw(n,npw,nwf,system,fftmap,cfir,vir,x,Hx,Sx)
      use modmain
      use modxs
      Use modfvsystem
      use modmpi
      implicit none
      integer ::n,nwf,npw
      complex(8) :: x(n,nwf),Sx(n,nwf),Hx(n,nwf)
      complex(8), intent(in) :: cfir(*)
      complex(8), intent(in) :: vir(*)
      Type (fftmap_type) :: fftmap
      Type (evsystem) :: system

      complex(8), allocatable :: zax(:,:),zax2(:,:),zalo(:,:),salo(:,:),zlox(:,:),buf(:,:)
      complex(8), allocatable :: zfft(:),sx2(:,:),wf(:,:),bufsend(:),bufrecv(:),fftsx(:,:)

      integer :: i,ig,ix,LOoffset,if1,if3,ilo,l1,l,m,io1,ilo2,maxnlo,nxy,j,nwf_local

      integer :: is,ia,ias
      real(8) :: ta,tb,tc,td,ts1,ts2
      Real(8) :: alpha,a2,rscale
      real(8), allocatable :: vtmp(:) 
      Parameter (alpha=1d0 / 137.03599911d0, a2=0.5d0*alpha**2) 
      logical, allocatable :: fftpattern(:)
      integer, allocatable :: fftlist(:),streaklo(:),streakhi(:)
      integer :: fftnumber,nstreaks
      integer :: offset,gksize
      complex(8) :: zsum,zsum2
!stop
call timesec(ta)
Sx=zzero
Hx=zzero
if (associated(system%hamilton%za)) then
      call zgemm('N', &           ! TRANSA = 'C'  op( A ) = A**H.
                 'N', &           ! TRANSB = 'N'  op( B ) = B.
                  n, &          ! M ... rows of op( A ) = rows of C
                  nwf, &           ! N ... cols of op( B ) = cols of C
                  n, &          ! K ... cols of op( A ) = rows of op( B )
                  zone, &          ! alpha
                  system%hamilton%za, &           ! A
                  n,&           ! LDA ... leading dimension of A
                  x, &           ! B
                  n, &          ! LDB ... leading dimension of B
                  zzero, &          ! beta
                  Hx, &  ! C
                  n &      ! LDC ... leading dimension of C
                  )
      call zgemm('N', &           ! TRANSA = 'C'  op( A ) = A**H.
                 'N', &           ! TRANSB = 'N'  op( B ) = B.
                  n, &          ! M ... rows of op( A ) = rows of C
                  nwf, &           ! N ... cols of op( B ) = cols of C
                  n, &          ! K ... cols of op( A ) = rows of op( B )
                  zone, &          ! alpha
                  system%overlap%za, &           ! A
                  n,&           ! LDA ... leading dimension of A
                  x, &           ! B
                  n, &          ! LDB ... leading dimension of B
                  zzero, &          ! beta
                  Sx, &  ! C
                  n &      ! LDC ... leading dimension of C
                  )
else


    Sx=zzero
    Hx=zzero

    allocate(zax(mt_hscf%maxaa,nwf))
    allocate(zax2(mt_hscf%maxaa,nwf))
    allocate(buf(mt_hscf%maxaa,nwf))

    maxnlo=mt_hscf%maxnlo


    call timesec(tc)

      LOoffset=npw
      Do is = 1, nspecies
        Do ia = 1, natoms (is)
          ias = idxas (ia, is)
        call zgemm('N', &           ! TRANSA = 'C'  op( A ) = A**H.
                   'N', &           ! TRANSB = 'N'  op( B ) = B.
                    mt_hscf%maxaa, &          ! M ... rows of op( A ) = rows of C
                    nwf, &           ! N ... cols of op( B ) = cols of C
                    npw, &          ! K ... cols of op( A ) = rows of op( B )
                    zone, &          ! alpha
                    system%apwi(1,1,ias), &           ! A
                    mt_hscf%maxaa,&           ! LDA ... leading dimension of A
                    x(1,1), &           ! B
                    n, &          ! LDB ... leading dimension of B
                    zzero, &          ! beta
                    zax, &  ! C
                    mt_hscf%maxaa &      ! LDC ... leading dimension of C
                   )
#ifdef MPI
!        buf=zzero
!        call MPI_ALLREDUCE(zax, buf, haaijSize*nwf,MPI_DOUBLE_COMPLEX,MPI_SUM, MPI_COMM_WORLD,ierr)
!        zax=buf
#endif
        call zgemm('C', &           ! TRANSA = 'C'  op( A ) = A**H.
                   'N', &           ! TRANSB = 'N'  op( B ) = B.
                    npw, &          ! M ... rows of op( A ) = rows of C
                    nwf, &           ! N ... cols of op( B ) = cols of C
                    mt_hscf%maxaa, &  ! K ... cols of op( A ) = rows of op( B )
                    zone, &          ! alpha
                    system%apwi(1,1,ias), &           ! A
                    mt_hscf%maxaa,&           ! LDA ... leading dimension of A
                    zax, &           ! B
                    mt_hscf%maxaa, &          ! LDB ... leading dimension of B
                    zone, &          ! beta
                    Sx(1,1), &  ! C
                    n &      ! LDC ... leading dimension of C
                   )


        zax2=0d0
        if3=0
        do ilo = 1, nlorb (is)
          l=lorbl (ilo, is)
! LO-APW
          if1=0
          do l1=0,l-1
            if1=if1+apword(l1,is)*(2*l1+1)
          enddo
          do m=-l,l
            do io1=1,apword(l,is)
              if1=if1+1
              Sx(LOoffset+if3+m+l+1,:) =Sx(LOoffset+if3+m+l+1,:)+oalo(io1,ilo,ias)*zax(if1,:)
            enddo
          enddo
          if3=if3+2*l+1
        enddo

! Now the Hamiltonian
! APW-APW
          call zgemm('N', &           ! TRANSA = 'N'  op( A ) = A.
                     'N', &           ! TRANSB = 'N'  op( B ) = B.
                      mt_hscf%maxaa, &          ! M ... rows of op( A ) = rows of C
                      nwf, &           ! N ... cols of op( B ) = cols of C
                      mt_hscf%maxaa, &          ! K ... cols of op( A ) = rows of op( B )
                      zone, &          ! alpha
                      mt_hscf%main%aa(1,1,ias), &        ! A
                      mt_hscf%maxaa,&           ! LDA ... leading dimension of A
                      zax, &           ! B
                      mt_hscf%maxaa, &          ! LDB ... leading dimension of B
                      zzero, &          ! beta
                      zax2, &  ! C
                      mt_hscf%maxaa &      ! LDC ... leading dimension of C
                      )

          call zgemm('C', &           ! TRANSA = 'C'  op( A ) = A**H.
                     'N', &           ! TRANSB = 'N'  op( B ) = B.
                      npw, &          ! M ... rows of op( A ) = rows of C
                      nwf, &           ! N ... cols of op( B ) = cols of C
                      mt_hscf%maxaa, &          ! K ... cols of op( A ) = rows of op( B )
                      zone, &          ! alpha
                      system%apwi(1,1,ias), &           ! A
                      mt_hscf%maxaa,&           ! LDA ... leading dimension of A
                      zax2, &           ! B
                      mt_hscf%maxaa, &          ! LDB ... leading dimension of B
                      zone, &          ! beta
                      Hx(1,1), &  ! C
                      n &      ! LDC ... leading dimension of C
                     )

          if (nlorb(is).ne.0) then
! LO-APW
          call zgemm('N', &           ! TRANSA = 'N'  op( A ) = A.
                     'N', &           ! TRANSB = 'N'  op( B ) = B.
                      mt_hscf%losize(is), &          ! M ... rows of op( A ) = rows of C
                      nwf, &           ! N ... cols of op( B ) = cols of C
                      mt_hscf%maxaa, &          ! K ... cols of op( A ) = rows of op( B )
                      zone, &          ! alpha
                      mt_hscf%main%loa(1,1,ias), &        ! A
                      maxnlo,&           ! LDA ... leading dimension of A
                      zax, &           ! B
                      mt_hscf%maxaa, &          ! LDB ... leading dimension of B
                      zone, &          ! beta
                      Hx(LOoffset+1,1), &  ! C
                      n &      ! LDC ... leading dimension of C
                      )
           endif

          LOoffset=LOoffset+if3
        enddo
     enddo


     deallocate(zax,zax2,buf)

!write(*,*) Sx(1,1:10)

!stop

if (.false.) then
     call timesec(td)

! Overlap 
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,ig,zfft)
allocate(zfft(ngrtot))
!$OMP DO
do i=1,nwf

    zfft=0d0
    do ig=1,npw
     zfft(igfft(current_igkig(ig)))=x(ig,i)
    enddo
    Call zfftifc (3, ngrid,1, zfft)
    do ig=1,ngrtot
     zfft(ig)=zfft(ig)*cfunir(ig)
    enddo
    Call zfftifc (3, ngrid,-1, zfft)

    do ig=1,npw
     Sx(ig,i)=Sx(ig,i)+zfft(igfft(current_igkig(ig)))
    enddo
enddo
!$OMP END DO
deallocate(zfft)
!$OMP END PARALLEL 

! Kinetic energy
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,ix,ig,zfft)
allocate(zfft(ngrtot))
!$OMP DO
do i=1,nwf
  do ix=1,3   
    zfft=0d0
    do ig=1,npw
      zfft(igfft(current_igkig(ig)))=x(ig,i)*current_vgkc(ix,ig)
    enddo

    Call zfftifc (3, ngrid,1, zfft)

    if (input%groundstate%ValenceRelativity.ne."none") then
      do ig=1,ngrtot
       zfft(ig)=zfft(ig)/(1d0-veffir(ig)*a2)*cfunir(ig)
      enddo
    else
      do ig=1,ngrtot
       zfft(ig)=zfft(ig)*cfunir(ig)
      enddo
    endif

    Call zfftifc (3, ngrid,-1, zfft)
    do ig=1,npw
     Hx(ig,i)=Hx(ig,i)+0.5d0*zfft(igfft(current_igkig(ig)))*current_vgkc(ix,ig)
    enddo

  enddo
enddo
!$OMP END DO
deallocate(zfft)
!$OMP END PARALLEL

! potential energy
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,ig,zfft)
allocate(zfft(ngrtot))
!$OMP DO
do i=1,nwf

    zfft=0d0
    do ig=1,npw
     zfft(igfft(current_igkig(ig)))=x(ig,i)
    enddo
    Call zfftifc (3, ngrid,1, zfft)
    do ig=1,ngrtot
     zfft(ig)=zfft(ig)*veffir(ig)*cfunir(ig)
    enddo
    Call zfftifc (3, ngrid,-1, zfft)

    do ig=1,npw
     Hx(ig,i)=Hx(ig,i)+zfft(igfft(current_igkig(ig)))
    enddo
enddo
!$OMP END DO
deallocate(zfft)
!$OMP END PARALLEL

else

! Overlap 
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,ig,zfft)
allocate(zfft(fftmap%ngrtot))
!$OMP DO
do i=1,nwf

    zfft=0d0
    do ig=1,npw
     zfft(fftmap%igfft(current_igkig(ig)))=x(ig,i)
    enddo
    Call zfftifc (3, fftmap%ngrid,1, zfft)
    do ig=1,fftmap%ngrtot
     zfft(ig)=zfft(ig)*cfir(ig)
    enddo
    Call zfftifc (3, fftmap%ngrid,-1, zfft)

    do ig=1,npw
     Sx(ig,i)=Sx(ig,i)+zfft(fftmap%igfft(current_igkig(ig)))
    enddo
enddo
!$OMP END DO
deallocate(zfft)
!$OMP END PARALLEL 



! Kinetic energy
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,ix,ig,zfft)
allocate(zfft(fftmap%ngrtot))
!$OMP DO
do i=1,nwf
  do ix=1,3   
    zfft=0d0
    do ig=1,npw
      zfft(fftmap%igfft(current_igkig(ig)))=x(ig,i)*current_vgkc(ix,ig)
    enddo

    Call zfftifc (3, fftmap%ngrid,1, zfft)

    if (input%groundstate%ValenceRelativity.ne."none") then
      write(*,*) 'davidson.f90/hapwsapw.f90: the relativistic correction needs to be specified properly in the code'
      stop
      do ig=1,fftmap%ngrtot
! vir includes cfir already, but it should not. Should we use meffir?
       zfft(ig)=zfft(ig)/(1d0-vir(ig)*a2)*cfir(ig)
      enddo
    else
      do ig=1,fftmap%ngrtot
       zfft(ig)=zfft(ig)*cfir(ig)
      enddo
    endif

    Call zfftifc (3, fftmap%ngrid,-1, zfft)
    do ig=1,npw
     Hx(ig,i)=Hx(ig,i)+0.5d0*zfft(fftmap%igfft(current_igkig(ig)))*current_vgkc(ix,ig)
    enddo

  enddo
enddo
!$OMP END DO
deallocate(zfft)
!$OMP END PARALLEL

! potential energy
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,ig,zfft)
allocate(zfft(fftmap%ngrtot))
!$OMP DO
do i=1,nwf

    zfft=0d0
    do ig=1,npw
     zfft(fftmap%igfft(current_igkig(ig)))=x(ig,i)
    enddo
    Call zfftifc (3, fftmap%ngrid,1, zfft)
    do ig=1,fftmap%ngrtot
     zfft(ig)=zfft(ig)*vir(ig)
    enddo
    Call zfftifc (3, fftmap%ngrid,-1, zfft)

    do ig=1,npw
     Hx(ig,i)=Hx(ig,i)+zfft(fftmap%igfft(current_igkig(ig)))
    enddo
enddo
!$OMP END DO
deallocate(zfft)
!$OMP END PARALLEL

endif !false

endif !constructHS


!write(*,*) 'sum(Hx)',sum(Hx)
!write(*,*) 'sum(Sx)',sum(Sx)
call timesec(tb)

#ifdef TIMINGS
write(*,*) 'HapwSapw',tb-ta,td-tc
!stop
#endif
end subroutine HapwSapw

