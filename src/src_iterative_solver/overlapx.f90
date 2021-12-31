subroutine OverlapX(n,npw,nwf,system,fftmap,cfir,x,Sx)
      use modmain
      use modxs
      Use modfvsystem
      use modmpi
      implicit none
      integer ::n,nwf,npw,nxy
      complex(8) :: x(n,nwf),Sx(n,nwf)
      complex(8), intent(in) :: cfir(*)
      Type (fftmap_type) :: fftmap
      Type (evsystem) :: system

      complex(8), allocatable :: zax(:,:),fftsx(:,:)
      complex(8), allocatable :: zfft(:),buf(:,:),wf(:,:),bufsend(:),bufrecv(:)

      integer :: i,ig,j,k

      integer :: is,ia,ias
      real(8) :: ta,tb,tc,td,rscale,ts1,ts2
      logical, allocatable :: fftpattern(:)
      integer, allocatable :: fftlist(:),streaklo(:),streakhi(:)
      integer :: fftnumber,nstreaks
      integer :: nwf_local,wflo,wfhi,nsteps,blockwf
      complex(8) :: zsum,zsum2

      call timesec(ta)
      Sx=zzero
      if (associated(system%overlap%za)) then
! The overlap matrix is available
        call zgemm('N', &           ! TRANSA = 'C'  op( A ) = A**H.
                   'N', &           ! TRANSB = 'N'  op( B ) = B.
                    npw, &          ! M ... rows of op( A ) = rows of C
                    nwf, &           ! N ... cols of op( B ) = cols of C
                    npw, &          ! K ... cols of op( A ) = rows of op( B )
                    zone, &          ! alpha
                    system%overlap%za, &           ! A
                    system%overlap%rank,&           ! LDA ... leading dimension of A
                    x, &           ! B
                    npw, &          ! LDB ... leading dimension of B
                    zzero, &          ! beta
                    Sx, &  ! C
                    npw &      ! LDC ... leading dimension of C
                    )
      else
! No overlap matrix
        call timesec(tc)

        allocate(zax(mt_hscf%maxaa,nwf))
        allocate(buf(mt_hscf%maxaa,nwf))

! APW-APW
        Do is = 1, nspecies
          Do ia = 1, natoms (is)
            ias = idxas (ia, is)
            zax=zzero
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
            call zgemm('C', &           ! TRANSA = 'C'  op( A ) = A**H.
                       'N', &           ! TRANSB = 'N'  op( B ) = B.
                        npw, &          ! M ... rows of op( A ) = rows of C
                        nwf, &           ! N ... cols of op( B ) = cols of C
                        mt_hscf%maxaa,&  ! K ... cols of op( A ) = rows of op( B )
                        zone, &          ! alpha
                        system%apwi(1,1,ias), &           ! A
                        mt_hscf%maxaa,&           ! LDA ... leading dimension of A
                        zax, &           ! B
                        mt_hscf%maxaa, &          ! LDB ... leading dimension of B
                        zone, &          ! beta
                        Sx(1,1), &  ! C
                        n &      ! LDC ... leading dimension of C
                       )
          enddo
       enddo

       deallocate(zax,buf)
       call timesec(td)
if (.false.) then
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ig,zfft)
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
else
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ig,zfft)
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

endif



     endif

     call timesec(tb)
#ifdef TIMINGS
     write(*,*) 'OverlapX', nwf,tb-ta,td-tc
#endif
end subroutine OverlapX


