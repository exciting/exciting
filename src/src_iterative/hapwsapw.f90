! Copyright (C) 2015-2023 exciting team (Berlin and Riga)
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine HapwSapw(n,npw,nwf,system,fftmap,cfir,vir,mir,x,Hx,Sx)
!> The subroutine calculates H|psi> and S|psi> applied to the APW part of the trial wavefunction. 
!>  
      Use physical_constants, Only : alpha
      use modmain
      use modxs
      Use modfvsystem
      Use modmpi, only: terminate_mpi_env, mpiglobal 
      implicit none
      integer, intent(in) :: n              ! leading dimension of |psi>
      integer, intent(in) :: nwf            ! number of trial wavefunctions |psi> 
      integer, intent(in) :: npw            ! size of the APW part
      Type (evsystem) :: system             ! eigensystem datastructure containing either S and H or APW matching coefficients 
      Type (fftmap_type) :: fftmap          ! FFT grid and the mappings
      complex(8), intent(in) :: cfir(*)     ! MT step function on the FFT grid defined in fftmap
      complex(8), intent(in) :: mir(*)      ! Relativistic correction term on the FFT grid defined in fftmap
      complex(8), intent(in) :: vir(*)      ! Potential times the step function on the FFT grid defined in fftmap
      complex(8), intent(in) :: x(n,nwf)    ! |psi>
      complex(8), intent(out) :: Sx(n,nwf)  ! S|psi>
      complex(8), intent(out) :: Hx(n,nwf)  ! H|psi>

! local variables

      complex(8), allocatable :: zax(:,:),zax2(:,:),buf(:,:)
      complex(8), allocatable :: zfft(:)
      integer :: i,ig,ix,LOoffset,if1,if3,ilo,l1,l,m,io1,maxnlo
      integer :: is,ia,ias
      real(8) :: ta,tb,tc,td
      Real(8) :: a2
      Parameter (a2=0.5d0*alpha**2) 

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
      else ! the matrix-free version

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
! -----Applies the MT part of the overlap: S|psi>-----
! Calculates  A|psi>, in other words it expresses the trial wavefunction in terms of radial basis functions.
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

! Applies A^H A|psi>, APW part of the result
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

! The LO part of the result 
            zax2=0d0
            if3=0
            do ilo = 1, nlorb (is)
              l=lorbl (ilo, is)
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

!  -----Applies the MT part of the Hamiltonian: H|psi>-----
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

!----The interstitial part----

! This bit is commented out, because it is a slow implementation with full FFT grids.
! I am reluctant to remove this, because it may be useful for debugging.
!     call timesec(td)
!
!     ! Overlap 
!!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,ig,zfft)
!     allocate(zfft(ngrtot))
!!$OMP DO
!     do i=1,nwf
!
!       zfft=0d0
!       do ig=1,npw
!         zfft(igfft(current_igkig(ig)))=x(ig,i)
!       enddo
!       Call zfftifc (3, ngrid,1, zfft)
!       do ig=1,ngrtot
!         zfft(ig)=zfft(ig)*cfunir(ig)
!       enddo
!       Call zfftifc (3, ngrid,-1, zfft)
!
!       do ig=1,npw
!         Sx(ig,i)=Sx(ig,i)+zfft(igfft(current_igkig(ig)))
!       enddo
!     enddo
!!$OMP END DO
!     deallocate(zfft)
!!$OMP END PARALLEL 
!
!! Kinetic energy
!!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,ix,ig,zfft)
!     allocate(zfft(ngrtot))
!!$OMP DO
!     do i=1,nwf
!       do ix=1,3   
!         zfft=0d0
!         do ig=1,npw
!           zfft(igfft(current_igkig(ig)))=x(ig,i)*current_vgkc(ix,ig)
!         enddo
!
!         Call zfftifc (3, ngrid,1, zfft)
!
!         if (input%groundstate%ValenceRelativity.ne."none") then
!           do ig=1,ngrtot
!             zfft(ig)=zfft(ig)/(1d0-veffir(ig)*a2)*cfunir(ig)
!           enddo
!         else
!           do ig=1,ngrtot
!             zfft(ig)=zfft(ig)*cfunir(ig)
!           enddo
!         endif
!
!         Call zfftifc (3, ngrid,-1, zfft)
!         do ig=1,npw
!           Hx(ig,i)=Hx(ig,i)+0.5d0*zfft(igfft(current_igkig(ig)))*current_vgkc(ix,ig)
!         enddo
!
!       enddo
!     enddo
!!$OMP END DO
!     deallocate(zfft)
!!$OMP END PARALLEL
!
!! potential energy
!!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,ig,zfft)
!     allocate(zfft(ngrtot))
!!$OMP DO
!     do i=1,nwf
!
!       zfft=0d0
!       do ig=1,npw
!         zfft(igfft(current_igkig(ig)))=x(ig,i)
!       enddo
!       Call zfftifc (3, ngrid,1, zfft)
!       do ig=1,ngrtot
!         zfft(ig)=zfft(ig)*veffir(ig)*cfunir(ig)
!       enddo
!       Call zfftifc (3, ngrid,-1, zfft)
!
!       do ig=1,npw
!         Hx(ig,i)=Hx(ig,i)+zfft(igfft(current_igkig(ig)))
!       enddo
!     enddo
!!$OMP END DO
!     deallocate(zfft)
!!$OMP END PARALLEL
!------------- End of the commented out section -------------------

! This is the bit that actually is used.
! FFTs are applied on a smaller grid.

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

! Kinetic energy
!$OMP DO
     do i=1,nwf
       do ix=1,3   
         zfft=0d0
         do ig=1,npw
           zfft(fftmap%igfft(current_igkig(ig)))=x(ig,i)*current_vgkc(ix,ig)
         enddo

         Call zfftifc (3, fftmap%ngrid,1, zfft)

         if (input%groundstate%ValenceRelativity.eq."zora") then
           do ig=1,fftmap%ngrtot
             zfft(ig)=zfft(ig)*mir(ig)
           enddo
         elseif (input%groundstate%ValenceRelativity.eq."none") then
           do ig=1,fftmap%ngrtot
             zfft(ig)=zfft(ig)*cfir(ig)
           enddo
         elseif (input%groundstate%ValenceRelativity.eq."iora*") then
           write(*,*) "Davidson/hapwsapw: matrixless iora* is currently not supported"
           call terminate_mpi_env(mpiglobal) 
         endif

         Call zfftifc (3, fftmap%ngrid,-1, zfft)
         do ig=1,npw
           Hx(ig,i)=Hx(ig,i)+0.5d0*zfft(fftmap%igfft(current_igkig(ig)))*current_vgkc(ix,ig)
         enddo

       enddo
     enddo
!$OMP END DO

! potential energy
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

endif !constructHS


call timesec(tb)

#ifdef TIMINGS
 write(*,*) 'HapwSapw',tb-ta,td-tc
#endif
end subroutine HapwSapw


