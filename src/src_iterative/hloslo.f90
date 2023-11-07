! Copyright (C) 2015-2023 exciting team (Riga and Berlin)
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
subroutine HloSlo(n,npw,nwf,system,x,Hx,Sx)
!> The subroutine calculates H|psi> and S|psi> applied to the LO part of the trial wavefunction. 
!>  
      Use physical_constants, Only : alpha
      use modmain
      use modxs
      Use modfvsystem
      use modmpi
      implicit none
      integer ::n,nwf,npw
      complex(8) :: x(n,nwf),Sx(n,nwf),Hx(n,nwf)
      Type (evsystem) :: system

      complex(8), allocatable :: zax(:,:),zax2(:,:),zalo(:,:),salo(:,:),zlox(:,:)

      integer :: i,ig,ix,LOoffset,if1,if3,ilo,l1,l,m,io1,ilo2,maxnlo,offset

      integer :: is,ia,ias,gksize,napw
      real(8) :: ta,tb,summa
      Real(8) :: a2
      Parameter (a2=0.5d0*alpha**2) 

      call timesec(ta)

      if (associated(system%hamilton%za)) then
        if (n-npw.ne.0) then
              call zgemm('N', &           ! TRANSA = 'C'  op( A ) = A**H.
                         'N', &           ! TRANSB = 'N'  op( B ) = B.
                          n, &          ! M ... rows of op( A ) = rows of C
                          nwf, &           ! N ... cols of op( B ) = cols of C
                          n-npw, &          ! K ... cols of op( A ) = rows of op( B )
                          zone, &          ! alpha
                          system%hamilton%za(1,npw+1), &           ! A
                          n,&           ! LDA ... leading dimension of A
                          x(npw+1,1), &           ! B
                          n, &          ! LDB ... leading dimension of B
                          zzero, &          ! beta
                          Hx, &  ! C
                          n &      ! LDC ... leading dimension of C
                          )
              call zgemm('N', &           ! TRANSA = 'C'  op( A ) = A**H.
                         'N', &           ! TRANSB = 'N'  op( B ) = B.
                          n, &          ! M ... rows of op( A ) = rows of C
                          nwf, &           ! N ... cols of op( B ) = cols of C
                          n-npw, &      ! K ... cols of op( A ) = rows of op( B )
                          zone, &          ! alpha
                          system%overlap%za(1,npw+1), &           ! A
                          n,&           ! LDA ... leading dimension of A
                          x(npw+1,1), &           ! B
                          n, &          ! LDB ... leading dimension of B
                          zzero, &          ! beta
                          Sx, &  ! C
                          n &      ! LDC ... leading dimension of C
                          )
        endif

      else

        Sx=zzero
        Hx=zzero
        allocate(zax2(mt_hscf%maxaa,nwf))
        maxnlo=mt_hscf%maxnlo

        call timesec(ta)

        LOoffset=npw
        Do is = 1, nspecies
          if (nlorb(is).ne.0) then

!Determine the number of radial functions for given species
            napw=0
            l=input%groundstate%lmaxmat
            do l1=0,l
              napw=napw+(2*l1+1)*apword(l1,is)
            enddo     

            Do ia = 1, natoms (is)
              ias = idxas (ia, is)
              zax2=0d0
              if3=0
              do ilo = 1, nlorb (is)
                l=lorbl (ilo, is)
! APW-LO
                if1=0
                do l1=0,l-1
                  if1=if1+apword(l1,is)*(2*l1+1)
                enddo
                do m=-l,l
                  do io1=1,apword(l,is)
                    if1=if1+1
                    zax2(if1,:)=zax2(if1,:)+oalo(io1,ilo,ias)*x(LOoffset+if3+m+l+1,:)
                  enddo
                enddo
! LO-LO
                if1=0
                Do ilo2 = 1, nlorb (is)
                  If (lorbl(ilo2, is) .Eq. l) Then
                    do m=-l,l
                      Sx(LOoffset+if3+m+l+1,:) =Sx(LOoffset+if3+m+l+1,:)+ololo(ilo, ilo2, ias)*x(LOoffset+if1+m+l+1,:)
                    enddo
                  endif
                  if1=if1+2*lorbl(ilo2, is)+1
                enddo

                if3=if3+2*l+1
              enddo

! more of APW-LO
              call zgemm('C', &           ! TRANSA = 'C'  op( A ) = A**H.
                         'N', &           ! TRANSB = 'N'  op( B ) = B.
                         npw, &          ! M ... rows of op( A ) = rows of C
                         nwf, &           ! N ... cols of op( B ) = cols of C
                         mt_hscf%maxaa,&  ! K ... cols of op( A ) = rows of op( B )
                         zone, &          ! alpha
                         system%apwi(1,1,ias), &           ! A
                         mt_hscf%maxaa,&           ! LDA ... leading dimension of A
                         zax2, &           ! B
                         mt_hscf%maxaa, &          ! LDB ... leading dimension of B
                         zone, &          ! beta
                         Sx(1,1), &  ! C
                         n &      ! LDC ... leading dimension of C
                        )

        
! Now the Hamiltonian
! APW-LO
              zax2=0d0
              call zgemm('C', &           ! TRANSA = 'N'  op( A ) = A.
                         'N', &           ! TRANSB = 'N'  op( B ) = B.
                         napw, &          ! M ... rows of op( A ) = rows of C
                         nwf, &           ! N ... cols of op( B ) = cols of C
                         mt_hscf%losize(is), &          ! K ... cols of op( A ) = rows of op( B )
                         zone, &          ! alpha
                         mt_hscf%main%loa(1,1,ias), &        ! A
                         mt_hscf%maxnlo,&           ! LDA ... leading dimension of A
                         x(LOoffset+1,1), &           ! B
                         n, &          ! LDB ... leading dimension of B
                         zzero, &          ! beta
                         zax2, &  ! C
                         mt_hscf%maxaa &      ! LDC ... leading dimension of C
                        )
              call zgemm('C', &           ! TRANSA = 'C'  op( A ) = A**H.
                         'N', &           ! TRANSB = 'N'  op( B ) = B.
                         npw, &          ! M ... rows of op( A ) = rows of C
                         nwf, &           ! N ... cols of op( B ) = cols of C
                         napw, &         ! K ... cols of op( A ) = rows of op( B )
                         zone, &          ! alpha
                         system%apwi(1,1,ias), &           ! A
                         mt_hscf%maxaa,&           ! LDA ... leading dimension of A
                         zax2, &           ! B
                         mt_hscf%maxaa, &          ! LDB ... leading dimension of B
                         zone, &          ! beta
                         Hx(1,1), &  ! C
                         n &      ! LDC ... leading dimension of C
                        )
! LO-LO
              call zgemm('N', &           ! TRANSA = 'N'  op( A ) = A.
                         'N', &           ! TRANSB = 'N'  op( B ) = B.
                         mt_hscf%losize(is), &          ! M ... rows of op( A ) = rows of C
                         nwf, &           ! N ... cols of op( B ) = cols of C
                         mt_hscf%losize(is), &          ! K ... cols of op( A ) = rows of op( B )
                         zone, &          ! alpha
                         mt_hscf%main%lolo(1,1,ias), &        ! A
                         maxnlo, &           ! LDA ... leading dimension of A
                         x(LOoffset+1,1), &           ! B
                         n, &          ! LDB ... leading dimension of B
                         zone, &          ! beta
                         Hx(LOoffset+1,1), &  ! C
                         n &      ! LDC ... leading dimension of C
                        )

              LOoffset=LOoffset+if3
            enddo
          endif
       enddo

       deallocate(zax2)

     endif


     call timesec(tb)
#ifdef TIMINGS
     if (rank.eq.0) write(*,*) 'HloSlo',tb-ta
#endif
end subroutine HloSlo


