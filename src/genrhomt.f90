!
!BOP
! !ROUTINE: genrhomt
! !INTERFACE:
!
!
Subroutine genrhomt (basis1,basis2,densitymt)
! !USES:
      Use modinput
      Use modmain
      Use modmpi
! !INPUT/OUTPUT PARAMETERS:
!   ik     : k-point number (in,integer)
!   evecfv : first-variational eigenvectors (in,complex(nmatmax,nstfv,nspnfv))
!   evecsv : second-variational eigenvectors (in,complex(nstsv,nstsv))
! !DESCRIPTION:
!   Generates the muffin tin part of the charge density from the density matrix.
!
! !REVISION HISTORY:
!   Created May 2014 (Andris)
!EOP
!BOC
      Implicit None
! arguments
      type(apw_lo_basis_type) :: basis1,basis2
      real (8) :: densitymt(lmmaxvr,nrmtmax,natmtot)
! local variables
      Integer :: nsd, ispn, jspn, is, ia, ias, ist
      Integer :: ir, irc, itp, igk, ifg, i, j, n
      Real (8) :: t1, t2, t3, t4
      Real (8) :: ts0, ts1
      Complex (8) zt1, zt2, zt3
      real(8) :: fr(nrmtmax)
      integer :: l3,lm3,m3,io1,io2,if1,if3,l1,m1,lm1,l2,m2,lm2,if3old,if1old,maxi
      integer :: ilo1,ilo2,maxnlo,maxaa,wfsize
      integer, pointer :: losize(:)
      complex(8), allocatable :: buf(:,:,:)
      real(8), allocatable :: factors(:),rho(:,:,:),rholoc(:,:),factorsnew(:,:),frnew(:,:)

      complex(8), allocatable :: rm(:),rmang(:)
      real(8) :: alpha,a
      parameter (alpha=1d0 / 137.03599911d0)

!write(*,*) sum(mt_dm%main%aa)

      a=0.5d0*alpha**2

      maxnlo=mt_dm%maxnlo
      losize=>mt_dm%losize
      maxaa=mt_dm%maxaa
      wfsize=maxaa+maxnlo

#ifdef MPI
!write(*,*) 'all reduce'
!      allocate(buf(wfsize,wfsize,natmtot))
!      buf(:,:,:)= dcmplx(0d0,0d0)
!      call MPI_ALLREDUCE(mt_dm%main%ff, buf, 2*wfsize*wfsize*natmtot,MPI_DOUBLE_PRECISION,MPI_SUM, MPI_COMM_WORLD,ierr)
!      mt_dm%main%ff=buf
!      deallocate(buf)


!write(*,*) 'all reduce done'
#endif

      allocate(factors(lmmaxvr))


      Call timesec (ts0)
call timesec(t1)


      densitymt=0d0

      Do is = 1, nspecies
        n = lmmaxvr * nrcmt (is)
        allocate(rho(nrmt(is),lmmaxvr,1))

        Do ia = 1, natoms (is)
          ias = idxas (ia, is)
          rho=0d0

! APW-APW part
!if (.true.) then
#ifdef xUSEOMP
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(ilo1,ilo2,if1old,l1,l3,io2,frnew,if1,if3old,factorsnew,m1,lm1,if3,lm3,m3,i,lm2) SHARED(nrmtmax,basis1,basis2,input,is,apword,indgnt,idxlm,gnt,mt_dm,ia,lmmaxvr,nrmt,ias,nlorb,lorbl,losize,maxaa,maxnlo) REDUCTION(+:rho)
#endif
          allocate(frnew(nrmtmax,1))
          allocate(factorsnew(lmmaxvr,2))
          rho=0d0
          if1=0
          Do l1 = 0, input%groundstate%lmaxapw
            Do io1 = 1, apword (l1, is)
              if3=0
              if1old=if1

             Do l3 = 0, input%groundstate%lmaxapw
                Do io2 = 1, apword (l3, is)
                  frnew(:,1)=basis1%apwfr (:, 1, io1, l1, ias) * basis1%apwfr (:, 1, io2, l3, ias)
                  if1=if1old
                  if3old=if3
                  factorsnew=0d0
                  Do m1 = - l1, l1
                    lm1 = idxlm (l1, m1)
                    if1=if1+1
                    if3=if3old
                    Do m3 = - l3, l3
                      lm3 = idxlm (l3, m3)
                      if3=if3+1
                      i=1
                      do while(indgnt(i,lm3,lm1).ne.0) 
                        lm2=indgnt(i,lm3,lm1)
                        factorsnew(lm2,1)=factorsnew(lm2,1)+dble(mt_dm%main%ff(if1,if3,ias)*conjg(listgnt(i,lm3,lm1)))
!                        write(*,*) "indgnt",lm2,conjg(listgnt(i,lm3,lm1))
                        i=i+1
                      enddo                      
                    enddo
                  enddo
#ifdef xUSEOMP
!$OMP DO
#endif
                    do lm2=1,lmmaxvr
                      if (factorsnew(lm2,1).ne.0d0) then
                        rho(1:nrmt(is),lm2,1)=rho(1:nrmt(is),lm2,1)+frnew(1:nrmt(is),1)*factorsnew(lm2,1)
                      endif
                      
                    enddo
#ifdef xUSEOMP
!$OMP END DO NOWAIT
#endif

                End Do
              End Do


            End Do
          End Do

if (losize(is).gt.0) then
!write(*,*) 'howdy'
!APW-LO part

          if1=0
          Do l1 = 0, input%groundstate%lmaxapw
            Do io1 = 1, apword (l1, is)
              if3=0
              if1old=if1
              Do ilo2 = 1, nlorb (is)
                l3 = lorbl (ilo2, is)

                
                frnew(:,1)=basis1%apwfr (:, 1, io1, l1, ias) * basis1%lofr (:, 1, ilo2, ias)
                  if1=if1old
                  if3old=if3
                  factorsnew=0d0
                  Do m1 = - l1, l1
                    lm1 = idxlm (l1, m1)
                    if1=if1+1
                    if3=if3old
                    Do m3 = - l3, l3
                      lm3 = idxlm (l3, m3)
                      if3=if3+1

                      i=1
!                      do lm2=1,lmmaxvr
!                        factorsnew(lm2,1)=factorsnew(lm2,1)+2d0*conjg(gntryy(lm2,lm1,lm3))*dble(mt_dm%main%ff(if1,maxaa+if3,ias))   
!                      enddo                      
                      do while(indgnt(i,lm3,lm1).ne.0)
                        lm2=indgnt(i,lm3,lm1)
                        factorsnew(lm2,1)=factorsnew(lm2,1)+2d0*dble(mt_dm%main%ff(if1,maxaa+if3,ias)*conjg(listgnt(i,lm3,lm1)))
                        i=i+1
                      enddo


                    enddo
                  enddo

#ifdef xUSEOMP
!$OMP DO
#endif
                  do lm2=1,lmmaxvr
                    if (factorsnew(lm2,1).ne.0d0) then
                      rho(1:nrmt(is),lm2,1)=rho(1:nrmt(is),lm2,1)+frnew(1:nrmt(is),1)*factorsnew(lm2,1)
                    endif
                  enddo
#ifdef xUSEOMP
!$OMP END DO NOWAIT
#endif
 
              End Do
            End Do
          End Do


!LO-LO part
            if1=0
            Do ilo1 = 1, nlorb (is)
              l1 = lorbl (ilo1, is)
              if3=0 
              if1old=if1
              Do ilo2 = 1, nlorb (is)
                l3 = lorbl (ilo2, is)
                frnew(:,1)=basis1%lofr (:, 1, ilo1, ias) * basis1%lofr (:, 1, ilo2, ias)
                if1=if1old
                if3old=if3
                factorsnew=0d0
                Do m1 = - l1, l1
                  lm1 = idxlm (l1, m1)
                  if1=if1+1
                  if3=if3old
                  Do m3 = - l3, l3
                    lm3 = idxlm (l3, m3)
                    if3=if3+1

                      i=1
                      do while(indgnt(i,lm3,lm1).ne.0)
                        lm2=indgnt(i,lm3,lm1)
                        factorsnew(lm2,1)=factorsnew(lm2,1)+dble(mt_dm%main%ff(maxaa+if1,maxaa+if3,ias)*conjg(listgnt(i,lm3,lm1)))
                        i=i+1
                      enddo
                      
!                       do lm2=1,lmmaxvr
!                        factorsnew(lm2,1)=factorsnew(lm2,1)+conjg(gntryy(lm2,lm1,lm3))*dble(mt_dm%main%ff(maxaa+if1,maxaa+if3,ias))   
!                      enddo

                  enddo
                enddo

#ifdef xUSEOMP
!$OMP DO
#endif

                do lm2=1,lmmaxvr
                  if (factorsnew(lm2,1).ne.0d0) then
                    rho(1:nrmt(is),lm2,1)=rho(1:nrmt(is),lm2,1)+frnew(1:nrmt(is),1)*factorsnew(lm2,1)
                  endif
                enddo
#ifdef xUSEOMP
!$OMP END DO NOWAIT
#endif

              End Do
            End Do



!          read(*,*) 
endif
          deallocate(frnew,factorsnew)
#ifdef xUSEOMP
!$OMP END PARALLEL
#endif
 



          do lm2=1,lmmaxvr
            densitymt(lm2,1:nrmt(is),ias)=rho(1:nrmt(is),lm2,1)
          enddo

        End Do
        deallocate(rho)
      End Do

call timesec(t2)



!      deallocate(mt_dm%main%ff)
      
      Call timesec (ts1)
!      write(*,*) 'genrhomt',ts1-ts0
!      read(*,*)
      timerho = timerho + ts1 - ts0
      Return
End Subroutine
!EOC

