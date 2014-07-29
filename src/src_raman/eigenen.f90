! Copyright (C) 2014 exciting team 
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
! History:
! created 1994, Claudia Draxl 
! adapted for exciting 2014, Stefan Kontur
!
!
subroutine EIGENen
!
! solve oscillator eigen problem
!
   use modinput
   use raman_trmat
   use raman_inter, only : xa,h
   use raman_ew
   use raman_fij, only : fi
   use mod_misc, only: filext
   implicit none
   integer :: iidim,i2,i3,i4,ifail,ip,iii,mmm,ind,i,ii,iint,iv,jj
   real(8) :: wnorm
   real(8), allocatable :: zz(:,:),z(:,:)
   real(8), allocatable :: alfi(:),alfr(:),beta(:),work(:)
   real(8) :: u(4),v(4)
   character(3) :: ext
!
!
! allocate local arrays
allocate( z(2*input%properties%raman%ninter,2*input%properties%raman%ninter) )
allocate( zz(2*input%properties%raman%ninter,2*input%properties%raman%ninter) )
allocate( alfi(2*input%properties%raman%ninter) )
allocate( alfr(2*input%properties%raman%ninter) )
allocate( beta(2*input%properties%raman%ninter) )
allocate( work(16*input%properties%raman%ninter) )
!
if (input%properties%raman%writefunc) then
   do i = 1, input%properties%raman%nstate
      write(ext,'(i3.3)') i
      open(unit=250+i,file='RAMAN_EIGFUNC_'//trim(ext)//trim(filext),status='unknown',form='formatted')
   enddo
   open(unit=250,file='RAMAN_EIGVAL'//trim(filext),status='unknown',form='formatted')
endif
      iidim = 2*input%properties%raman%ninter
      write(66,'(//,52("*"),"   EIGEN   ",53("*"),// &
     & " Number of meshpoints: ",I5, 8x,"Matrix size:",I5)') input%properties%raman%ninter+1,iidim
!
      call FIJK
      call matrix_tr(iidim)
      i2 = 2*input%properties%raman%ninter
      i3 = 2*input%properties%raman%ninter
      i4 = 2*input%properties%raman%ninter
!
!  call LAPACK routine DGGEV which computes eigenvalues
!  and right eigenvectors for A * v(j) = lambda(j) * B * v(j)
!  where A and B are N-by-N real nonsymmetric matrices
!  eigenvalues lambda are alpha(re,im)/beta (both can be zero)
!  A = T (= S), B = R (= O)
!  wave functions are obtained as psi(xi) = Sum z_i P_i(xi)
!  
      call DGGEV('N','V',iidim,T,i2,R,i3,alfr,alfi,beta,zz,i4,z,i4, &
     &                     work,16*input%properties%raman%ninter,ifail)
      if (ifail .ne. 0) then
         write (66,*)
         write (66, '("Error(Eigenen): Error in DGGEV routine. INFO = ")') ifail
         stop
      else
         ip = 0
         do iii = 1, iidim
            if (beta(iii) .eq. 0.0d0) then
               write (66, '("Error(Eigenen): LAMBDA is infinite.")')
            else
               if (alfi(iii) .eq. 0.0d0) then       ! if j-th eigenvalue is real, then v(j) = VR(:,j)
                  eigen(iii) = alfr(iii)/beta(iii)
!                 write(66,*) alfr(iii),beta(iii),eigen(iii) 
                  do mmm = 1,iidim
                     zz(mmm,iii) = z(mmm,iii)
                  enddo
               else                                 ! if j-th and (j+1)-th eigenvalues form complex conj pair, then v(j) = VR(:,j) + i*VR(:,j+1)
!                                                   ! [ and v(j+1) = VR(:,j) - i*VR(:,j+1) ]
                  eigen(iii) = alfr(iii)/beta(iii)
!                 write(66,*)eigen(iii) 
                  do mmm = 1,iidim
                     zz(mmm,iii) = z(mmm,iii-ip)
                  enddo
                  ip = 1 - ip
               endif
            endif
         enddo
      endif
!
      call sort(iidim)                              !  sorts eigenvalues and functions ascending
      do i = 1, input%properties%raman%nstate                              !  loop over wave functions
         e1(i) = eigen(i)*1000.0d0                    ! eigenvalue in mHa
         e2(i) = eigen(i)*fhaev                       ! in eV
         e3(i) = eigen(i)*fhawn                       ! in cm-1
         if (i .gt. 1) de(i) = e3(i) - e3(i-1)
         z1(1,i) = 0.d0                             !  z1 refers to real part of the eigen function, z2 to its derivative, according to
         z2(1,i) = zz(1,indexi(i))                  !  definition of eigenfunctions Psi(0) = u1, Psi'(0) = u2, Psi(1) = u3, Psi'(1) = u4
         do iv = 1,input%properties%raman%ninter-1
            z1(iv+1,i) = zz(2*iv,indexi(i))         !  indexi stems from subroutine sort
            z2(iv+1,i) = zz(2*iv+1,indexi(i))
         enddo 
         z1(input%properties%raman%ninter+1,i) = 0.d0
         z2(input%properties%raman%ninter+1,i) = zz(2*input%properties%raman%ninter,indexi(i))
!                                                   ! Normierung der Eigenfunktion:
         wnorm = 0.d0                               !  compute integrals    < i | i >
         do iint = 1,input%properties%raman%ninter  !  loop over intervals
            u(1) = z1(iint,i)                       !  note that u3(int) = u1(int+1), and u4(int) = u2(int+1)
            u(2) = z2(iint,i)
            u(3) = z1(iint+1,i)
            u(4) = z2(iint+1,i)
            v(1) = z1(iint,i)
            v(2) = z2(iint,i)
            v(3) = z1(iint+1,i)
            v(4) = z2(iint+1,i)
            do ii = 1,4
               do jj = 1,4
                  wnorm = wnorm + u(ii)*v(jj)*FI(ii,jj,0)
               enddo
            enddo
         enddo
         wnorm = wnorm*h
         wnorm = 1.d0/dsqrt(wnorm)
         if (input%properties%raman%writefunc) then
            write(66,*)
            write(66,*) 'eigen function ',i
            if (input%properties%raman%writefunc) write(250,*) eigen(i)
         endif
         do ind = 1,input%properties%raman%ninter+1
            z1(ind,i) = z1(ind,i)*wnorm
            z2(ind,i) = z2(ind,i)*wnorm
            if (input%properties%raman%writefunc .and. (ind .eq. input%properties%raman%ninter+1)) then
               write(66,43) xa(input%properties%raman%ninter)+h,z1(ind,i),z2(ind,i)
               if (input%properties%raman%writefunc) &
                      &  write(250+i,43) xa(input%properties%raman%ninter)+h,z1(ind,i)*0.001d0,z2(ind,i)*0.001d0
            elseif (input%properties%raman%writefunc) then
               write(66,43) xa(ind),z1(ind,i),z2(ind,i)
               if (input%properties%raman%writefunc) write(250+i,43) xa(ind),z1(ind,i)*0.001d0,z2(ind,i)*0.001d0
            endif
43          format(3(f24.12))
         enddo
      enddo
!40    continue
!
!     OUTPUT for all dimensions
!
      write(66, '(//,23x,"EIGENVALUES",22x,"TRANSITION ENERGY", &
                & //,7x,"[ mHa ]",11x,"[ eV ]",12x,"[ cm-1 ]",12x,"[ cm-1 ]",/)')
      write(66,'(f14.3,5x,f12.5,8x,f12.2,8x)') e1(1),e2(1),e3(1)
      do i = 2, input%properties%raman%nstate
         write(66,'(f14.3,5x,f12.5,8x,2(f12.2,8x))') e1(i),e2(i),e3(i),de(i)
      enddo
      write(66,*)
!     
deallocate( z,zz,alfi,alfr,beta,work )
if (input%properties%raman%writefunc) then
   do i = 1, input%properties%raman%nstate
      close(250+i)
   enddo
   close(250)
endif
!
      return
      end
!
!
