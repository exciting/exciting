
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

Subroutine emat_wannier( ik, ikq, ist1, nst1, ist2, nst2, emat)
      Use modmain
      Use modinput
      Implicit None
! arguments
      Integer, Intent (In) :: ik, ikq, ist1, nst1, ist2, nst2
      Complex (8), Intent (Out) :: emat( nst1, nst2)
! local variables
      Integer :: l1, l2, l3, m1, m2, m3, o1, o2, ilo1, ilo2, idxlo1, idxlo2
      Integer :: lm1, lm2, lm3, ist, jst, lmo, idxlmo1, idxlmo2
      Integer :: is, ia, ias, i1, i2, i3
      Integer :: ngkq, igk, ifg, ir
      Integer :: i 
      integer :: nr
      integer :: idxlostart
      Real (8) :: vecqc(3), qc, tp(2)
      Real (8) :: vkqc(3), x, t1, gnt
      Real (8) :: v1(3), v2(3), v3(3), vecql(3)
      complex(8) :: expqr
! automatic arrays
      Complex (8) ylm( lmmaxvr)
! allocatable arrays
      Integer, Allocatable :: igkqig(:)
      integer, allocatable :: idxgnt(:,:,:)
      integer, allocatable :: lmo2l(:), lmo2m(:), lmo2o(:), lm2l(:)
      Real (8), Allocatable :: jlqr(:,:)
      real(8), allocatable :: fr(:), gf(:), cf(:,:)
      real(8), allocatable :: vgkql(:,:), vgkqc(:,:), gkqc(:), tpgkqc(:,:)
      real(8), allocatable :: iraa(:,:,:,:,:), iral(:,:,:,:), irll(:,:,:)
      complex(8), allocatable :: listgnt(:,:,:)
      complex(8), allocatable :: pref(:)
      complex(8), allocatable :: mataa(:,:), matal(:,:), matla(:,:), matll(:,:)
      complex(8), allocatable :: irgntaa(:,:), irgntal(:,:), irgntla(:,:), irgntll(:,:)
      complex(8), allocatable :: blockmt(:,:)
      complex(8), allocatable :: match_combined1(:,:), match_combined2(:,:)
      complex(8), allocatable :: auxmat(:,:), auxvec(:)
      Complex (8), Allocatable :: apwalm1(:,:,:,:), apwalm2(:,:,:,:)
      Complex (8), Allocatable :: evecfv1(:,:), evecfv2(:,:)
      Complex (8), Allocatable :: evecsv1(:,:), evecsv2(:,:)
      Complex (8), Allocatable :: wfir(:)
      Complex (8), Allocatable :: zfir1(:), zfir2(:)
      Complex (8), Allocatable :: em(:,:)
      Logical :: first_call = .true.
      Save :: first_call
! external functions
      Real (8) :: gaunt
      Complex (8) zdotc
      External gaunt, zdotc
! check if q-vector is zero
      vecql(:) = vkl( :, ikq) - vkl( :, ik)
      t1 = vecql(1)**2 + vecql(2)**2 + vecql(3)**2
      If( t1 .Lt. input%structure%epslat) Then
         emat(:,:) = 0.d0
         Do i = 1, min( nst1, nst2)
            emat( i, i) = 1.d0
         End Do
         Return
      End If
! check q-vector is commensurate with k-point grid
      v1 (:) = dble( input%groundstate%ngridk(:))*vecql(:)
      v2 (:) = abs( v1(:)-Nint(v1(:)))
      If( (v2(1) .Gt. input%structure%epslat) .Or. &
        & (v2(2) .Gt. input%structure%epslat) .Or. &
        & (v2(3) .Gt. input%structure%epslat)) Then
         Write (*,*)
         Write (*, '("Error(emat_wannier): q-vector incommensurate with k-point grid")')
         Write (*, '(" ngridk : ", 3I6)') input%groundstate%ngridk
         Write (*, '(" vecql : ", 3G18.10)') vecql
         Write (*,*)
         Stop
      End If
! allocate local arrays
      Allocate( igkqig( ngkmax))
      Allocate( jlqr( 0:input%groundstate%lmaxvr, nrmtmax))
      Allocate( apwalm1( ngkmax, apwordmax, lmmaxapw, natmtot))
      Allocate( apwalm2( ngkmax, apwordmax, lmmaxapw, natmtot))
      Allocate( evecfv1( nmatmax, nstfv))
      Allocate( evecfv2( nmatmax, nstfv))
      If (input%groundstate%tevecsv) Then
         Allocate( evecsv1( nstsv, nstsv))
         Allocate( evecsv2( nstsv, nstsv))
      End If
      Allocate( wfir( ngkmax))
      Allocate( zfir1( ngrtot), zfir2( ngrtot))
      Allocate( em( nstfv, nstfv))
! q-vector in Cartesian coordinates
      Call r3mv( bvec, vecql, vecqc)
! length and spherical coordinates of q-vector
      Call sphcrd( vecqc, qc, tp)
! generate the conjugate spherical harmonics of the q-vector
      Call genylm( input%groundstate%lmaxvr, tp, ylm)
      ylm(:) = conjg( ylm(:))
! get the eigenvector for k-point k
      Call getevecfv( vkl( :, ik), vgkl( :, :, :, ik), evecfv1)
! get the eigenvector for k-point k+q
      Call getevecfv( vkl( :, ikq), vgkl( :, :, :, ikq), evecfv2)
! find the matching coefficients for k-point k
      Call match( ngk( 1, ik), gkc( :, 1, ik), tpgkc( :, :, 1, ik), sfacgk(:, :, 1, ik), apwalm1)
! find the matching coefficients for k-point k+q
      Call match( ngk( 1, ikq), gkc( :, 1, ikq), tpgkc( :, :, 1, ikq), sfacgk(:, :, 1, ikq), apwalm2)
! generate radial functions
      call genapwfr
      call genlofr
! set the first-variational matrix element array to zero
      em(:,:) = zzero
! generate Gaunt coefficients
      allocate( idxgnt( input%groundstate%lmaxapw+1, lmmaxapw, lmmaxapw))
      allocate( listgnt( input%groundstate%lmaxapw+1, lmmaxapw, lmmaxapw))
      idxgnt(:,:,:) = 0
      listgnt(:,:,:) = 0.d0 
      do l3 = 1, input%groundstate%lmaxapw
        do m3 = -l3, l3
          lm3 = idxlm( l3, m3)
          
          do l1 = 1, input%groundstate%lmaxapw
            do m1 = -l1, l1
              lm1 = idxlm( l1, m1)
              i = 0
            
              do l2 = 1, input%groundstate%lmaxvr
                do m2 = -l2, l2
                  lm2 = idxlm( l2, m2)
                  gnt = gaunt( l1, l2, l3, m1, m2, m3)
                  if( .not.(gnt .eq. 0.d0)) then
                    i = i + 1
                    listgnt( i, lm1, lm3) = gnt
                    idxgnt( i, lm1, lm3) = lm2
                  end if
                end do
              end do

            end do
          end do
      
        end do
      end do
! compute prefactors of Reyleigh expansion
      allocate( pref( (input%groundstate%lmaxvr + 1)**2))
      allocate( lm2l( (input%groundstate%lmaxvr + 1)**2))
      do l3 = 0, input%groundstate%lmaxvr
        do m3 = -l3, l3
          lm3 = idxlm( l3, m3)
          pref( lm3) = fourpi*zil( l3)*ylm( lm3)
          lm2l( lm3) = l3
        end do
      end do

!------------------------------------!
!     muffin-tin matrix elements     !
!------------------------------------!

! allocate matrices
      ! radial integral
      allocate( iraa( 0:input%groundstate%lmaxvr, 0:input%groundstate%lmaxapw, apwordmax, 0:input%groundstate%lmaxapw, apwordmax), &
                iral( 0:input%groundstate%lmaxvr, 0:input%groundstate%lmaxapw, apwordmax, nlomax), &  
                irll( 0:input%groundstate%lmaxvr, nlomax, nlomax))
      ! sum of different contributions over all muffin tins
      allocate( mataa( ngk( 1, ik), ngk( 1, ikq)), &
                matal( ngk( 1, ik), nlotot), &
                matla( nlotot, ngk( 1, ikq)), &
                matll( nlotot, nlotot))
      ! radial integrals times Gaunt times prefactor and sum over
      allocate( irgntaa( lmmaxapw*apwordmax, lmmaxapw*apwordmax), &
                irgntal( lmmaxapw*apwordmax, nlotot), &
                irgntla( nlotot, lmmaxapw*apwordmax), &
                irgntll( nlotot, nlotot))
      ! matching coefficients with combined l, m, apword index
      allocate( match_combined1( lmmaxapw*apwordmax, ngkmax), &
                match_combined2( lmmaxapw*apwordmax, ngkmax))
      ! block matrix for calculation of final total muffin tin contribution
      allocate( blockmt( ngk( 1, ik)+nlotot, ngk( 1, ikq)+nlotot))
      ! auxiliary matrix and vector
      allocate( auxmat( lmmaxapw*apwordmax, ngkmax))
      allocate( auxvec( ngk( 1, ik)+nlotot))
      allocate( fr( nrmtmax), gf( nrmtmax), cf( 3, nrmtmax))
      ! index maps
      allocate( lmo2l( lmmaxapw*apwordmax), &
                lmo2m( lmmaxapw*apwordmax), &
                lmo2o( lmmaxapw*apwordmax))

      mataa(:,:) = zzero
      matal(:,:) = zzero
      matla(:,:) = zzero
      matll(:,:) = zzero
      idxlostart = 0
! start loops over atoms and species
      do is = 1, nspecies
        nr = nrmt( is)
! compute the spherical Bessel functions
        do ir = 1, nr
          x = qc*spr( ir, is)
          call sbessel( input%groundstate%lmaxvr, x, jlqr( :, ir))
        end do
        do ia = 1, natoms( is)
          expqr = exp( zi*dot_product( vecqc, atposc( :, ia, is))) 
          ias = idxas( ia, is)
! calculate radial integrals
! APW-APW
          do l2 = 0, input%groundstate%lmaxapw
            do o2 = 1, apword( l2, is)
              do l1 = 0, input%groundstate%lmaxapw
                do o1 = 1, apword( l1, is)
                  do l3 = 0, input%groundstate%lmaxvr
                    do ir = 1, nr
                      fr( ir) = apwfr( ir, 1, o1, l1, ias)*jlqr( l3, ir)*apwfr( ir, 1, o2, l2, ias)*spr( ir, is)**2
                    end do
                    call fderiv( -1, nr, spr( :, is), fr, gf, cf)
                    iraa( l3, l1, o1, l2, o2) = gf( nr)
                  end do
                end do
              end do
            end do
          end do
! APW-LO
          do ilo1 = 1, nlorb( is) 
            do l1 = 0, input%groundstate%lmaxapw
              do o1 = 1, apword( l1, is)
                do l3 = 0, input%groundstate%lmaxvr
                  do ir = 1, nr
                    fr( ir) = apwfr( ir, 1, o1, l1, ias)*jlqr( l3, ir)*lofr( ir, 1, ilo1, ias)*spr( ir, is)**2
                  end do
                  call fderiv( -1, nr, spr( :, is), fr, gf, cf)
                  iral( l3, l1, o1, ilo1) = gf( nr)
                end do
              end do
            end do
          end do
! LO-LO
          do ilo2 = 1, nlorb( is) 
            do ilo1 = 1, nlorb( is)
              do l3 = 0, input%groundstate%lmaxvr
                do ir = 1, nr
                  fr( ir) = lofr( ir, 1, ilo1, ias)*jlqr( l3, ir)*lofr( ir, 1, ilo2, ias)*spr( ir, is)**2
                end do
                call fderiv( -1, nr, spr( :, is), fr, gf, cf)
                irll( l3, ilo1, ilo2) = gf( nr)
              end do
            end do
          end do
! write combined matching coefficients
          match_combined1(:,:) = zzero
          match_combined2(:,:) = zzero
          lmo = 0
          do l1 = 0, input%groundstate%lmaxapw
            do o1 = 1, apword( l1, is)
              do m1 = -l1, l1
                lmo = lmo + 1
                lmo2l( lmo) = l1
                lmo2m( lmo) = m1
                lmo2o( lmo) = o1
                lm1 = idxlm( l1, m1)
                match_combined1( lmo, :) = apwalm1( :, o1, lm1, ias)
                match_combined2( lmo, :) = apwalm2( :, o1, lm1, ias)
              end do
            end do
          end do
! multiply radial integrals with prefactor and Gaunt and sum over
          irgntaa(:,:) = zzero
          irgntal(:,:) = zzero
          irgntla(:,:) = zzero
          irgntll(:,:) = zzero
! APW-APW
          do idxlmo1 = 1, lmo
            l1 = lmo2l( idxlmo1)
            m1 = lmo2m( idxlmo1)
            o1 = lmo2o( idxlmo1)
            lm1 = idxlm( l1, m1)
            do idxlmo2 = 1, lmo
              l2 = lmo2l( idxlmo2)
              m2 = lmo2m( idxlmo2)
              o2 = lmo2o( idxlmo2)
              lm2 = idxlm( l2, m2)
              i = 1
              do while( idxgnt( i, lm2, lm1) .ne. 0)
                lm3 = idxgnt( i, lm2, lm1)
                irgntaa( idxlmo1, idxlmo2) = irgntaa( idxlmo1, idxlmo2) + conjg( pref( lm3)*expqr)*conjg( listgnt( i, lm2, lm1))*iraa( lm2l( lm3), l1, o1, l2, o2)
                i = i + 1
              end do
            end do
          end do
! APW-LO and LO-APW
          do idxlmo1 = 1, lmo
            l1 = lmo2l( idxlmo1)
            m1 = lmo2m( idxlmo1)
            o1 = lmo2o( idxlmo1)
            idxlo1 = 0
            do ilo1 = 1, nlorb( is)
              l2 = lorbl( ilo1, is)
              do m2 = -l2, l2
                lm2 = idxlm( l2, m2)
                idxlo1 = idxlo1 + 1
                ! APW-LO
                i = 1
                do while( idxgnt( i, lm2, lm1) .ne. 0)
                  lm3 = idxgnt( i, lm2, lm1)
                  irgntal( idxlmo1, idxlostart+idxlo1) = irgntal( idxlmo1, idxlostart+idxlo1) + conjg( pref( lm3)*expqr)*conjg( listgnt( i, lm2, lm1))*iral( lm2l( lm3), l1, o1, ilo1)
                  i = i + 1
                end do
                ! LO-APW
                i = 1
                do while( idxgnt( i, lm1, lm2) .ne. 0)
                  lm3 = idxgnt( i, lm1, lm2)
                  irgntla( idxlostart+idxlo1, idxlmo1) = irgntla( idxlostart+idxlo1, idxlmo1) + conjg( pref( lm3)*expqr)*conjg( listgnt( i, lm1, lm2))*iral( lm2l( lm3), l1, o1, ilo1)
                  i = i + 1
                end do
              end do
            end do
          end do
! LO-LO
          idxlo1 = 0
          do ilo1 = 1, nlorb( is)
            l1 = lorbl( ilo1, is)
            do m1 = -l1, l1
              lm1 = idxlm( l1, m1)
              idxlo1 = idxlo1 + 1
              idxlo2 = 0
              do ilo2 = 1, nlorb( is)
                l2 = lorbl( ilo2, is)
                do m2 = -l2, l2
                  lm2 = idxlm( l2, m2)
                  idxlo2 = idxlo2 + 1
                  i = 1
                  do while( idxgnt( i, lm2, lm1) .ne. 0)
                    lm3 = idxgnt( i, lm2, lm1)
                    irgntll( idxlostart+idxlo1, idxlostart+idxlo2) = irgntll( idxlostart+idxlo1, idxlostart+idxlo2) + conjg( pref( lm3)*expqr)*conjg( listgnt( i, lm2, lm1))*irll( lm2l( lm3), ilo1, ilo2)
                    i = i + 1
                  end do
                end do
              end do
            end do
          end do
! update contribution muffin tin sums
          ! APW-APW
          auxmat(:,:) = zzero
          call ZGEMM( 'N', 'N', lmo, ngk( 1, ikq), lmo, zone, irgntaa, lmo, match_combined2( 1:lmo, 1:ngk( 1, ikq)), lmo, zzero, auxmat, lmo)
          call ZGEMM( 'C', 'N', ngk( 1, ik), ngk( 1, ikq), lmo, zone, match_combined1( 1:lmo, 1:ngk( 1, ik)), lmo, auxmat, lmo, zone, mataa, ngk( 1, ik))
          ! APW-LO
          call ZGEMM( 'C', 'N', ngk( 1, ik), nlotot, lmo, zone, match_combined1( 1:lmo, 1:ngk( 1, ik)), lmo, irgntal, lmo, zone, matal, ngk( 1, ik))
          ! LO-APW
          call ZGEMM( 'N', 'N', nlotot, ngk( 1, ikq), lmo, zone, irgntla, nlotot, match_combined2( 1:lmo, 1:ngk( 1, ikq)), lmo, zone, matla, nlotot)
          ! LO-LO
          matll(:,:) = matll(:,:) + irgntll(:,:)
! update starting index for local orbitals
          idxlostart = idxlostart + idxlo1
! end loops over atoms and species
        end do
      end do
! final calculation of total muffin tin contribution
      blockmt(:,:) = zzero
      blockmt( 1:ngk( 1, ik), 1:ngk( 1, ikq)) = mataa(:,:)
      blockmt( 1:ngk( 1, ik), (ngk( 1, ikq)+1):(ngk( 1, ikq)+nlotot)) = matal(:,:)
      blockmt( (ngk( 1, ik)+1):(ngk( 1, ik)+nlotot), 1:ngk( 1, ikq)) = matla(:,:)
      blockmt( (ngk( 1, ik)+1):(ngk( 1, ik)+nlotot), (ngk( 1, ikq)+1):(ngk( 1, ikq)+nlotot)) = matll(:,:)
      if( (ik .eq. 1) .and. (ikq .eq. 2)) then
        write(*,*) '---------------------------------------------------------'
        write(*,*) shape( blockmt)
        do l1 = 1, ngk( 1, ik)+nlotot
          do l2 = 1, ngk( 1, ikq)+nlotot-1
            !write( *, '(E23.12,"+1i*(",E23.12,"), ")', advance='no') matll( l1, l2)
            write( *, '(E23.12,", ")', advance='no') abs( blockmt( l1, l2))
          end do
          !write( *, '(E23.12,"+1i*(",E23.12,"); ")', advance='no') matll( l1, l2)
          write( *, '(E23.12,"; ")') abs( blockmt( l1, l2))
        end do
        write(*,*)
      end if
      do jst = ist2, ist2 + nst2 - 1
        do ist = ist1, ist1 + nst1 - 1
          call ZGEMV( 'N', ngk( 1, ik)+nlotot, ngk( 1, ikq)+nlotot, zone, blockmt, ngk( 1, ik)+nlotot, evecfv2( 1:(ngk( 1, ikq)+nlotot), jst), 1, zzero, auxvec, 1)
          em( ist-ist1+1, jst-ist2+1) = dot_product( conjg( evecfv1( 1:(ngk( 1, ik)+nlotot), ist)), auxvec)
        end do
      end do

!--------------------------------------!
!     interstitial matrix elements     !
!--------------------------------------!
     allocate( vgkql( 3, ngkmax), vgkqc( 3, ngkmax), gkqc( ngkmax), tpgkqc( 2, ngkmax))
! generate the G+k+q-vectors
     call gengpvec (vkl( 1, ikq), vkqc, ngkq, igkqig, vgkql, vgkqc, gkqc, tpgkqc)
! store q+k-k', where k' is the k+q-vector mapped to [0,1)
      v1(:) = vecqc(:) + vkc( :, ik) - vkqc(:)
! compute exp(i(q+k-k').r) times by the characteristic function
      ir = 0
      Do i3 = 0, ngrid(3) - 1
        v2 (3) = dble( i3)/dble( ngrid(3))
        Do i2 = 0, ngrid(2) - 1
          v2 (2) = dble( i2)/dble( ngrid(2))
          Do i1 = 0, ngrid(1) - 1
            v2 (1) = dble( i1)/dble( ngrid(1))
            ir = ir + 1
            Call r3mv( input%structure%crystal%basevect, v2, v3)
            t1 = dot_product( v1(:), v3(:))
            zfir1( ir) = cfunir( ir)*cmplx( Cos(t1), Sin(t1), 8)
          End Do
        End Do
      End Do
! compute interstitial wavefunctions for k-point k
      Do jst = ist1, ist1 + nst1 - 1
        zfir2 (:) = 0.d0
        Do igk = 1, ngk( 1, ik)
          ifg = igfft( igkig( igk, 1, ik))
          zfir2( ifg) = evecfv1( igk, jst)
        End Do
! Fourier transform wavefunction to real-space
        Call zfftifc( 3, ngrid, 1, zfir2)
! multiply with the phase and characteristic function
        zfir2(:) = zfir2(:)*zfir1(:)
! Fourier transform back to G-space
        Call zfftifc (3, ngrid,-1, zfir2)
! store in wfir
        Do igk = 1, ngkq
          ifg = igfft( igkqig( igk))
          wfir( igk) = zfir2( ifg)
        End Do
! add to the first-variational matrix elements
        Do ist = ist2, ist2 + nst2 - 1
          em( jst-ist1+1, ist-ist2+1) = em( jst-ist1+1, ist-ist2+1) + zdotc( ngkq, wfir, 1, evecfv2( :, ist), 1)
        End Do
      End Do
!-------------------------------------------!
!     second-variational matrix elements    !
!-------------------------------------------!
!      If (input%groundstate%tevecsv) Then
! get the second-variational eigenvectors
!         Call getevecsv( vkl( :, ik), evecsv1)
!         Call getevecsv( vkql, evecsv2)
!         Do i = ist1, ist1 + nst1 - 1
!            Do j = ist2, ist2 + nst2 - 1
!               zsum = 0.d0
!               k = 0
!               Do ispn = 1, nspinor
!                  Do ist = 1, nstfv
!                     k = k + 1
!                     l = (ispn-1)*nstfv
!                     Do jst = 1, nstfv
!                        l = l + 1
!                        zsum = zsum + em( ist, jst)*conjg( evecsv2( k, i))*evecsv1( l, j)
!                     End Do
!                  End Do
!               End Do
!               emat( i-ist1+1, j-ist2+1) = zsum
!            End Do
!         End Do
!      Else
         emat(:,:) = em(:,:)
!      End If
      Deallocate( igkqig, jlqr)
      Deallocate( apwalm1, apwalm2, evecfv1, evecfv2)
      If( input%groundstate%tevecsv) deallocate( evecsv1, evecsv2)
      Deallocate( wfir, zfir1, zfir2, em)
      Return
End Subroutine
