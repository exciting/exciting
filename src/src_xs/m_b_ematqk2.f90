! Copyright(C) 2006-2008 S. Sagmeister and Claudia Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
module m_b_ematqk2

  implicit none

  contains

    !BOP
    ! !ROUTINE: b_ematqk2
    ! !INTERFACE:
    subroutine b_ematqk2(iq, ik, emat, bc)
    ! !USES:
      use modinput, only: input
      use mod_misc, only: task, filext
      use mod_constants, only: zzero, zone
      use mod_eigenvalue_occupancy, only: nstfv
      use mod_muffin_tin, only: idxlm
      use mod_Gkvector, only: ngk, vgkl, igkig, gkmax
      use mod_APW_LO, only: nlotot, apword, nlorb, lorbl
      use mod_Gvector, only: gc, ngvec, cfunig, ivg, ivgig
use mod_Gvector, only: intgv
      use mod_kpoint, only: vkl
      use mod_qpoint, only: vql
      use mod_atoms, only: nspecies, natmtot, natoms, idxas
      use modxs, only: bcbs, ikmapikq, msg, ngk0, xiohalo,&
                     & xiuhloa, evecfv,&
                     & vkl0, vgkl0, evecfv0, evalsv0,&
                     & occsv0, apwmaxsize, apwsize, losize,&
                     & lomaxsize, cmtfun0, cmtfun, apwcmt0,&
                     & apwcmt, ngq, igqig,&
                     & fnetim, fftmap_type, igkig0,&
                     & cpumtaa, cpumtalo, cpumtloa, cpumtlolo,&
                     & filext0, iqmt0, iqmt1
      use modxs, only: randphases
      use summations, only: doublesummation_simple_cz
      use m_getapwcmt
      use m_getlocmt
      use m_putemat
      use m_emattim
      use m_getunit
      use m_genfilname
#ifdef USEOMP
      use omp_lib
#endif
    ! !DESCRIPTION:
    ! Calculates plane wave elements between state ranges defined
    ! in (bc%n1,bc%il1,bc%iu1) and (bc%n2,bc%il2,bc%iu2) and saves them in 
    ! emat.
    !
    ! !REVISION HISTORY:
    !   Added to documentation scheme. (Aurich)
    !   Substituted the modxs:xiou reference with the pointer modxs:emat that
    !   can point to any of the modxs:xiXY.
    !EOP
    !BOC

      implicit none
          
      ! Arguments
      integer, intent(in) :: iq, ik
      type(bcbs), intent(in) :: bc
      complex(8), intent(inout) :: emat(:,:,:)

      ! Local variables
      character(*), parameter :: thisnam = 'b_ematqk'
      ! Allocatable arrays
      complex(8), allocatable :: evecfvo0(:, :)
      complex(8), allocatable :: evecfvu(:, :)
      complex(8), allocatable :: evecfvo20(:, :)
      complex(8), allocatable :: evecfvu2(:, :)
      complex(8), allocatable :: zfft0(:, :), zfft(:), zfftres(:), zfftcf(:)
      complex(8), allocatable :: xihir(:, :)
      ! Expansion coefficients of apw and lo functions
      complex(8), allocatable :: integrals(:, :, :)
      integer :: ikq, igq, n, n0, is, l, m, io, naug, ia, ias, lm, ilo
      integer :: whichthread, ig, igk, ix, igs, ist2, ist1
      real(8) :: cpuini, cpuread, cpumain, cpuwrite, cpuall
      real(8) :: cpugnt, cpumt, cpuir, cpufft
      real(8) :: cpumalores, cpumloares
      real(8) :: cpumlolores, cpumirres, cpudbg
      real(8) :: cpu0, cpu1, cpu00, cpu01, cpugntlocal, cpumtlocal
      real(8) :: vkql(3), vkkpq(3)
      integer :: shift(3), iv(3)
      type(fftmap_type) :: fftmap
      real(8) :: emat_gmax
      character(256) :: filename
      real(8), parameter :: epslat = 1.0d-6

      integer :: ist, i
      logical :: shiftcheck

      ! If task 330 is 'writeemat'
      if(task .eq. 330) then
        call chkpt(3, (/ task, iq, ik /),&
          & 'b_ematqk2: task, q - point index, k - point index; q - dependent matrix elements')
      end if

      call timesec(cpu0)

      ! Find k'=-k-q-point
      ikq = ikmapikq(ik, iq)
      !write(*,*) "ematqk2: ik=",ik," iq=",iq,"ikq=",ikq

      ! Check for stop statement
      write(msg, *) 'for q-point', iq, ': k-point:', ik - 1, ' finished'
      !call xschkstop

      ! Timing variables
      cpumtaa = 0.d0
      cpumtalo = 0.d0
      cpumtloa = 0.d0
      cpumtlolo = 0.d0
      cpugnt = 0.d0
      cpumt = 0.d0
      cpuir = 0.d0
      cpumalores = 0.d0
      cpumloares = 0.d0
      cpumlolores = 0.d0
      cpumirres = 0.d0
      cpudbg = 0.d0

      ! Get number of G+k and G+k' vectors
      n0 = ngk0(1, ik)
      n = ngk(1, ikq)
      !write(*,*) "ematqk2: n0, n", n0, n

      ! Allocate matrix elements array
      if(allocated(xiohalo)) deallocate(xiohalo)
      allocate(xiohalo(bc%n1, nlotot))
      if(allocated(xiuhloa)) deallocate(xiuhloa)
      allocate(xiuhloa(nlotot, bc%n2))

      ! Allocate temporary arrays
      allocate(evecfvo0(nlotot, bc%n1))
      allocate(evecfvu(nlotot, bc%n2))

      ! Zero arrays
      xiohalo(:, :) = zzero
      xiuhloa(:, :) = zzero

      call timesec(cpu1)
      cpuini = cpu1 - cpu0

      ! Read eigenvectors for k'
      !   Read first variational eigenvectors from EVECFV_QMTXXX.OUT 
      !   (file extension needs to be set by calling routine)
      call getevecfv(vkl(1, ikq), vgkl(1, 1, 1, ikq), evecfv)
      !if(task == 440 .or. task == 441) then 
      !  do ist = 1, nstfv
      !    evecfv(:,ist,:) = evecfv(:,ist,:)*randphases(ist,ikq)
      !  end do
      !end if

      ! Save local orbital coefficients
      evecfvu(:, :) = evecfv(ngk(1, ikq)+1:ngk(1, ikq)+nlotot, bc%il2:bc%iu2, 1)

      ! Read eigenvectors for k
      !   Read first variational eigenvectors from EVECFV_QMT000.OUT 
      call getevecfv0(vkl0(1, ik), vgkl0(1, 1, 1, ik), evecfv0)
      !if(task == 440 .or. task == 441) then 
      !  do ist = 1, nstfv
      !    evecfv0(:,ist,:) = evecfv0(:,ist,:)*randphases(ist,ik)
      !  end do
      !end if

      ! Save local orbital coefficients
      evecfvo0(:, :) = evecfv0(ngk0(1, ik)+1:ngk0(1, ik)+nlotot, bc%il1:bc%iu1, 1)

      ! Determine maximum number of LAPW basis
      ! functions over all species used in the 
      ! xs gs runs.
      apwmaxsize=0
      allocate(apwsize(nspecies))
      do is=1, nspecies
        naug=0
        ! Note lmaxapwwf is set to xs%lmaxmat = 5 by default
        do l=0, input%xs%lmaxapwwf
          naug=naug+(2*l+1)*apword(l, is)
        end do
        apwsize(is)=naug
        apwmaxsize=max(apwmaxsize, naug)
      end do

      ! Determine maximum number of LO basis functions
      ! over all species. 
      allocate(losize(nspecies))
      lomaxsize=0
      losize=0
      do is=1, nspecies
        do ilo = 1, nlorb(is)
          losize(is)=losize(is)+2*lorbl(ilo, is)+1
        end do
        lomaxsize=max(lomaxsize, losize(is))
      end do

      allocate(cmtfun0(bc%n1, apwmaxsize+lomaxsize, natmtot))
      cmtfun0=zzero
      allocate(cmtfun(bc%n2, apwmaxsize+lomaxsize, natmtot))
      cmtfun=zzero
      apwcmt0=zzero
      apwcmt=zzero

      ! Get C_{m,G+k}*A^{l,m,p,a}_{G+k}, i.e. the eigenvector coefficients
      ! corresponing to the APW part of the basis multiplied by the matching
      ! coefficients. 
      filename = 'APWCMT'//trim(adjustl(filext0))
      call getapwcmt(iqmt0, ik, 1, nstfv, input%xs%lmaxapwwf, apwcmt0,&
        & tbra=.true., fname=trim(adjustl(filename)))

      ! Repack (C*A)^* for APW basis function and C^* for LO basis functions
      ! into one array 
      ilo=0
      ! Loop over atoms
      do is=1, nspecies
        do ia=1, natoms(is)
          ias = idxas(ia, is)
          ! naug = combinded l,m,p index for atom ias
          naug=0
          ! l,m combinded index
          do l=0, input%xs%lmaxapwwf
            do m=-l, l
              lm=idxlm(l, m)
              ! APW order p
              do io=1, apword(l, is)
                naug=naug+1
                ! Repack (C*A)^*
                cmtfun0(1:bc%n1, naug, ias)=conjg(apwcmt0(bc%il1:bc%iu1, io, lm, ias))
              end do
            end do
          end do
          ! Append LO eigenvector coefficients (C)^* for that atom
          cmtfun0(1:bc%n1, naug+1:naug+losize(is), ias)=&
            & conjg(transpose(evecfvo0(ilo+1:ilo+losize(is), 1:bc%n1)))
          ilo=ilo+losize(is)
        end do
      end do

      ! Get C_{m,G+k'}*A^{l,m,p,a}_{G+k'} 
      filename = 'APWCMT'//trim(adjustl(filext))
      call getapwcmt(iqmt1, ikq, 1, nstfv, input%xs%lmaxapwwf, apwcmt,&
        & tbra=.false., fname=trim(adjustl(filename)))

      ! Repack (C*A)^* for APW basis function and C^* for LO basis functions
      ! into one array as above 
      ilo=0
      do is=1, nspecies
        do ia=1, natoms(is)
          ias = idxas(ia, is)
          naug=0
          do l=0, input%xs%lmaxapwwf
            do m=-l, l
              do io=1, apword(l, is)
                lm=idxlm(l, m)
                naug=naug+1
                cmtfun(1:bc%n2, naug, ias)=conjg(apwcmt(bc%il2:bc%iu2, io, lm, ias))
              end do
            end do
          end do
          cmtfun(1:bc%n2, naug+1:naug+losize(is), ias)=&
            & conjg(transpose(evecfvu(ilo+1:ilo+losize(is), 1:bc%n2)))
          ilo=ilo+losize(is)
        end do
      end do

      call timesec(cpu0)
      cpuread = cpu0 - cpu1

      ! Zero matrix elements array
      emat(:, :, :) = zzero

      whichthread=0

      ! Loop over G+q vectors
      cpugntlocal=0.0d0
      cpumtlocal=0.0d0

#ifdef USEOMP
    !$omp parallel default(shared) private(igq, integrals, cpu00, cpu01, whichthread)
#endif
      ! Allocate radial integrals
      allocate(integrals(apwmaxsize+lomaxsize, apwmaxsize+lomaxsize, natmtot))
#ifdef USEOMP
      whichthread=omp_get_thread_num()
    !$omp do
#endif
      do igq = 1, ngq(iq)
        call timesec(cpu00)
        ! Summation of gaunt coefficients w.r.t. radial integrals
        call b_ematgntsum2(iq, igq, integrals)

        call timesec(cpu01)
        if(whichthread.eq.0) cpugnt = cpugnt + cpu01 - cpu00

        ! Muffin-tin contribution
        call b_ematqkgmt2(iq, ik, igq, integrals, emat(:,:,igq), bc)
        call timesec(cpu00)

        if(whichthread.eq.0) cpumt = cpumt + cpu00 - cpu01
      end do ! igq
#ifdef USEOMP
    !$omp end do
#endif
      deallocate(integrals)        

#ifdef USEOMP
    !$omp end parallel
#endif
      deallocate(apwsize, losize, cmtfun, cmtfun0)
      deallocate(evecfvu, evecfvo0)
     
      allocate(evecfvo20(n0, bc%n1))
      allocate(evecfvu2(n, bc%n2))

      ! Save plane wave coefficients of ket state
      evecfvu2(:, :) = evecfv(1:ngk(1, ikq), bc%il2:bc%iu2, 1)
      ! Save plane wave coefficients of bra state
      evecfvo20(:, :) = evecfv0(1:ngk0(1, ik), bc%il1:bc%iu1, 1)
     
      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
      ! Interstitial contribution                               !
      !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

      ! Fourier transforms

      ! What is done here:
      ! Consider the integral over the interstitial region of <mk|e^{-i(G+q)r}|(nk)^*>.
      ! Writing it using the Lattice theta function as an integral over
      ! the whole volume: \int_V \theta(r) \Phi^*_{mk}(r) e^{-i(G+q)r} \Phi^*_{nk}(r)
      !
      ! \int_V \theta(r) \Phi^*_{mk}(r) e^{-i(G+q)r} \Phi^*_{nk}(r) 
      ! Insert lattice FT expressions -->
      ! =\int_V (\Sum_{G1} \theta(G1) e^{iG1r}) (\Sum_{G2} \Phi_{mk}(G2) e^{i(G2+k)r})^*
      !        e^{-i(G+q)r}(\Sum_{G3} \Phi_{nk}(G3) e^{i(G3+k')r})^*
      ! =\int_V (\Sum_{G1} \theta(G1) e^{iG1r}) (\Sum_{G2} \Phi_{mk}(G2) e^{iG2r})^*
      !         (\Sum_{G3} \Phi_{nk}(G3) e^{iG3)r})^* e^{i(-k'-k-q)r} e^{-iGr}
      ! =\int_V (\Sum_{G1} \theta(G1) e^{iG1r}) (\Sum_{G2} \Phi_{mk}(G2) e^{iG2r})^*
      !         (\Sum_{G3} \Phi_{nk}(G3) e^{iG3)r})^* e^{-i(G-Gshift)r}
      !
      ! By design it is e^{i(-k'-k-q)r} = e^{i Gshift r} with some Gshift lattice vector.
      ! All sums in the round brackets are directly computable as
      ! inverse FFT's using the eigenvector coefficients (LAPW part) and \theta(G).
      ! The corresponding real space quantities are then multiplied and subsequently
      ! transformed back to reciprocal space.
      ! If Gshift is non zero the final back transform needs to be shifted.
      ! 
      ! So 3 inverse FFT's are computed, the real space results multiplied and
      ! subsequently back transformed with another FFT to get the integral at all G

      call timesec(cpu00)
   
      ! Index link between k and k'=-(k+q) grids
      ikq = ikmapikq(ik, iq)

      ! Determine G vector that maps -k'-k-q back to [0,1)
      vkkpq(:) = -vkl(:,ikq)-vkl0(:,ik)-vql(:,iq)
      call r3frac(epslat, vkkpq, shift)
      if(any(abs(vkkpq)>epslat)) then 
        write(*,*) "Error: -k'-k-q does not map back to Gamma"
        write(*,*) "ikq, vkl", ikq, vkl(:,ikq)
        write(*,*) "ik, vkl0", ik, vkl0(:,ik)
        write(*,*) "iq, vql", iq, vql(:,iq)
        call terminate
      end if
        
      emat_gmax=2*gkmax +input%xs%gqmax

      call genfftmap(fftmap, emat_gmax)
      allocate(zfft0(fftmap%ngrtot+1, bc%n1))
      zfft0=zzero

      allocate(zfftcf(fftmap%ngrtot+1))
      zfftcf=zzero

      ! Inverse FT of characteristic lattice function (0 inside MT, 1 outside)
      ! form reciprocal to real space, incorporate Gshift here.
      do ig=1, ngvec
        if(gc(ig).lt.emat_gmax) then
          zfftcf(fftmap%igfft(ig))=cfunig(ig)
        end if
      end do
      call zfftifc(3, fftmap%ngrid, 1, zfftcf)

#ifdef USEOMP
  !$omp parallel default(shared) private(ist1, igk)
  !$omp do
#endif
      ! Inverse FT of interstitial part of bra state
      ! form reciprocal to real space for each band
      do ist1=1, bc%n1

        ! Map plane wave coefficients from G+k onto G fftmap 
        do igk=1, ngk0(1, ik)
          zfft0(fftmap%igfft(igkig0(igk, 1, ik)), ist1)=evecfvo20(igk, ist1)
        end do
        ! Do the iFT
        call zfftifc(3, fftmap%ngrid, 1, zfft0(:, ist1))
        ! Conjugate and multiply by characteristic lattice function
        zfft0(:, ist1)=conjg(zfft0(:, ist1))*zfftcf(:) 

      end do
#ifdef USEOMP
  !$omp end do
  !$omp end parallel
#endif
     
#ifdef USEOMP
  !$omp parallel default(shared) private(ist1, ist2, igk, shiftcheck, i, iv, igs, zfft, zfftres, igq)
#endif
      allocate(zfftres(fftmap%ngrtot+1))
      allocate(zfft(fftmap%ngrtot+1))
#ifdef USEOMP
  !$omp do
#endif
      ! For each band used for the ket state do 
      do ist2=1, bc%n2

        zfft=zzero

        ! Map plane wave coefficients form G+k' onto G fftmap
        do igk=1, ngk(1, ikq)
          ! Map to fftmap G+k' --> G
          zfft(fftmap%igfft(igkig(igk, 1, ikq)))=evecfvu2(igk, ist2)
        end do
        ! Do the iFT of the ket state
        call zfftifc(3, fftmap%ngrid, 1, zfft)

        ! For all bra state bands do
        do ist1=1, bc%n1

          ! Take the real space product:
          ! eveck^*_ist1(r)*\Theta(r)*eveckq^*_ist2(r)
          do ig=1, fftmap%ngrtot
            zfftres(ig)=zfft0(ig, ist1)*conjg(zfft(ig)) 
          end do

          ! Back transform the result to G space
          call zfftifc(3, fftmap%ngrid, -1, zfftres)

          ! Add the interstitial result to emat for all G+q within the cutoff
          do igq=1, ngq(iq)
            if(any(shift /= 0)) then 
              iv = ivg(:,igqig(igq,iq)) - shift
              shiftcheck=.true.
              do i=1,3
                shiftcheck = (iv(i) >= intgv(i,1) .and. iv(i) <= intgv(i,2) .and. shiftcheck)
              end do
              if(shiftcheck) then 
                igs = ivgig(iv(1),iv(2),iv(3))
                if(gc(igs) .lt. emat_gmax) then 
                  emat(ist1, ist2, igq)=emat(ist1, ist2, igq)&
                    &+ zfftres(fftmap%igfft(igs))
                end if
              end if
            else
              emat(ist1, ist2, igq)=emat(ist1, ist2, igq)&
                &+ zfftres(fftmap%igfft(igqig(igq, iq)))
            endif
          end do

        end do
      end do
#ifdef USEOMP
  !$omp end do
#endif
      deallocate(zfftres, zfft)
#ifdef USEOMP
!$omp end parallel
#endif

      deallocate(fftmap%igfft)
      deallocate(zfft0, zfftcf)
      call timesec(cpu01)
      cpufft=cpu01-cpu00

      !write(*,*) "ematqk: emat2 = "
      !write(*,'(2E12.5)') emat

      call timesec(cpu1)
      cpumain = cpu1 - cpu0

      call timesec(cpu0)
      cpuwrite = cpu0 - cpu1
      cpuall = cpuini + cpuread + cpumain + cpuwrite

      ! Write timing information
      ! Tasks: 440=scrcoulint, 441=exccoulint, 450=kernxc_bse, 451=kernxc_bse3
      if( (task.ne.440) .and. (task .ne. 441) .and. (task .ne. 450)&
        & .and. (task .ne. 451)) then
        call emattim(iq, ik, trim(fnetim), cpuini, cpuread, cpumain,&
          & cpuwrite, cpuall, cpugnt, cpumt, cpuir, cpumalores,&
          & cpumloares, cpumlolores, cpumirres, cpudbg, cpumtaa,&
          & cpumtalo, cpumtloa, cpumtlolo, cpufft)
      end if

    end subroutine b_ematqk2
    !EOC

    subroutine b_ematqkgmt2(iq, ik, igq, integrals, emat, bc)
      use modinput, only: input
      use mod_constants, only: zzero, zone, fourpi
      use mod_atoms, only: natmtot, nspecies, natoms, idxas
      use modxs, only: bcbs, apwmaxsize, lomaxsize, &
                     & sfacgq,&
                     & apwsize, losize, cmtfun, cmtfun0,&
                     & cpumtaa
#ifdef USEOMP
    use omp_lib
#endif

      implicit none

      ! Arguments
      integer, intent(in) :: iq, ik, igq
      type(bcbs), intent(in) :: bc
      complex(8), intent(inout) :: emat(:,:)

      complex(8) :: integrals(apwmaxsize+lomaxsize,apwmaxsize+lomaxsize,natmtot)

      ! Local variables
      character(*), parameter :: thisnam = 'b_ematqkgmt'
      integer :: is, ia, ias
      integer :: lmax1, lmax3, ikt, zmsize, whichthread
      complex(8), allocatable :: zm(:,:)
      complex(8) :: prefactor
      real(8) :: cmt0, cmt1

#ifdef USEOMP
      whichthread=omp_get_thread_num()
#else
      whichthread=0
#endif

      ikt = ik
      lmax1 = input%xs%lmaxapwwf
      lmax3 = lmax1
      zmsize=apwmaxsize+lomaxsize

      allocate(zm(1:bc%iu2-bc%il2+1,zmsize))

      emat(:, :) = zzero

      ! Loop over species and atoms
      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas(ia, is)

          call timesec(cmt0)
          !---------------------------!
          !     apw-apw contribution  !
          !---------------------------!
          prefactor=fourpi*conjg(sfacgq(igq, ias, iq))

          call zgemm('n', 'n', bc%n2, apwsize(is)+losize(is), apwsize(is)+losize(is),&
            & zone, cmtfun(1,1,ias), bc%n2, integrals(1,1,ias), apwmaxsize+lomaxsize,&
            & zzero, zm, bc%n2)
          call zgemm('n', 't', bc%n1, bc%n2, apwsize(is)+losize(is),&
            & prefactor, cmtfun0(1,1,ias), bc%n1, zm, bc%n2, zone, emat, bc%n1)

          call timesec(cmt1)
          if(whichthread.eq.0) then
            cpumtaa = cpumtaa + cmt1 - cmt0
          endif

        ! End loop over species and atoms
        end do ! ia
      end do ! is

      deallocate(zm)
    end subroutine b_ematqkgmt2

    subroutine b_ematgntsum2(iq, igq, integrals)
      use modinput, only: input
      use mod_misc, only: filext
      use mod_constants, only: zzero, zil
      use mod_muffin_tin, only: idxlm
      use mod_atoms, only: natmtot, nspecies, natoms, idxas
      use mod_APW_LO, only: lolmax, apwordmax, nlomax, apword,&
                          & nlorb, lorbl
      use modxs, only: apwmaxsize, apwsize, lomaxsize, ylmgq,&
                     & xsgnt, rilolo, riaa, riloa
      use m_findgntn0, only: l1map, l2map, l3map,&
                           & m1map, m2map, m3map,&
                           & l1shape, l2shape, l3shape,&
                           & m1shape, m2shape, m3shape
      use m_getunit

      implicit none

      ! Arguments
      integer, intent(in) :: iq, igq
      complex(8) :: integrals(apwmaxsize+lomaxsize, apwmaxsize+lomaxsize, natmtot)

      ! Local variables
      integer :: is, ia, ias, iaug1, iaug2
      integer :: l1, l2, l3, m2, lm2
      integer :: ilo, ilo1, ilo2, io, io1, io2
      integer :: lmax1, lmax2, lmax3, lmmax1, lmmax2, lmmax3
      integer :: u1, u2, u3, u4
      integer :: m1, m3, lm1, lm3, cl1, cm1, cl2, cm2, cl3, cm3
      complex(8), dimension(:,:,:,:), allocatable :: intrgaa, intrgloa, intrglolo, intrgalo

      ! Set lm related local variables
      lmax1 = max(input%xs%lmaxapwwf, lolmax)
      lmax2 = input%xs%lmaxemat
      ! lmax1 and lmax3 should be the same
      lmax3 = lmax1
      lmmax1 = (lmax1+1) ** 2
      lmmax2 = (lmax2+1) ** 2
      lmmax3 = (lmax3+1) ** 2

      !write(*,*) "lmax1, lmax2, lmax3", lmax1, lmax2, lmax3
      !write(*,*) "lmmax1, lmmax2, lmmax3", lmmax1, lmmax2, lmmax3

      ! Allocate arrays for radial integrals and Bessel functions
      allocate(intrgaa(lmmax1, apwordmax, lmmax3, apwordmax))
      allocate(intrgloa(-lolmax:lolmax, nlomax, lmmax3, apwordmax))
      allocate(intrgalo(-lolmax:lolmax, nlomax, lmmax3, apwordmax))
      allocate(intrglolo(-lolmax:lolmax, nlomax, -lolmax:lolmax, nlomax))

      integrals=zzero

      ! Debug output
      if(input%xs%dbglev .gt. 2) then
        ! Apw-apw
        call getunit(u1)
        open(unit=u1, file='IRADGAUNTaa'//filext, form='formatted', action='write', status='replace')
        write(u1, '(a)') 'igq, ias, lm1, io1, lm3, io2,   intrgaa'
        write(u1, '(a)') '------------------------------------------------------'
        ! Apw-lo
        call getunit(u2)
        open(unit=u2, file='IRADGAUNTalo'//filext, form='formatted', action='write', status='replace')
        write(u2, '(a)') 'igq, ias, m3, ilo, lm1, io,     intrgalo'
        write(u2, '(a)') '------------------------------------------------------'
        ! Lo-apw
        call getunit(u3)
        open(unit=u3, file='IRADGAUNTloa'//filext, form='formatted', action='write', status='replace')
        write(u3, '(a)') 'igq, ias, m1, ilo, lm3, io,     intrgloa'
        write(u3, '(a)') '------------------------------------------------------'
        ! Lo-lo
        call getunit(u4)
        open(unit=u4, file='IRADGAUNTlolo'//filext, form='formatted', action='write', status='replace')
        write(u4, '(a)') 'igq, ias, m1, ilo1, m3, ilo2,   intrglolo'
        write(u4, '(a)') '------------------------------------------------------'
      end if

      ! Begin loop over species
      do is = 1, nspecies
        ! Begin loop over atoms
        do ia = 1, natoms(is)

          intrgaa(:, :, :, :) = zzero
          intrgloa(:, :, :, :) = zzero
          intrgalo(:, :, :, :) = zzero
          intrglolo(:, :, :, :) = zzero

          ias = idxas(ia, is)
          !---------------------------!
          !     apw-apw integrals     !
          !---------------------------!
          ! Loop over l,m,p combinations of bra apw
          ! for which there are non zero gaunt coefficients
          do cl1 = 1, l1shape
            l1 = l1map(cl1)
            do cm1 = 1, m1shape(l1)
              m1 = m1map(l1, cm1)
              lm1 = idxlm(l1, m1)
              do io1 = 1, apword(l1, is)

                ! Loop over l',m',p' combinations of ket apw
                ! for which there are non-zero gaunt coefficients 
                ! given l,m of the bra apw state
                do cl2 = 1, l2shape(l1, m1)
                  l3 = l2map(l1, m1, cl2)
                  do cm2 = 1, m2shape(l1, m1, l3)
                    m3 = m2map(l1, m1, l3, cm2)
                    lm3 = idxlm(l3, m3)
                    do io2 = 1, apword(l3, is)

                      ! Loop over l'',m'' combinations of the 
                      ! Reighley expansion of e^{-i(G+q)r}, for which
                      ! there are non-zero gaunt coefficients given l,m and l',m'
                      do cl3 = 1, l3shape(l1, m1, l3, m3)
                        l2 = l3map(l1, m1, l3, m3, cl3)
                        do cm3 = 1, m3shape(l1, m1, l3, m3, l2)
                          m2 = m3map(l1, m1, l3, m3, l2, cm3)
                          lm2 = idxlm(l2, m2)

                          intrgaa(lm1, io1, lm3, io2) =&
                            & intrgaa(lm1, io1, lm3, io2)&
                            &+ conjg(zil(l2)) * riaa(l1, io1, l3, io2, l2, ias, igq)&
                            &* conjg(ylmgq(lm2, igq, iq)) * xsgnt(lm1, lm2, lm3)

                        end do
                      end do

                      ! Debug output
                      if(input%xs%dbglev .gt. 2) then
                         write(u1, '(6i5, 2g18.10)') igq, ias, &
                           & lm1, io1, lm3, io2, intrgaa(lm1, io1, lm3, io2)
                      end if

                    end do
                  end do
                end do

              end do
            end do
          end do

          ! Write APW-APW part in total intrgals array. 
          !   Note: bra and ket index swapped.
          ! Changed for calculation of \int \Phi^*_mk e^{-i(G+q)r} \Phi^*_nk'

          ! Combined l',m',p' index 
          iaug2=0
          do l3=0, input%xs%lmaxapwwf
            do m3=-l3, l3
              ! Change 1: m' --> -m'
              lm3=idxlm(l3, -m3)
              do io2=1, apword(l3, is)
                iaug2=iaug2+1

                ! Combined l,m,p index
                iaug1=0
                do l1=0, input%xs%lmaxapwwf
                  do m1=-l1, l1
                      lm1=idxlm(l1, m1)
                      do io1=1, apword(l1, is)
                        iaug1=iaug1+1

                        ! Change 2: *(-1)^m'
                        integrals(iaug2, iaug1, ias) =&
                          & intrgaa(lm1, io1, lm3, io2)*(-1.0d0)**m3

                      end do
                    end do
                  end do

              end do
            end do
          end do

          !-------------------------------------!
          !     apw-local-orbital integrals     !
          !-------------------------------------!
          ! Loop over l,m,p combinations of bra apw
          ! for which there are non zero gaunt coefficients
          do cl1 = 1, l1shape
            l1 = l1map(cl1)
            do cm1 = 1, m1shape(l1)
              m1 = m1map(l1, cm1)
              lm1 = idxlm(l1, m1)
              do io = 1, apword(l1, is)

                ! Loop over all local orbitals for that species
                ! with l',m' combinations for which there are non-zero
                ! gaunt coefficients given the l,m of the apw
                do ilo = 1, nlorb(is)
                  l3 = lorbl(ilo, is)
                  do cm2 = 1, m2shape(l1, m1, l3)
                    m3 = m2map(l1, m1, l3, cm2)
                    lm3 = idxlm(l3, m3)

                    ! Loop over l'',m'' combinations of the 
                    ! Reighley expansion of e^{-i(G+q)r}, for which
                    ! there are non-zero gaunt coefficients given l,m and l',m'
                    do cl3 = 1, l3shape(l1, m1, l3, m3)
                      l2 = l3map(l1, m1, l3, m3, cl3)
                      do cm3 = 1, m3shape(l1, m1, l3, m3, l2)
                        m2 = m3map(l1, m1, l3, m3, l2, cm3)
                        lm2 = idxlm(l2, m2)

                        intrgalo(m3, ilo, lm1, io) = intrgalo(m3, ilo, lm1, io)&
                          &+ conjg(zil(l2)) * riloa(ilo, l1, io, l2, ias, igq)&
                          &* conjg(ylmgq(lm2, igq, iq)) * xsgnt(lm1, lm2, lm3)

                      end do
                    end do

                    ! Debug output
                    if(input%xs%dbglev .gt. 2) then
                       write(u2, '(6i5, 2g18.10)') igq, ias, &
                         & m3, ilo, lm1, io, intrgalo(m3, ilo, lm1, io)
                    end if

                  end do
                end do

              end do
            end do
          end do

          ! Write APW-LO part in total intrgals array. 
          !   Note: bra and ket index swapped.
          ! Changed for calculation of \int \Phi^*_mk e^{-i(G+q)r} \Phi^*_nk'

          ! Combined l',m' index (LO)
          iaug2=0
          do ilo = 1, nlorb(is)
            l3 = lorbl(ilo, is)
            do m3=-l3, l3
              ! Change 1: m' --> -m'
              lm3=idxlm(l3, -m3)
              iaug2=iaug2+1
             
              ! Combined l,m,p index (APW)
              iaug1=0
              do l1=0, input%xs%lmaxapwwf
                do m1=-l1, l1
                  lm1=idxlm(l1, m1)
                  do io=1, apword(l1, is)
                    iaug1=iaug1+1

                    ! Change 2: *(-1)^m'
                    integrals(apwsize(is)+iaug2, iaug1, ias) =&
                      & intrgalo(-m3, ilo, lm1, io)*(-1.0d0)**m3

                  end do
                end do
               end do

            end do
          end do

          !-------------------------------------!
          !     local-orbital-apw integrals     !
          !-------------------------------------!
          ! Loop over all local orbitals for that species (bra)
          ! with l,m combinations for which there are non zero gaunt coefficients
          do ilo = 1, nlorb(is)
            l1 = lorbl(ilo, is)
            do cm1 = 1, m1shape(l1)
              m1 = m1map(l1, cm1)
              lm1 = idxlm(l1, m1)

              ! Loop over l',m',p' combinations of (ket) apw
              ! for which there are non zero gaunt coefficients
              ! given l,m of the LO
              do cl2 = 1, l2shape(l1, m1)
                l3 = l2map(l1, m1, cl2)
                do cm2 = 1, m2shape(l1, m1, l3)
                  m3 = m2map(l1, m1, l3, cm2)
                  lm3 = idxlm(l3, m3)
                  do io = 1, apword(l3, is)

                    ! Loop over l'',m'' combinations of the 
                    ! Reighley expansion of e^{-i(G+q)r}, for which
                    ! there are non-zero gaunt coefficients given l,m and l',m'
                    do cl3 = 1, l3shape(l1, m1, l3, m3)
                      l2 = l3map(l1, m1, l3, m3, cl3)
                      do cm3 = 1, m3shape(l1, m1, l3, m3, l2)
                        m2 = m3map(l1, m1, l3, m3, l2, cm3)
                        lm2 = idxlm(l2, m2)

                        intrgloa(m1, ilo, lm3, io) = intrgloa(m1, ilo, lm3, io)&
                          &+ conjg(zil(l2)) * riloa(ilo, l3, io, l2, ias, igq)&
                          &* conjg(ylmgq(lm2, igq, iq)) * xsgnt(lm1, lm2, lm3)

                      end do
                    end do

                    ! Debug output
                    if(input%xs%dbglev .gt. 2) then
                       write(u3, '(6i5, 2g18.10)') igq, ias,&
                         & m1, ilo, lm3, io, intrgloa(m1, ilo, lm3, io)
                    end if

                  end do
                end do
              end do

            end do
          end do

          ! Write LO-APW part in total intrgals array. 
          !   Note: bra and ket index swapped.
          ! Changed for calculation of \int \Phi^*_mk e^{-i(G+q)r} \Phi^*_nk'

          ! Combined l,m index (LO)
          iaug2=0
          do ilo = 1, nlorb(is)
            l3 = lorbl(ilo, is)
            do m3=-l3, l3
              lm3=idxlm(l3, m3)
              iaug2=iaug2+1

              ! Combined l',m',p' index (APW)
              iaug1=0
              do l1=0, input%xs%lmaxapwwf
                do m1=-l1, l1
                  ! Change 1: m' --> -m'
                  lm1=idxlm(l1, -m1)
                  do io=1, apword(l1, is)
                    iaug1=iaug1+1

                    integrals(iaug1, apwsize(is)+iaug2, ias) =&
                      & intrgloa(m3, ilo, lm1, io)*(-1.0d0)**m1

                  end do
                end do
               end do

            end do
          end do

          !-----------------------------------------------!
          !     local-orbital-local-orbital integrals     !
          !-----------------------------------------------!
          ! Loop over all local orbitals for that species (bra)
          ! with l,m combinations for which there are non zero gaunt coefficients
          do ilo1 = 1, nlorb(is)
            l1 = lorbl(ilo1, is)
            do cm1 = 1, m1shape(l1)
              m1 = m1map(l1, cm1)
              lm1 = idxlm(l1, m1)

              ! Loop over all local orbitals for that species (ket)
              ! with l',m' combinations for which there are non zero gaunt coefficients
              ! given l,m of the other LO
              do ilo2 = 1, nlorb(is)
                l3 = lorbl(ilo2, is)
                do cm2 = 1, m2shape(l1, m1, l3)
                  m3 = m2map(l1, m1, l3, cm2)
                  lm3 = idxlm(l3, m3)

                  ! Loop over l'',m'' combinations of the 
                  ! Reighley expansion of e^{-i(G+q)r}, for which
                  ! there are non-zero gaunt coefficients given l,m and l',m'
                  do cl3 = 1, l3shape(l1, m1, l3, m3)
                    l2 = l3map(l1, m1, l3, m3, cl3)
                    do cm3 = 1, m3shape(l1, m1, l3, m3, l2)
                      m2 = m3map(l1, m1, l3, m3, l2, cm3)
                      lm2 = idxlm(l2, m2)

                      intrglolo(m1, ilo1, m3, ilo2) = intrglolo(m1, ilo1, m3, ilo2)&
                        &+ conjg(zil(l2)) * rilolo(ilo1, ilo2, l2, ias, igq)&
                        &* conjg(ylmgq(lm2, igq, iq)) * xsgnt(lm1, lm2, lm3)

                    end do
                  end do

                  ! Debug output
                  if(input%xs%dbglev .gt. 2) then
                     write(u4, '(6i5, 2g18.10)') igq, ias, m1, &
                       & ilo1, m3, ilo2, intrglolo(m1, ilo1, m3, ilo2)
                  end if

                end do
              end do
            end do
          end do

          ! Write LO-LO part in total intrgals array. 
          !   Note: bra and ket index swapped.
          ! Changed for calculation of \int \Phi^*_mk e^{-i(G+q)r} \Phi^*_nk'

          ! Combined l',m' index (LO)
          iaug2=0
          do ilo2 = 1, nlorb(is)
            l3 = lorbl(ilo2, is)
            do m3=-l3, l3
              iaug2=iaug2+1

              ! Combined l,m index (LO)
              iaug1=0
              do ilo1 = 1, nlorb(is)
                l1 = lorbl(ilo1, is)
                do m1=-l1, l1
                  iaug1=iaug1+1

                  ! Change: m'-->-m' & *(-1)^m'
                  integrals(apwsize(is)+iaug2, apwsize(is)+iaug1, ias)=&
                    & intrglolo(m1, ilo1, -m3, ilo2)*(-1.0d0)**m3

                end do
              end do

            end do
          end do

        ! End loops over atoms and species
        end do
      end do

      ! Deallocate
      if(allocated(intrgaa)) deallocate(intrgaa)
      if(allocated(intrgloa)) deallocate(intrgloa)
      if(allocated(intrgalo)) deallocate(intrgalo)
      if(allocated(intrglolo)) deallocate(intrglolo)

      ! Debug code
      if(input%xs%dbglev .gt. 2) then
        ! Close files
        close(u1)
        close(u2)
        close(u3)
        close(u4)
      end if
          
    end subroutine b_ematgntsum2

end module m_b_ematqk2
