! Copyright(C) 2006-2008 S. Sagmeister and Claudia Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
module m_b_ematqk

  implicit none

  contains

    !BOP
    ! !ROUTINE: b_ematqk
    ! !INTERFACE:
    subroutine b_ematqk(iq, ik, emat, bc)
    ! !USES:
      use modinput, only: input
      use mod_misc, only: task
      use mod_constants, only: zzero, zone
      use mod_eigenvalue_occupancy, only: evalsv, occsv, nstfv
      use mod_muffin_tin, only: idxlm
      use mod_Gkvector, only: ngk, vgkl, igkig, gkmax
      use mod_APW_LO, only: nlotot, apword, nlorb, lorbl
      use mod_Gvector, only: gc, ngvec, cfunig, ivg, ivgig
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
                     & cpumtaa, cpumtalo, cpumtloa, cpumtlolo
      use summations, only: doublesummation_simple_cz
      use m_getapwcmt
      use m_getlocmt
      use m_putemat
      use m_emattim
      use m_getunit
      use m_genfilname

use m_writecmplxparts
use modxs, only: ngqmax
use mod_Gkvector, only: ngkmax
use mod_APW_LO, only: lolmax

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
      integer :: whichthread, ig, ix, igs, ist2, ist1
      real(8) :: cpuini, cpuread, cpumain, cpuwrite, cpuall
      real(8) :: cpugnt, cpumt, cpuir, cpufft
      real(8) :: cpumalores, cpumloares
      real(8) :: cpumlolores, cpumirres, cpudbg
      real(8) :: cpu0, cpu1, cpu00, cpu01, cpugntlocal, cpumtlocal
      real(8) :: vkql(3)
      integer :: shift(3), iv(3)
      type(fftmap_type) :: fftmap
      real(8) :: emat_gmax

      ! If task 330 is 'writeemat'
      if(task .eq. 330) then
        call chkpt(3, (/ task, iq, ik /),&
          & 'b_ematqk: task, q - point index, k - point index; q - dependent matrix elements')
      end if

      call timesec(cpu0)

      ! Find k+q-point
      ikq = ikmapikq(ik, iq)

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

      ! Get number of G+k and G+k+q vectors
      n0 = ngk0(1, ik)
      n = ngk(1, ikq)

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

      ! Read eigenvectors, eigenvalues and occupancies for k+q
      !   Read first variational eigenvectors from EVECFV_QMTXXX.OUT 
      !   (file extension needs to be set by calling routine)
      call getevecfv(vkl(1, ikq), vgkl(1, 1, 1, ikq), evecfv)

      ! Save local orbital coefficients
      evecfvu(:, :) = evecfv(ngk(1, ikq)+1:ngk(1, ikq)+nlotot, bc%il2:bc%iu2, 1)

      ! Read eigenvectors, eigenvalues and occupancies for k (q=0)
      !   Read first variational eigenvectors from EVECFV_QMT000.OUT 
      call getevecfv0(vkl0(1, ik), vgkl0(1, 1, 1, ik), evecfv0)

      ! Save local orbital coefficients
      evecfvo0(:, :) = evecfv0(ngk0(1, ik)+1:ngk0(1, ik)+nlotot, bc%il1:bc%iu1, 1)

!write(*,*) "Sagmeister grid vars"
!write(*,*) "  gkmax = ", gkmax
!write(*,*) "  gqmax = ", input%xs%gqmax
!write(*,*) "  gmaxvr = ", input%groundstate%gmaxvr
!write(*,*) "  ngkmax = ", ngkmax
!write(*,*) "  ngqmax = ", ngqmax
!write(*,*) "  (used) input%xs%lmaxapwwf = ", input%xs%lmaxapwwf
!write(*,*) "  input%xs%lmaxapw = ", input%xs%lmaxapw
!write(*,*) "  input%xs%lmaxmat = ", input%xs%lmaxmat
!write(*,*) "  input%groundstate%lmaxmat = ", input%groundstate%lmaxmat
!write(*,*) "  (used if larger that lmaxapwwf) lolmax = ", lolmax
!write(*,*) "  (used) input%xs%lmaxemat = ", input%xs%lmaxemat
!write(*,*) "Writing eigenvectors sagmeister"
!write(*,*) "ik=", ik
!write(*,*) "vkl=", vkl0(1:3, ik)
!write(*,*) "iq=", iq
!write(*,*) "vql=", vql(1:3, iq)
!write(*,*) "ikq=", ikq
!write(*,*) "vkql=", vkl(1:3, ikq)
!write(*,*) "ngk=", ngk0(1,ik) 
!write(*,*) "ngkq=", ngk(1,ikq)
!write(*,*) "ngq=", ngq(iq)

      ! Determine number of radial functions used in APW 
      ! basis functions per species
      apwmaxsize=0
      allocate(apwsize(nspecies))
      do is=1, nspecies
        naug=0
        do l=0, input%xs%lmaxapwwf
          naug=naug+(2*l+1)*apword(l, is)
        end do
        apwsize(is)=naug
        apwmaxsize=max(apwmaxsize, naug)
      end do

      ! Determine number of radial functions used in 
      ! LO basis functions per species
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

      ! Get APW expansion coefficients for k
      call getapwcmt(0, ik, 1, nstfv, input%xs%lmaxapwwf, apwcmt0)
      ilo=0
      do is=1, nspecies
        do ia=1, natoms(is)
          naug=0
          ias = idxas(ia, is)
          do l=0, input%xs%lmaxapwwf
            do m=-l, l
              do io=1, apword(l, is)
                naug=naug+1
                lm=idxlm(l, m)
                cmtfun0(1:bc%n1, naug, ias)=conjg(apwcmt0(bc%il1:bc%iu1, io, lm, ias))
              end do
            end do
          end do
          cmtfun0(1:bc%n1, naug+1:naug+losize(is), ias)=&
            & conjg(transpose(evecfvo0(ilo+1:ilo+losize(is), 1:bc%n1)))
          ilo=ilo+losize(is)
        end do
      end do

      ! Get APW expansion coefficients for k+q
      call getapwcmt(iq, ikq, 1, nstfv, input%xs%lmaxapwwf, apwcmt)
      ilo=0
      do is=1, nspecies
        do ia=1, natoms(is)
          naug=0
          ias = idxas(ia, is)
          do l=0, input%xs%lmaxapwwf
            do m=-l, l
              do io=1, apword(l, is)
                naug=naug+1
                lm=idxlm(l, m)
                cmtfun(1:bc%n2, naug, ias)=apwcmt(bc%il2:bc%iu2, io, lm, ias)
              end do
            end do
          end do
          cmtfun(1:bc%n2, naug+1:naug+losize(is), ias)=&
            & transpose(evecfvu(ilo+1:ilo+losize(is), 1:bc%n2))
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
        call b_ematgntsum(iq, igq, integrals)
        call timesec(cpu01)
        if(whichthread.eq.0) cpugnt = cpugnt + cpu01 - cpu00
        ! Muffin-tin contribution
        call b_ematqkgmt(iq, ik, igq, integrals, emat(:,:,igq), bc)
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
      evecfvu2(:, :) = evecfv(1:ngk(1, ikq), bc%il2:bc%iu2, 1)
      evecfvo20(:, :) = evecfv0(1:ngk0(1, ik), bc%il1:bc%iu1, 1)
     
      ! Interstitial contribution
      mm: if(input%xs%pwmat.eq.'mm') then
        ! Matrix-matrix multiplication
#ifdef USEOMP
    !$omp parallel default(shared) private(igq, xihir, cpu00, cpu01, whichthread)
        whichthread=omp_get_thread_num()
#endif 
        allocate(xihir(n0, n))
#ifdef USEOMP
    !$omp do
#endif
        do igq = 1, ngq(iq)
          call timesec(cpu00)
          ! Interstitial contribution
          call b_ematqkgir(iq, ik, igq, xihir, n0, n)
          call timesec(cpu01)
          if(whichthread.eq.0) cpuir = cpuir + cpu01 - cpu00

          ! Interstitial contribution
          call doublesummation_simple_cz(emat(:, :, igq), evecfvo20,&
            & xihir, evecfvu2, zone, zone, .true.)
          call timesec(cpu00)
          if(whichthread.eq.0) cpumirres = cpumirres + cpu00 - cpu01

          call timesec(cpu01)
          if(whichthread.eq.0) cpudbg = cpudbg + cpu01 - cpu00
        end do ! Igq
#ifdef USEOMP
    !$omp end do
#endif
        deallocate(xihir)
#ifdef USEOMP
    !$omp end parallel
#endif

      else mm
        ! Fourier transforms
        call timesec(cpu00)
     
        ikq = ikmapikq(ik, iq)
        vkql(:)=vkl(:, ik)+vql(:, iq)
          
        ! Umklapp treatment
        do ix=1, 3
          if(vkql(ix).ge.1d0-1d-13) then
            shift(ix)=-1
          else
            shift(ix)=0
          end if
        end do
          
        emat_gmax=2*gkmax +input%xs%gqmax

        call genfftmap(fftmap, emat_gmax)
        allocate(zfft0(fftmap%ngrtot+1, bc%n1))
        zfft0=zzero

        allocate(zfftcf(fftmap%ngrtot+1))
        zfftcf=zzero

        do ig=1, ngvec
          if(gc(ig).lt.emat_gmax) then
            zfftcf(fftmap%igfft(ig))=cfunig(ig)
          end if
        end do
        call zfftifc(3, fftmap%ngrid, 1, zfftcf)

#ifdef USEOMP
    !$omp parallel default(shared) private(ist1, ig)
    !$omp do
#endif
        do ist1=1, bc%n1
          do ig=1, ngk0(1, ik)
            zfft0(fftmap%igfft(igkig0(ig, 1, ik)), ist1)=evecfvo20(ig, ist1)
          end do
          call zfftifc(3, fftmap%ngrid, 1, zfft0(:, ist1))
          zfft0(:, ist1)=conjg(zfft0(:, ist1))*zfftcf(:) 
        end do
#ifdef USEOMP
    !$omp end do
    !$omp end parallel
#endif
       
#ifdef USEOMP
    !$omp parallel default(shared) private(ist1, ist2, ig, zfft, iv, igs, zfftres, igq)
#endif
        allocate(zfftres(fftmap%ngrtot+1))
        allocate(zfft(fftmap%ngrtot+1))
#ifdef USEOMP
    !$omp do
#endif
        do ist2=1, bc%n2
          zfft=zzero
          if(sum(shift).ne.0) then
            do ig=1, ngk(1, ikq)
              iv=ivg(:, igkig(ig, 1, ikq))+shift
              igs=ivgig(iv(1), iv(2), iv(3))
              zfft(fftmap%igfft(igs))=evecfvu2(ig, ist2)
            end do
          else
            do ig=1, ngk(1, ikq)
              zfft(fftmap%igfft(igkig(ig, 1, ikq)))=evecfvu2(ig, ist2)
            end do
          end if
          call zfftifc(3, fftmap%ngrid, 1, zfft)

          do ist1=1, bc%n1
            do ig=1, fftmap%ngrtot
              zfftres(ig)=zfft(ig)*zfft0(ig, ist1) 
            end do

            call zfftifc(3, fftmap%ngrid, -1, zfftres)
            do igq=1, ngq(iq)
              emat(ist1, ist2, igq)=emat(ist1, ist2, igq)&
                &+ zfftres(fftmap%igfft(igqig(igq, iq)))
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

      end if mm

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

    end subroutine b_ematqk
    !EOC

    subroutine b_ematqkgmt(iq, ik, igq, integrals, emat, bc)
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
    end subroutine b_ematqkgmt

    subroutine b_ematgntsum(iq, igq, integrals)
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
          do cl1 = 1, l1shape
            l1 = l1map(cl1)
            do cm1 = 1, m1shape(l1)
              m1 = m1map(l1, cm1)
              lm1 = idxlm(l1, m1)
              do io1 = 1, apword(l1, is)
                do cl2 = 1, l2shape(l1, m1)
                  l3 = l2map(l1, m1, cl2)
                  do cm2 = 1, m2shape(l1, m1, l3)
                    m3 = m2map(l1, m1, l3, cm2)
                    lm3 = idxlm(l3, m3)
                    do io2 = 1, apword(l3, is)
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

          iaug2=0
          do l3=0, input%xs%lmaxapwwf
            do m3=-l3, l3
              do io2=1, apword(l3, is)
                iaug2=iaug2+1
                lm3=idxlm(l3, m3)
                iaug1=0
                do l1=0, input%xs%lmaxapwwf
                  do m1=-l1, l1
                      do io1=1, apword(l1, is)
                        iaug1=iaug1+1
                        lm1=idxlm(l1, m1)
                        integrals(iaug2, iaug1, ias)=intrgaa(lm1, io1, lm3, io2)
                      end do
                    end do
                  end do
              end do
            end do
          end do

          !-------------------------------------!
          !     apw-local-orbital integrals     !
          !-------------------------------------!
          do cl1 = 1, l1shape
            l1 = l1map(cl1)
            do cm1 = 1, m1shape(l1)
              m1 = m1map(l1, cm1)
              lm1 = idxlm(l1, m1)
              do io = 1, apword(l1, is)
                do ilo = 1, nlorb(is)
                  l3 = lorbl(ilo, is)
                  do cm2 = 1, m2shape(l1, m1, l3)
                    m3 = m2map(l1, m1, l3, cm2)
                    lm3 = idxlm(l3, m3)
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

          iaug2=0
          do ilo = 1, nlorb(is)
            l3 = lorbl(ilo, is)
            do m3=-l3, l3
              iaug2=iaug2+1
              lm3=idxlm(l3, m3)
             
              iaug1=0
              do l1=0, input%xs%lmaxapwwf
                do m1=-l1, l1
                  do io=1, apword(l1, is)
                    iaug1=iaug1+1
                    lm1=idxlm(l1, m1)
                    integrals(apwsize(is)+iaug2, iaug1, ias)=intrgalo(m3, ilo, lm1, io)
                  end do
                end do
               end do

            end do
          end do

          !-------------------------------------!
          !     local-orbital-apw integrals     !
          !-------------------------------------!
          do ilo = 1, nlorb(is)
            l1 = lorbl(ilo, is)
            do cm1 = 1, m1shape(l1)
              m1 = m1map(l1, cm1)
              lm1 = idxlm(l1, m1)
              do cl2 = 1, l2shape(l1, m1)
                l3 = l2map(l1, m1, cl2)
                do cm2 = 1, m2shape(l1, m1, l3)
                  m3 = m2map(l1, m1, l3, cm2)
                  lm3 = idxlm(l3, m3)
                  do io = 1, apword(l3, is)
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

          iaug2=0
          do ilo = 1, nlorb(is)
            l3 = lorbl(ilo, is)
            do m3=-l3, l3
              iaug2=iaug2+1
              lm3=idxlm(l3, m3)

              iaug1=0
              do l1=0, input%xs%lmaxapwwf
                do m1=-l1, l1
                  do io=1, apword(l1, is)
                    iaug1=iaug1+1
                    lm1=idxlm(l1, m1)
                    integrals(iaug1, apwsize(is)+iaug2, ias)=intrgloa(m3, ilo, lm1, io)
                  end do
                end do
               end do

            end do
          end do

          !-----------------------------------------------!
          !     local-orbital-local-orbital integrals     !
          !-----------------------------------------------!
          do ilo1 = 1, nlorb(is)
            l1 = lorbl(ilo1, is)
            do cm1 = 1, m1shape(l1)
              m1 = m1map(l1, cm1)
              lm1 = idxlm(l1, m1)
              do ilo2 = 1, nlorb(is)
                l3 = lorbl(ilo2, is)
                do cm2 = 1, m2shape(l1, m1, l3)
                  m3 = m2map(l1, m1, l3, cm2)
                  lm3 = idxlm(l3, m3)
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

          iaug2=0
          do ilo2 = 1, nlorb(is)
            l3 = lorbl(ilo2, is)
            do m3=-l3, l3
              iaug2=iaug2+1

              iaug1=0
              do ilo1 = 1, nlorb(is)
                l1 = lorbl(ilo1, is)
                do m1=-l1, l1
                  iaug1=iaug1+1
                  integrals(apwsize(is)+iaug2, apwsize(is)+iaug1, ias)=&
                    & intrglolo(m1, ilo1, m3, ilo2)
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
          
    end subroutine b_ematgntsum

    subroutine b_ematqkgir(iq, ik, igq, xihir, n0, n)
      use mod_Gvector, only: ivg, ivgig, cfunig
      use mod_Gkvector, only: ngkmax, igkig, ngk
      use mod_qpoint, only: vql
      use mod_kpoint, only: vkl
      use modxs, only: ikmapikq, ngkmax0, vkl0, igkig0,&
                     & igqig, ngk0

      implicit none

      ! Arguments
      integer, intent(in) :: iq, ik, igq, n0,n
      complex(8), intent(out) :: xihir(n0,n) 

      ! Local variables
      character(*), parameter :: thisnam = 'b_ematqkgir'
      integer :: ikq, ig, ig1, ig2, ig3, igk0, igk, iv(3), iv1(3), iv3(3), ivu(3)
      integer, allocatable :: aigk0(:), aigk(:)

      ! Grid index for k+q point
      ikq = ikmapikq(ik, iq)
      allocate(aigk0(ngkmax0), aigk(ngkmax))

      ! Positive umklapp G-vector
      ivu(:) = nint(vkl0(:, ik)+vql(:, iq)-vkl(:, ikq))

      ! Precalculate for speedup
      aigk0(:) = igkig0(:, 1, ik)
      aigk(:) = igkig(:, 1, ikq)
      ig3 = igqig(igq, iq)
      iv3(:) = ivg(:, ig3)

      do igk0 = 1, ngk0(1, ik)
        ig1 = aigk0(igk0)
        iv1(:) = ivg(:, ig1) + iv3(:)
        do igk = 1, ngk(1, ikq)
          ig2 = aigk(igk)
          ! Umklapp of k+q vector included
          iv(:) = iv1(:) - (ivg(:, ig2)-ivu(:))
          ig = ivgig(iv(1), iv(2), iv(3))
          xihir(igk0, igk) = cfunig(ig)
        end do
      end do

      deallocate(aigk0, aigk)
    end subroutine b_ematqkgir

end module m_b_ematqk
