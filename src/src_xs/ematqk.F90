! Copyright(C) 2006-2008 S. Sagmeister and Claudia Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
subroutine ematqk(iq, ik)
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
  use modxs, only: ikmapikq, msg, ngk0, xiohalo,&
                 & xiuhloa, nst1, nst2, evecfv,&
                 & istl1, istl2, istu1, istu2,&
                 & vkl0, vgkl0, evecfv0, evalsv0,&
                 & occsv0, apwmaxsize, apwsize, losize,&
                 & lomaxsize, cmtfun0, cmtfun, apwcmt0,&
                 & apwcmt, xiou, ngq, igqig,&
                 & fnetim, fftmap_type, igkig0,&
                 & cpumtaa, cpumtalo, cpumtloa, cpumtlolo
  ! One routine modules:
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

  implicit none
      
  ! Arguments
  integer, intent(inout) :: iq, ik

  ! Local variables
  character(*), parameter :: thisnam = 'ematqk'
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
 
  ! Task 330 is 'writeemat'
  if(task .eq. 330) then
    call chkpt(3, (/ task, iq, ik /),&
      & 'ematqk: task, q - point index, k - point index; q - dependent matrix elements')
  end if

  call timesec(cpu0)

  ! Find k+q-point
  ikq = ikmapikq(ik, iq)

  ! Check for stop statement
  write(msg, *) 'for q-point', iq, ': k-point:', ik - 1, ' finished'
  call xschkstop

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

  n0 = ngk0(1, ik)
  n = ngk(1, ikq)

  ! Allocate matrix elements array
  if(allocated(xiohalo)) deallocate(xiohalo)
  allocate(xiohalo(nst1, nlotot))
  if(allocated(xiuhloa)) deallocate(xiuhloa)
  allocate(xiuhloa(nlotot, nst2))
  ! Allocate temporary arrays
  allocate(evecfvo0(nlotot, nst1))
  allocate(evecfvu(nlotot, nst2))
  ! Zero arrays
  xiohalo(:, :) = zzero
  xiuhloa(:, :) = zzero

  call timesec(cpu1)
  cpuini = cpu1 - cpu0

  ! Read eigenvectors, eigenvalues and occupancies for G+k+q
  call getevecfv(vkl(1, ikq), vgkl(1, 1, 1, ikq), evecfv)
  call getevalsv(vkl(1, ikq), evalsv(1, ikq))
  
  ! Read occupation numbers for G+k+q
  call getoccsv(vkl(1, ikq), occsv(1, ikq))
  evecfvu(:, :) = evecfv(ngk(1, ikq)+1:ngk(1, ikq)+nlotot, istl2:istu2, 1)

  ! Read eigenvectors, eigenvalues and occupancies for G+k(q=0)
  call getevecfv0(vkl0(1, ik), vgkl0(1, 1, 1, ik), evecfv0)
  call getevalsv0(vkl0(1, ik), evalsv0(1, ik))
  
  ! Read occupation numbers for G+k
  call getoccsv0(vkl0(1, ik), occsv0(1, ik))
  evecfvo0(:, :) = evecfv0(ngk0(1, ik)+1:ngk0(1, ik)+nlotot, istl1:istu1, 1)

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

  allocate(losize(nspecies))
  lomaxsize=0
  losize=0
  do is=1, nspecies
    do ilo = 1, nlorb(is)
      losize(is)=losize(is)+2*lorbl(ilo, is)+1
    end do
    lomaxsize=max(lomaxsize, losize(is))
  end do

  allocate(cmtfun0(nst1, apwmaxsize+lomaxsize, natmtot))
  cmtfun0=zzero
  allocate(cmtfun(nst2, apwmaxsize+lomaxsize, natmtot))
  cmtfun=zzero
  apwcmt0=zzero
  apwcmt=zzero

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
            cmtfun0(1:nst1, naug, ias)=conjg(apwcmt0(istl1:istu1, io, lm, ias))
          end do
        end do
      end do
      cmtfun0(1:nst1, naug+1:naug+losize(is), ias)=&
        & conjg(transpose(evecfvo0(ilo+1:ilo+losize(is), 1:nst1)))
      ilo=ilo+losize(is)
    end do
  end do

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
            cmtfun(1:nst2, naug, ias)=apwcmt(istl2:istu2, io, lm, ias)
          end do
        end do
      end do
      cmtfun(1:nst2, naug+1:naug+losize(is), ias)=&
        & transpose(evecfvu(ilo+1:ilo+losize(is), 1:nst2))
      ilo=ilo+losize(is)
    end do
  end do

  call timesec(cpu0)
  cpuread = cpu0 - cpu1

  ! Zero matrix elements array
  xiou(:, :, :) = zzero

  whichthread=0

  ! Loop over G+q vectors
  cpugntlocal=0.0d0
  cpumtlocal=0.0d0

#ifdef USEOMP
!$omp parallel default(shared) private(igq, integrals, cpu00, cpu01, whichthread)
!print *, whichthread
#endif
  allocate(integrals(apwmaxsize+lomaxsize, apwmaxsize+lomaxsize, natmtot))
#ifdef USEOMP
  whichthread=omp_get_thread_num()
!$omp do
#endif
  do igq = 1, ngq(iq)
    call timesec(cpu00)
    ! Summation of gaunt coefficients wrt radial integrals
    call ematgntsum(iq, igq, integrals)
    call timesec(cpu01)
    if(whichthread.eq.0) cpugnt = cpugnt + cpu01 - cpu00
    ! Muffin-tin contribution
    call ematqkgmt(iq, ik, igq, integrals)
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
 
  allocate(evecfvo20(n0, nst1))
  allocate(evecfvu2(n, nst2))
  evecfvu2(:, :) = evecfv(1:ngk(1, ikq), istl2:istu2, 1)
  evecfvo20(:, :) = evecfv0(1:ngk0(1, ik), istl1:istu1, 1)
 
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
      call ematqkgir(iq, ik, igq, xihir, n0, n)
      call timesec(cpu01)
      if(whichthread.eq.0) cpuir = cpuir + cpu01 - cpu00

      ! Interstitial contribution
      call doublesummation_simple_cz(xiou(:, :, igq), evecfvo20,&
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
    allocate(zfft0(fftmap%ngrtot+1, nst1))
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
    do ist1=1, nst1
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
    do ist2=1, nst2
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

      do ist1=1, nst1
        do ig=1, fftmap%ngrtot
          zfftres(ig)=zfft(ig)*zfft0(ig, ist1) 
        end do

        call zfftifc(3, fftmap%ngrid, -1, zfftres)
        do igq=1, ngq(iq)
          xiou(ist1, ist2, igq)=xiou(ist1, ist2, igq)&
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

end subroutine ematqk
