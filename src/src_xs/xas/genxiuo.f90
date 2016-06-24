!
!
!
! Copyright (C) 2006-2008 S. Sagmeister and Claudia Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine genxiuo (iq, ik, xi)
      Use mod_constants
      Use mod_eigenvalue_occupancy
      Use mod_misc
      Use mod_Gkvector
      Use mod_APW_LO
      Use mod_Gvector
      Use mod_kpoint
      Use modinput
      Use modmpi
      Use modxs
      Use summations
      Use m_getapwcmt
      Use m_getlocmt
      Use m_putemat
      Use m_emattim
      Use m_getunit
      Use m_genfilname
      use modmain
      Use modxas
#ifdef USEOMP
      use omp_lib
#endif

      Implicit None
!      type(mtints_type) :: integrals
      
  ! arguments
      Integer, Intent (InOut) :: iq, ik
      Complex(8), Intent(Out) :: xi(xasstop-xasstart+1, sto2-sta2+1, ngq(iq))
  ! local variables
      Character (*), Parameter :: thisnam = 'ematqk'
  ! allocatable arrays
      Complex (8), Allocatable :: evecfvo0 (:, :)
      Complex (8), Allocatable :: evecfvu (:, :)
      Complex (8), Allocatable :: evecfvo20 (:, :)
      Complex (8), Allocatable :: evecfvu2 (:, :)
      Complex (8), Allocatable :: apwalmt (:, :, :, :), apwalmt0 (:, :, :, :)
      Complex(8), allocatable :: wfmt(:,:,:), wfmt0(:,:,:)
      Complex (8), Allocatable :: zfft0(:,:),zfft(:),zfftres(:),zfftcf(:)
      Complex (8), Allocatable :: xihir (:, :)
  ! expansion coefficients of APW and LO functions
      Complex (8), Allocatable :: integral(:,:,:,:)

      Integer :: ikq, igq, n, n0, is,l,m,io,naug,ia,ias,lm,ilo,whichthread,ig,ix,igs,ist2,ist1
      Real (8) :: cpuini, cpuread, cpumain, cpuwrite, cpuall
      Real (8) :: cpugnt, cpumt, cpuir, cpufft
      Real (8) :: cpumalores, cpumloares
      Real (8) :: cpumlolores, cpumirres, cpudbg
      Real (8) :: cpu0, cpu1, cpu00, cpu01, cpugntlocal, cpumtlocal
      Real (8) :: vkql(3)
      integer :: shift(3),iv(3)
      logical :: umklapp
      type (fftmap_type) :: fftmap
      Real (8) :: emat_gmax,ta,tb


!      write(*,*) 'nst',nst1,nst2
!      write(*,*) 'istl/istu',istl1,istu1,istl2,istu2

!      write(*,*) 'howdy'
!      read(*,*) emat_gmax
!      call genfftmap(fftmap,emat_gmax)

!      write(*,*) size(ngk)
!      write(*,*) size(ngk0)
!      stop
!
!      ik=64
!      iq=8

!     write(*,*) ik,iq
 
      If (task .Eq. 330) Call chkpt (3, (/ task, iq, ik /), 'ematqk: ta&
     &sk, q - point index, k - point index; q - dependent matrix elemen&
     &ts')
      Call timesec (cpu0)
  ! find k+q-point
      ikq = ikmapikq (ik, iq)
  ! check for stop statement
      Write (msg,*) 'for q-point', iq, ': k-point:', ik - 1, ' finished&
     &'
      Call xschkstop
!
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
!
  ! allocate temporary arrays
      n0 = ngk0 (1, ik)
      n = ngk (1, ikq)
  ! allocate matrix elements array
      If (allocated(xiohalo)) deallocate (xiohalo)
      Allocate (xiohalo(nst1, nlotot))
      If (allocated(xiuhloa)) deallocate (xiuhloa)
      Allocate (xiuhloa(nlotot, nst2))
  ! allocate temporary arrays
      Allocate (evecfvo0(nlotot, nst1))
      Allocate (evecfvu(nlotot, nst2))
      Allocate (apwalmt(ngkmax, apwordmax, lmmaxapw, natmtot))
      Allocate (apwalmt0(ngkmax, apwordmax, lmmaxapw, natmtot))
      Allocate(wfmt(lmmaxapw,nrcmtmax,nstfv))
      Allocate(wfmt0(lmmaxapw,nrcmtmax,nstfv))
!      Allocate (evecfvo20(n0, nst1))
!      Allocate (evecfvu2(n, nst2))
!      Allocate (xihir(n0, n))
  ! zero arrays
      xiohalo (:, :) = zzero
      xiuhloa (:, :) = zzero
!
      Call timesec (cpu1)
      cpuini = cpu1 - cpu0
!
  ! read eigenvectors, eigenvalues and occupancies for G+k+q
      Call getevecfv (vkl(1, ikq), vgkl(1, 1, 1, ikq), evecfv)
      Call getevalsv (vkl(1, ikq), evalsv(1, ikq))
  ! read occupation numbers for G+k+q
      Call getoccsv (vkl(1, ikq), occsv(1, ikq))
      evecfvu (:, :) = evecfv (ngk(1, ikq)+1:ngk(1, ikq)+nlotot, istl2:istu2, 1)
  ! read matching coefficients and  radial wavefunction for G+k+q
	  Call match (ngk(1, ikq), gkc(1, 1, ikq), tpgkc(1, 1, 1, ikq), &
        & sfacgk(1, 1, 1, ikq), apwalmt)
      Call wavefmt(input%groundstate%lradstep, &
        &  input%groundstate%lmaxapw,input%xs%bse%xasspecies,input%xs%bse%xasatom,ngk(1, ikq),apwalmt, &
        &  evecfv,lmmaxapw,wfmt)

  ! read eigenvectors, eigenvalues and occupancies for G+k (q=0)
      Call getevecfv0 (vkl0(1, ik), vgkl0(1, 1, 1, ik), evecfv0)
      Call getevalsv0 (vkl0(1, ik), evalsv0(1, ik))
  ! read occupation numbers for G+k
      Call getoccsv0 (vkl0(1, ik), occsv0(1, ik))
      evecfvo0 (:, :) = evecfv0 (ngk0(1, ik)+1:ngk0(1, ik)+nlotot, &
     & istl1:istu1, 1)
  ! read matching coefficients and  radial wavefunction for G+k (q=0)
	  Call match (ngk(1, ik), gkc(1, 1, ik), tpgkc(1, 1, 1, ik), &
        & sfacgk(1, 1, 1, ik), apwalmt0)
      Call wavefmt(input%groundstate%lradstep, &
        &  input%groundstate%lmaxapw,input%xs%bse%xasspecies,input%xs%bse%xasatom,ngk(1, ik),apwalmt0, &
        &  evecfv0,lmmaxapw,wfmt0)


      apwmaxsize=0
      allocate(apwsize(nspecies))
      do is=1,nspecies
        naug=0
        do l=0, input%xs%lmaxapwwf
          naug=naug+(2*l+1)*apword(l,is)
        enddo
        apwsize(is)=naug
        apwmaxsize=max(apwmaxsize,naug)
      enddo

      allocate(losize(nspecies))
      lomaxsize=0
      losize=0
      do is=1,nspecies
        Do ilo = 1, nlorb (is)
          losize(is)=losize(is)+2*lorbl (ilo, is)+1
        enddo
        lomaxsize=max(lomaxsize,losize(is))
      enddo



      allocate(cmtfun0(nst1,apwmaxsize+lomaxsize,natmtot))
      cmtfun0=zzero
      allocate(cmtfun(nst2,apwmaxsize+lomaxsize,natmtot))
      cmtfun=zzero
      apwcmt0=zzero
      apwcmt=zzero


      Call timesec (cpu0)
      cpuread = cpu0 - cpu1
!
  ! zero matrix elements array
      xi (:, :, :) = zzero
if (.true.) then
!
      whichthread=0
  ! loop over G+q vectors
      cpugntlocal=0d0
      cpumtlocal=0d0
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(igq,integral,cpu00,cpu01,whichthread)
!print *,whichthread
#endif
      
#ifdef USEOMP
      whichthread=omp_get_thread_num()
!$OMP DO
#endif
      Do igq = 1, ngq (iq)
         Call timesec (cpu00)
     ! summation of Gaunt coefficients wrt radial integrals
		 Allocate (integral(input%xs%lmaxemat+1,lmmaxapw,xasstop-xasstart+1,sto2-sta2+1))
         Call ematradou (iq, ik, igq,wfmt0,integral)
         Call timesec (cpu01)
         if (whichthread.eq.0) cpugnt = cpugnt + cpu01 - cpu00
     ! muffin-tin contribution
         Call ematsumuo (iq, ik, igq,integral, xi)
         Call timesec (cpu00)
         if (whichthread.eq.0) cpumt = cpumt + cpu00 - cpu01
		 deallocate(integral) 
!
      End Do ! igq
#ifdef USEOMP
!$OMP END DO
#endif
           
!      deallocate(integrals%aa,integrals%alo,integrals%loa,integrals%lolo)
     
#ifdef USEOMP
!$OMP END PARALLEL
#endif
      deallocate (apwsize,losize,cmtfun,cmtfun0)
      Deallocate (evecfvu, evecfvo0)




      Allocate (evecfvo20(n0, nst1))
      Allocate (evecfvu2(n, nst2))
      evecfvu2 (:, :) = evecfv (1:ngk(1, ikq), istl2:istu2, 1)
      evecfvo20 (:, :) = evecfv0 (1:ngk0(1, ik), istl1:istu1, 1)

! interstitial contribution
if (input%xs%pwmat.eq.'mm') then
! matrix-matrix multiplication

#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(igq,xihir,cpu00,cpu01,whichthread)
      whichthread=omp_get_thread_num()
#endif 
      Allocate(xihir(n0,n))
!#ifdef USEOMP
!!$OMP DO
!#endif
!      Do igq = 1, ngq (iq)
!!         Call terminateqry ('ematqk')
!         Call timesec (cpu00)
!     ! interstitial contribution
!         Call ematqkgir (iq, ik, igq,xihir,n0,n)
!         Call timesec (cpu01)
!         if (whichthread.eq.0) cpuir = cpuir + cpu01 - cpu00

!     ! interstitial contribution
!         Call doublesummation_simple_cz (xiou(:, :, igq), evecfvo20, xihir, evecfvu2, zone, zone, .True.)
!         Call timesec (cpu00)
!         if (whichthread.eq.0) cpumirres = cpumirres + cpu00 - cpu01

!         Call timesec (cpu01)
!         if (whichthread.eq.0) cpudbg = cpudbg + cpu01 - cpu00
!      End Do ! igq
!#ifdef USEOMP
!!$OMP END DO
!#endif
!      deallocate(xihir)
#ifdef USEOMP
!$OMP END PARALLEL
#endif

! elseif (input%xs%pwmat.eq.'pwmat') then
else
! fourier transforms

      call timesec(cpu00)
 
      ikq = ikmapikq (ik, iq)
      vkql(:)=vkl(:,ik)+vql(:,iq)
      
! umklapp treatment
      do ix=1,3
        if (vkql(ix).ge.1d0-1d-13) then
          shift(ix)=-1
        else
          shift(ix)=0
        endif
      enddo
      
      emat_gmax=2*gkmax +input%xs%gqmax
        

      call genfftmap(fftmap,emat_gmax)
      allocate(zfft0(fftmap%ngrtot+1,nst1))
      zfft0=zzero

        allocate(zfftcf(fftmap%ngrtot+1))
        zfftcf=zzero
        do ig=1,ngvec
         if (gc(ig).lt.emat_gmax) then
          zfftcf(fftmap%igfft(ig))=cfunig(ig)
         endif
        enddo
        Call zfftifc (3, fftmap%ngrid,1, zfftcf)

#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ist1,ig)
!$OMP DO
#endif
      do ist1=1,nst1
        do ig=1,ngk0(1, ik)
          zfft0(fftmap%igfft(igkig0(ig,1,ik)),ist1)=evecfvo20(ig,ist1)
        enddo
        Call zfftifc (3, fftmap%ngrid,1, zfft0(:,ist1))
        zfft0(:,ist1)=conjg(zfft0(:,ist1))*zfftcf(:) 
      enddo
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
   
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ist1,ist2,ig,zfft,iv,igs,zfftres,igq)
#endif
       allocate(zfftres(fftmap%ngrtot+1))
       allocate(zfft(fftmap%ngrtot+1))
#ifdef USEOMP
!$OMP DO
#endif
      do ist2=1,nst2
        zfft=zzero
        if (sum(shift).ne.0) then
          do ig=1,ngk(1, ikq)
            iv=ivg (:, igkig(ig,1,ikq))+shift
            igs=ivgig(iv(1),iv(2),iv(3))
            zfft(fftmap%igfft(igs))=evecfvu2(ig,ist2)
          enddo
        else
          do ig=1,ngk(1, ikq)
            zfft(fftmap%igfft(igkig(ig,1,ikq)))=evecfvu2(ig,ist2)
          enddo
        endif
        Call zfftifc (3, fftmap%ngrid,1, zfft)

        do ist1=1,nst1
          do ig=1,fftmap%ngrtot
            zfftres(ig)=zfft(ig)*zfft0(ig,ist1) 
          enddo

          Call zfftifc (3, fftmap%ngrid,-1, zfftres)
          do igq=1,ngq(iq)
            xiou(ist1,ist2,igq)=xiou(ist1,ist2,igq)+zfftres(fftmap%igfft(igqig (igq, iq)))
          enddo
        enddo
      enddo
#ifdef USEOMP
!$OMP END DO
#endif
      deallocate(zfftres,zfft)
#ifdef USEOMP
!$OMP END PARALLEL
#endif

      deallocate(fftmap%igfft)
      deallocate(zfft0,zfftcf)
      call timesec(cpu01)
      cpufft=cpu01-cpu00

endif
!
endif
      Call timesec (cpu1)
      cpumain = cpu1 - cpu0
!
!      Deallocate (evecfvu2, evecfvo20)

      Call timesec (cpu0)
      cpuwrite = cpu0 - cpu1
      cpuall = cpuini + cpuread + cpumain + cpuwrite
!
  ! write timing information
      If ((task.ne.440).and.(task .Ne. 441) &
     & .And. (task .Ne. 450) .And. (task .Ne. 451)) Then
         Call emattim (iq, ik, trim(fnetim), cpuini, cpuread, cpumain, &
        & cpuwrite, cpuall, cpugnt, cpumt, cpuir, cpumalores, &
        & cpumloares, cpumlolores, cpumirres, cpudbg, cpumtaa, &
        & cpumtalo, cpumtloa, cpumtlolo,cpufft)
      End If

!
End Subroutine genxiuo
