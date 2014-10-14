!
!
!
! Copyright (C) 2006-2008 S. Sagmeister and Claudia Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine ematqk (iq, ik)
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
#ifdef USEOMP
      use omp_lib
#endif

      Implicit None
      type(mtints_type) :: integrals
      
  ! arguments
      Integer, Intent (In) :: iq, ik
  ! local variables
      Character (*), Parameter :: thisnam = 'ematqk'
  ! allocatable arrays
      Complex (8), Allocatable :: evecfvo0 (:, :)
      Complex (8), Allocatable :: evecfvu (:, :)
      Complex (8), Allocatable :: evecfvo20 (:, :)
      Complex (8), Allocatable :: evecfvu2 (:, :)
      Integer :: ikq, igq, n, n0, is,l,m,io,naug,ia,ias,lm,ilo,whichthread
      Real (8) :: cpuini, cpuread, cpumain, cpuwrite, cpuall
      Real (8) :: cpugnt, cpumt, cpuir
      Real (8) :: cpumalores, cpumloares
      Real (8) :: cpumlolores, cpumirres, cpudbg
      Real (8) :: cpu0, cpu1, cpu00, cpu01, cpugntlocal, cpumtlocal
!
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
!      evecfvu2 (:, :) = evecfv (1:ngk(1, ikq), istl2:istu2, 1)
!
  ! read eigenvectors, eigenvalues and occupancies for G+k (q=0)
      Call getevecfv0 (vkl0(1, ik), vgkl0(1, 1, 1, ik), evecfv0)
      Call getevalsv0 (vkl0(1, ik), evalsv0(1, ik))
  ! read occupation numbers for G+k
      Call getoccsv0 (vkl0(1, ik), occsv0(1, ik))
      evecfvo0 (:, :) = evecfv0 (ngk0(1, ik)+1:ngk0(1, ik)+nlotot, &
     & istl1:istu1, 1)
!      evecfvo20 (:, :) = evecfv0 (1:ngk0(1, ik), istl1:istu1, 1)
  ! change back file extension
!
      apwmaxsize=0
      do is=1,nspecies
        naug=0
        do l=0, input%xs%lmaxapwwf
          naug=naug+(2*l+1)*apword(l,is)
        enddo
        apwmaxsize=max(apwmaxsize,naug)
      enddo

      allocate(losize(is))
      lomaxsize=0
      losize=0
      do is=1,nspecies
        Do ilo = 1, nlorb (is)
          losize(is)=losize(is)+2*lorbl (ilo, is)+1
        enddo
        lomaxsize=max(lomaxsize,losize(is))
      enddo

      Call getapwcmt (0, ik, 1, nstfv, input%xs%lmaxapwwf, apwcmt0)
      allocate(apwcmtfun0(nst1,apwmaxsize,natmtot))
      apwcmtfun0=zzero
      do is=1,nspecies
        do ia=1,natoms(is)
          naug=0
          ias = idxas (ia, is)
          do l=0,input%xs%lmaxapwwf
            do m=-l,l
              do io=1,apword(l,is)
                naug=naug+1
                lm=idxlm(l,m)
                apwcmtfun0(1:nst1,naug,ias)=conjg(apwcmt0(istl1:istu1,io,lm,ias))
              enddo
            enddo
          enddo
        enddo
      enddo

      Call getapwcmt (iq, ikq, 1, nstfv, input%xs%lmaxapwwf, apwcmt)
      allocate(apwcmtfun(nst2,apwmaxsize,natmtot))
      apwcmtfun=zzero
      do is=1,nspecies
        do ia=1,natoms(is)
          naug=0
          ias = idxas (ia, is)
          do l=0,input%xs%lmaxapwwf
            do m=-l,l
              do io=1,apword(l,is)
                naug=naug+1
                lm=idxlm(l,m)
                apwcmtfun(1:nst2,naug,ias)=apwcmt(istl2:istu2,io,lm,ias)
              enddo
            enddo
          enddo
        enddo
      enddo

      locmt=zzero
!      Call getlocmt (0, ik, 1, nstfv, locmt0)
!      Call getlocmt (iq, ikq, 1, nstfv, locmt)
      allocate(locmtfun0(nst1,lomaxsize,natmtot))
      locmtfun0=zzero
      naug=0
      do is=1,nspecies
        do ia=1,natoms(is)
          ias=idxas(ia,is)
          locmtfun0(1:nst1,1:losize(is),ias)=conjg(transpose(evecfvo0 (naug+1:naug+losize(is), 1:nst1)))
          naug=naug+losize(is)
        enddo
      enddo

      allocate(locmtfun(nst2,lomaxsize,natmtot))
      locmtfun=zzero
      naug=0
      do is=1,nspecies
        do ia=1,natoms(is)
          ias=idxas(ia,is)
          locmtfun(1:nst2,1:losize(is),ias)=transpose(evecfvu (naug+1:naug+losize(is), 1:nst2))
          naug=naug+losize(is)
        enddo
      enddo
      

      
      Call timesec (cpu0)
      cpuread = cpu0 - cpu1
!
  ! zero matrix elements array
      xiou (:, :, :) = zzero
!
      whichthread=0
  ! loop over G+q vectors
cpugntlocal=0d0
cpumtlocal=0d0
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(igq,integrals,cpu00,cpu01,whichthread)
!print *,whichthread
#endif
      Allocate (integrals%aa(apwmaxsize,apwmaxsize, natmtot))
      Allocate (integrals%alo(lomaxsize, apwmaxsize, natmtot))
      Allocate (integrals%loa(apwmaxsize, lomaxsize, natmtot))
      Allocate (integrals%lolo(lomaxsize, lomaxsize, natmtot))
#ifdef USEOMP
      whichthread=omp_get_thread_num()
!$OMP DO
#endif

      Do igq = 1, ngq (iq)
         Call timesec (cpu00)
     ! summation of Gaunt coefficients wrt radial integrals
         Call ematgntsum (iq, igq,integrals)
         Call timesec (cpu01)
         if (whichthread.eq.0) cpugnt = cpugnt + cpu01 - cpu00

     ! muffin-tin contribution
         Call ematqkgmt (iq, ik, igq,integrals)
         Call timesec (cpu00)
         if (whichthread.eq.0) cpumt = cpumt + cpu00 - cpu01
!
      End Do ! igq
#ifdef USEOMP
!$OMP END DO
#endif
         
      deallocate(integrals%aa,integrals%alo,integrals%loa,integrals%lolo)
     
#ifdef USEOMP
!$OMP END PARALLEL
#endif

      deallocate (apwcmtfun,apwcmtfun0,losize,locmtfun,locmtfun0)
      Deallocate (evecfvu, evecfvo0)

      Allocate (evecfvo20(n0, nst1))
      Allocate (evecfvu2(n, nst2))
      Allocate (xihir(n0, n))
      evecfvu2 (:, :) = evecfv (1:ngk(1, ikq), istl2:istu2, 1)
      evecfvo20 (:, :) = evecfv0 (1:ngk0(1, ik), istl1:istu1, 1)


      Do igq = 1, ngq (iq)
         Call terminateqry ('ematqk')
         Call timesec (cpu00)
     ! interstitial contribution
         Call ematqkgir (iq, ik, igq)
         Call timesec (cpu01)
         cpuir = cpuir + cpu01 - cpu00
!
     ! interstitial contribution
         Call doublesummation_simple_cz (xiou(:, :, igq), evecfvo20, &
        & xihir, evecfvu2, zone, zone, .True.)
!
         Call timesec (cpu00)
         cpumirres = cpumirres + cpu00 - cpu01
!
         Call timesec (cpu01)
         cpudbg = cpudbg + cpu01 - cpu00
      End Do ! igq



!
      Call timesec (cpu1)
      cpumain = cpu1 - cpu0
!
  ! deallocate
      Deallocate (xihir)
!      Deallocate (evecfvu, evecfvo0)
      Deallocate (evecfvu2, evecfvo20)
      Call timesec (cpu0)
      cpuwrite = cpu0 - cpu1
      cpuall = cpuini + cpuread + cpumain + cpuwrite
! write(*,*)
!
  ! write timing information
      If ((task .Ne. 441) &
     & .And. (task .Ne. 450) .And. (task .Ne. 451)) Then
         Call emattim (iq, ik, trim(fnetim), cpuini, cpuread, cpumain, &
        & cpuwrite, cpuall, cpugnt, cpumt, cpuir, cpumalores, &
        & cpumloares, cpumlolores, cpumirres, cpudbg, cpumtaa, &
        & cpumtalo, cpumtloa, cpumtlolo)
      End If
!
End Subroutine ematqk
