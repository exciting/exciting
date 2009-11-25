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
      Implicit None
  ! arguments
      Integer, Intent (In) :: iq, ik
  ! local variables
      Character (*), Parameter :: thisnam = 'ematqk'
  ! allocatable arrays
      Complex (8), Allocatable :: evecfvo0 (:, :)
      Complex (8), Allocatable :: evecfvu (:, :)
      Complex (8), Allocatable :: evecfvo20 (:, :)
      Complex (8), Allocatable :: evecfvu2 (:, :)
      Integer :: ikq, igq, n, n0
      Real (8) :: cpuini, cpuread, cpumain, cpuwrite, cpuall
      Real (8) :: cpugnt, cpumt, cpuir
      Real (8) :: cpumalores, cpumloares
      Real (8) :: cpumlolores, cpumirres, cpudbg
      Real (8) :: cpu0, cpu1, cpu00, cpu01
!
      If (task .Eq. 330) Call chkpt (3, (/ task, iq, ik /), 'ematqk: ta&
     &sk, q - point index, k - point index; q - dependent matrix elemen&
     &ts')
      Call cpu_time (cpu0)
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
      Allocate (evecfvo20(n0, nst1))
      Allocate (evecfvu2(n, nst2))
      Allocate (xihir(n0, n))
  ! zero arrays
      xiohalo (:, :) = zzero
      xiuhloa (:, :) = zzero
!
      Call cpu_time (cpu1)
      cpuini = cpu1 - cpu0
!
  ! read eigenvectors, eigenvalues and occupancies for G+k+q
      Call getevecfv (vkl(1, ikq), vgkl(1, 1, 1, ikq), evecfv)
      Call getevalsv (vkl(1, ikq), evalsv(1, ikq))
  ! read occupation numbers for G+k+q
      Call getoccsv (vkl(1, ikq), occsv(1, ikq))
      evecfvu (:, :) = evecfv (ngk(1, ikq)+1:ngk(1, ikq)+nlotot, &
     & istl2:istu2, 1)
      evecfvu2 (:, :) = evecfv (1:ngk(1, ikq), istl2:istu2, 1)
!
  ! read eigenvectors, eigenvalues and occupancies for G+k (q=0)
      Call getevecfv0 (vkl0(1, ik), vgkl0(1, 1, 1, ik), evecfv0)
      Call getevalsv0 (vkl0(1, ik), evalsv0(1, ik))
  ! read occupation numbers for G+k
      Call getoccsv0 (vkl0(1, ik), occsv0(1, ik))
      evecfvo0 (:, :) = evecfv0 (ngk0(1, ik)+1:ngk0(1, ik)+nlotot, &
     & istl1:istu1, 1)
      evecfvo20 (:, :) = evecfv0 (1:ngk0(1, ik), istl1:istu1, 1)
  ! change back file extension
!
      Call getapwcmt (0, ik, 1, nstfv, input%xs%lmaxapwwf, apwcmt0)
      Call getapwcmt (iq, ikq, 1, nstfv, input%xs%lmaxapwwf, apwcmt)
      Call getlocmt (0, ik, 1, nstfv, locmt0)
      Call getlocmt (iq, ikq, 1, nstfv, locmt)
!
      Call cpu_time (cpu0)
      cpuread = cpu0 - cpu1
!
  ! zero matrix elements array
      xiou (:, :, :) = zzero
!
  ! loop over G+q vectors
      Do igq = 1, ngq (iq)
         Call terminateqry ('ematqk')
         Call cpu_time (cpu00)
     ! summation of Gaunt coefficients wrt radial integrals
         Call ematgntsum (iq, igq)
         Call cpu_time (cpu01)
         cpugnt = cpugnt + cpu01 - cpu00
     ! muffin-tin contribution
         Call ematqkgmt (iq, ik, igq)
         Call cpu_time (cpu00)
         cpumt = cpumt + cpu00 - cpu01
     ! interstitial contribution
         Call ematqkgir (iq, ik, igq)
         Call cpu_time (cpu01)
         cpuir = cpuir + cpu01 - cpu00
!
         If (( .Not. input%xs%fastemat) .And. (nlotot .Gt. 0)) Then
        ! muffin-tin contributions
        ! APW-lo contribution
        ! multiplication xi = xiho * evecfvu
            Call zgemm ('n', 'n', nst1, nst2, nlotot, zone, xiohalo, &
           & nst1, evecfvu, nlotot, zone, xiou(1, 1, igq), nst1)
            Call cpu_time (cpu00)
            cpumalores = cpumalores + cpu00 - cpu01
        ! lo-APW contribution
        ! multiplication xi = evecfvo * xihu
            Call zgemm ('c', 'n', nst1, nst2, nlotot, zone, evecfvo0, &
           & nlotot, xiuhloa, nlotot, zone, xiou(1, 1, igq), nst1)
            Call cpu_time (cpu01)
            cpumloares = cpumloares + cpu01 - cpu00
        ! lo-lo contribution
            Call doublesummation_simple_cz (xiou(:, :, igq), evecfvo0, &
           & xih, evecfvu, zone, zone, .True.)
!
            Call cpu_time (cpu00)
            cpumlolores = cpumlolores + cpu00 - cpu01
            cpu01 = cpu00
         End If
!
     ! interstitial contribution
         Call doublesummation_simple_cz (xiou(:, :, igq), evecfvo20, &
        & xihir, evecfvu2, zone, zone, .True.)
!
         Call cpu_time (cpu00)
         cpumirres = cpumirres + cpu00 - cpu01
!
         Call cpu_time (cpu01)
         cpudbg = cpudbg + cpu01 - cpu00
      End Do ! igq
!
      Call cpu_time (cpu1)
      cpumain = cpu1 - cpu0
!
  ! deallocate
      Deallocate (xihir)
      Deallocate (evecfvu, evecfvo0)
      Deallocate (evecfvu2, evecfvo20)
      Call cpu_time (cpu0)
      cpuwrite = cpu0 - cpu1
      cpuall = cpuini + cpuread + cpumain + cpuwrite
!
  ! write timing information
      If ((task .Ne. 430) .And. (task .Ne. 440) .And. (task .Ne. 441) &
     & .And. (task .Ne. 450) .And. (task .Ne. 451)) Then
         Call emattim (iq, ik, trim(fnetim), cpuini, cpuread, cpumain, &
        & cpuwrite, cpuall, cpugnt, cpumt, cpuir, cpumalores, &
        & cpumloares, cpumlolores, cpumirres, cpudbg, cpumtaa, &
        & cpumtalo, cpumtloa, cpumtlolo)
      End If
!
End Subroutine ematqk
