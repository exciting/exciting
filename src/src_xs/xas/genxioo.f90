!
!
!
! Copyright (C) 2015 Christian Vorwerk and Claudia Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
!BOP
!
! !ROUTINE: genxioo
!
! !INTERFACE:
Subroutine genxioo (iq, ik, xi)
! !USES:
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
      use modxas
#ifdef USEOMP
      use omp_lib
#endif
! !INPUT/OUTPUT PARAMETERS:
!   ik  : k-point position (in,integer)
!   iq  : q-point position (in,integer)
!   xi  : planewave matrix element (complex(nxas, nxas, ngq(iq)), out)
! !DESCRIPTION:
!   Calculates Planewave Matrix elements $$M_{\mu \mu'}(\mathbf{q}+\mathbf{G})$$ between 
!   two core states. See Vorwerk's Master thesis for more details. 
!
! !REVISION HISTORY:
!  Based on the subroutine ematqk.F90
!  Created November 2015 (Christian Vorwerk)
!EOP
!BOC      

      Implicit None
!      type(mtints_type) :: integrals
      
  ! arguments
      Integer, Intent (InOut) :: iq, ik
      Complex(8), Intent(Out) :: xi(nxas, nxas, ngq(iq))
  ! local variables
      Character (*), Parameter :: thisnam = 'ematqk'
  ! allocatable arrays
      Complex (8), Allocatable :: evecfvo0 (:, :)
      Complex (8), Allocatable :: evecfvu (:, :)
      Complex (8), Allocatable :: evecfvo20 (:, :)
      Complex (8), Allocatable :: evecfvu2 (:, :)
      Complex (8), Allocatable :: zfft0(:,:),zfft(:),zfftres(:),zfftcf(:)
      Complex (8), Allocatable :: apwalmt (:, :, :, :), apwalmt0 (:, :, :, :)
      Complex(8), allocatable :: wfmt(:,:,:), wfmt0(:,:,:)
      Complex (8), Allocatable :: xihir (:, :)
  ! expansion coefficients of APW and LO functions
      Complex (8), Allocatable :: integral(:,:,:)

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
  ! read eigenvectors, eigenvalues and occupancies for G+k (q=0)
      Call getevecfv0 (vkl0(1, ik), vgkl0(1, 1, 1, ik), evecfv0)
      Call getevalsv0 (vkl0(1, ik), evalsv0(1, ik))
  ! read occupation numbers for G+k
      Call getoccsv0 (vkl0(1, ik), occsv0(1, ik))
      evecfvo0 (:, :) = evecfv0 (ngk0(1, ik)+1:ngk0(1, ik)+nlotot, &
     & istl1:istu1, 1)

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
		 Allocate (integral(input%xs%lmaxemat+1,nxas,nxas))
         Call ematradoo (iq, ik, igq,integral)
         Call timesec (cpu01)
         if (whichthread.eq.0) cpugnt = cpugnt + cpu01 - cpu00
     ! muffin-tin contribution
         Call ematsumoo (iq, ik, igq,integral, xi)
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
      Deallocate (evecfvu, evecfvo0)




      Allocate (evecfvo20(n0, nst1))
      Allocate (evecfvu2(n, nst2))
      evecfvu2 (:, :) = evecfv (1:ngk(1, ikq), istl2:istu2, 1)
      evecfvo20 (:, :) = evecfv0 (1:ngk0(1, ik), istl1:istu1, 1)

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
End Subroutine genxioo
! EOC
