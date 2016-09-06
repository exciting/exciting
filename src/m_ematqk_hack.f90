! Copyright(C) 2006-2008 S. Sagmeister and Claudia Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
module m_ematqk_hack
  use mod_kpointset
  use mod_atoms
  use mod_APW_LO
  use modinput, only: input

  !! TO BE INITIALIZED
  ! APW matching coefficients for k and k+q point
  ! lmaxapwwf  : maximum l over apws

  ! Integrate here
      !ematrad
      ! apwfr and lofr need to be precalculated
      !ematqk
      use modxs, only: ikmapikq, msg, ngk0, xiohalo,&
                     & xiuhloa, evecfv,&
                     & vkl0, vgkl0, evecfv0, evalsv0,&
                     & occsv0, apwmaxsize, apwsize, losize,&
                     & lomaxsize, cmtfun0, cmtfun, apwcmt0,&
                     & apwcmt, ngq, igqig,&
                     & fnetim, fftmap_type, igkig0
      ! ematgntsum
      use modxs, only: apwmaxsize, apwsize, lomaxsize, ylmgq,&
                     & xsgnt
      ! input%xs%lmaxemat ! Defines the l cutoff of e^-i(G+q)r

      ! ematqkgmt
      ! sfacgq needs to be supplied

      ! Substitute with mod_kpointset quantities
      use mod_Gkvector, only: ngk, vgkl, igkig, gkmax
      use mod_Gvector, only: gc, ngvec, cfunig, ivg, ivgig
      use mod_kpoint, only: vkl
      use mod_qpoint, only: vql

  implicit none

  private

  ! Cutoff l for APWs
  integer(4) :: lmaxapw
  ! Cutoff l for Rayleigh expansion of plane waves
  integer(4) :: lrayleigh 
  ! radial integrals (apw-apw)
  real(8), allocatable :: riaa(:, :, :, :, :, :, :)
  ! radial integrals (lo-apw)
  real(8), allocatable :: riloa(:, :, :, :, :, :)
  ! radial integrals (lo-lo)
  real(8), allocatable :: rilolo(:, :, :, :, :)

  ! Our lovely grids
  type(kq_set) :: kqgrid
  type(Gk_set) :: Gqgrid

  public :: ematqk_setup, ematdummy

  contains

    subroutine ematqk_setup(kqg, gqg, lray)
      implicit none
      type(kq_set), intent(in) :: kqg
      type(Gk_set), intent(in) :: gqg
      integer(4), intent(in), optional :: lray

      ! Set module grids
      kqgrid = kqg
      Gqgrid = gqg

      ! Generate radial functions
      call genapwfr
      call genlofr

      ! Check maximal l for apws
      if(associated(input%xs)) then
        lmaxapw = input%xs%lmaxapwwf
      else
        lmaxapw = input%groundstate%lmaxapw
      end if

      ! Set maximal l for Rayleigh expansion
      if(present(lray)) then
        lrayleigh = lray
      else if(associated(input%xs)) then
        lrayleigh = input%xs%lmaxemat
      else
        lrayleigh = 3
      end if

    end subroutine ematqk_setup
    
    subroutine ematdummy(iq)
      integer(4), intent(in) :: iq

      integer(4) :: riaashape(7), i,j,k,l

      write(*,*) "Before"
      write(*,*) shape(riaa)
      write(*,*) shape(riloa)
      write(*,*) shape(rilolo)
      call ematrad_hack(iq)
      write(*,*) "After"
      write(*,*) shape(riaa)
      write(*,*) shape(riloa)
      write(*,*) shape(rilolo)

      write(*,'(E10.3)') riaa(:,1,1,1,1,1,1)

    end subroutine ematdummy

!    !BOP
!    ! !ROUTINE: ematqk
!    ! !INTERFACE:
!    subroutine ematqk_hack(iq, ik, emat, bc)
!    ! !USES:
!      use mod_misc, only: task
!      use mod_constants, only: zzero, zone
!      use mod_eigenvalue_occupancy, only: evalsv, occsv, nstfv
!      use mod_muffin_tin, only: idxlm
!      use mod_atoms, only: nspecies, natmtot, natoms, idxas
!      use summations, only: doublesummation_simple_cz
!      use m_getapwcmt
!      use m_getlocmt
!      use m_putemat
!      use m_emattim
!      use m_getunit
!      use m_genfilname
!#ifdef USEOMP
!      use omp_lib
!#endif
!    ! !DESCRIPTION:
!    ! Calculates plane wave elements between state ranges defined
!    ! in (bc%n1,bc%il1,bc%iu1) and (bc%n2,bc%il2,bc%iu2) and saves them in 
!    ! emat.
!    !
!    ! !REVISION HISTORY:
!    !   Added to documentation scheme. (Aurich) \\
!    !   Substituted the modxs:xiou reference with the pointer modxs:emat that
!    !   can point to any of the modxs:xiXY.
!    !EOP
!    !BOC
!
!      implicit none
!          
!      ! Arguments
!      integer, intent(in) :: iq, ik
!      type(bcbs), intent(in) :: bc
!      complex(8), intent(inout) :: emat(:,:,:)
!
!      ! Local variables
!      character(*), parameter :: thisnam = 'ematqk_hack'
!      ! Allocatable arrays
!      complex(8), allocatable :: evecfvo0(:, :)
!      complex(8), allocatable :: evecfvu(:, :)
!      complex(8), allocatable :: evecfvo20(:, :)
!      complex(8), allocatable :: evecfvu2(:, :)
!      complex(8), allocatable :: zfft0(:, :), zfft(:), zfftres(:), zfftcf(:)
!      complex(8), allocatable :: xihir(:, :)
!      ! Expansion coefficients of apw and lo functions
!      complex(8), allocatable :: integrals(:, :, :)
!      integer :: ikq, igq, n, n0, is, l, m, io, naug, ia, ias, lm, ilo
!      integer :: whichthread, ig, ix, igs, ist2, ist1
!      real(8) :: cpuini, cpuread, cpumain, cpuwrite, cpuall
!      real(8) :: cpugnt, cpumt, cpuir, cpufft
!      real(8) :: cpumalores, cpumloares
!      real(8) :: cpumlolores, cpumirres, cpudbg
!      real(8) :: cpu0, cpu1, cpu00, cpu01, cpugntlocal, cpumtlocal
!      real(8) :: vkql(3)
!      integer :: shift(3), iv(3)
!      type(fftmap_type) :: fftmap
!      real(8) :: emat_gmax
!
!     ! ! If task 330 is 'writeemat'
!     ! if(task .eq. 330) then
!     !   call chkpt(3, (/ task, iq, ik /),&
!     !     & 'ematqk: task, q - point index, k - point index; q - dependent matrix elements')
!     ! end if
!
!      ! Find k+q-point
!      ikq = ikmapikq(ik, iq)
!
!      ! Get number of G+k and G+k+q vectors
!      n0 = ngk0(1, ik)
!      n = ngk(1, ikq)
!
!      ! Allocate matrix elements array
!      if(allocated(xiohalo)) deallocate(xiohalo)
!      allocate(xiohalo(bc%n1, nlotot))
!      if(allocated(xiuhloa)) deallocate(xiuhloa)
!      allocate(xiuhloa(nlotot, bc%n2))
!
!      ! Allocate temporary arrays
!      allocate(evecfvo0(nlotot, bc%n1))
!      allocate(evecfvu(nlotot, bc%n2))
!
!      ! Zero arrays
!      xiohalo(:, :) = zzero
!      xiuhloa(:, :) = zzero
!
!      ! Read eigenvectors, eigenvalues and occupancies for k+q
!      call getevecfv(vkl(1, ikq), vgkl(1, 1, 1, ikq), evecfv)
!      ! Get local orbital coefficients for range 2 and spin 1 
!      evecfvu(:, :) = evecfv(ngk(1, ikq)+1:ngk(1, ikq)+nlotot, bc%il2:bc%iu2, 1)
!
!      ! Read eigenvectors, eigenvalues and occupancies for k (q=0)
!      call getevecfv0(vkl0(1, ik), vgkl0(1, 1, 1, ik), evecfv0)
!      ! Get local orbital coefficients for range 2 and spin 1 
!      evecfvo0(:, :) = evecfv0(ngk0(1, ik)+1:ngk0(1, ik)+nlotot, bc%il1:bc%iu1, 1)
!
!      ! Determine number of radial functions used in APW 
!      ! basis functions per species
!      apwmaxsize=0
!      allocate(apwsize(nspecies))
!      do is=1, nspecies
!        naug=0
!        do l=0, lmaxapw
!          naug=naug+(2*l+1)*apword(l, is)
!        end do
!        apwsize(is)=naug
!        apwmaxsize=max(apwmaxsize, naug)
!      end do
!
!      ! Determine number of radial functions used in 
!      ! LO basis functions per species
!      allocate(losize(nspecies))
!      lomaxsize=0
!      losize=0
!      do is=1, nspecies
!        do ilo = 1, nlorb(is)
!          losize(is)=losize(is)+2*lorbl(ilo, is)+1
!        end do
!        lomaxsize=max(lomaxsize, losize(is))
!      end do
!
!!!! Combined array ??
!      allocate(cmtfun0(bc%n1, apwmaxsize+lomaxsize, natmtot))
!      cmtfun0=zzero
!      allocate(cmtfun(bc%n2, apwmaxsize+lomaxsize, natmtot))
!      cmtfun=zzero
!      apwcmt0=zzero
!      apwcmt=zzero
!
!      ! Get APW expansion coefficients for k
!      call getapwcmt(0, ik, 1, nstfv, lmaxapw, apwcmt0)
!      ilo=0
!      do is=1, nspecies
!        do ia=1, natoms(is)
!          naug=0
!          ! Combined mufftin-species index (alpha)
!          ias = idxas(ia, is)
!
!          do l=0, lmaxapw
!            do m=-l, l
!              do io=1, apword(l, is)
!                naug=naug+1
!                lm=idxlm(l, m)
!                ! Save apw matching coefficients
!                cmtfun0(1:bc%n1, naug, ias)=conjg(apwcmt0(bc%il1:bc%iu1, io, lm, ias))
!              end do
!            end do
!          end do
!          ! Save local orbital coefficients 
!          cmtfun0(1:bc%n1, naug+1:naug+losize(is), ias)=&
!            & conjg(transpose(evecfvo0(ilo+1:ilo+losize(is), 1:bc%n1)))
!          ilo=ilo+losize(is)
!
!        end do
!      end do
!
!      ! Get APW expansion coefficients for k+q
!      call getapwcmt(iq, ikq, 1, nstfv, lmaxapw, apwcmt)
!      ilo=0
!      do is=1, nspecies
!        do ia=1, natoms(is)
!          naug=0
!          ias = idxas(ia, is)
!          do l=0, lmaxapw
!            do m=-l, l
!              do io=1, apword(l, is)
!                naug=naug+1
!                lm=idxlm(l, m)
!                cmtfun(1:bc%n2, naug, ias)=apwcmt(bc%il2:bc%iu2, io, lm, ias)
!              end do
!            end do
!          end do
!          cmtfun(1:bc%n2, naug+1:naug+losize(is), ias)=&
!            & transpose(evecfvu(ilo+1:ilo+losize(is), 1:bc%n2))
!          ilo=ilo+losize(is)
!        end do
!      end do
!
!      ! Zero matrix elements array
!      emat(:, :, :) = zzero
!
!      whichthread=0
!
!#ifdef USEOMP
!    !$omp parallel default(shared) private(igq, integrals, cpu00, cpu01, whichthread)
!#endif
!      ! Allocate radial integrals
!      allocate(integrals(apwmaxsize+lomaxsize, apwmaxsize+lomaxsize, natmtot))
!#ifdef USEOMP
!      whichthread=omp_get_thread_num()
!    !$omp do
!#endif
!      do igq = 1, Gqgrid%ngk(1, iq)
!        ! Summation of gaunt coefficients w.r.t. radial integrals
!        call ematgntsum_hack(iq, igq, integrals)
!
!        ! Muffin-tin contribution
!        call ematqkgmt_hack(iq, ik, igq, integrals, emat(:,:,igq), bc)
!
!      end do ! igq
!#ifdef USEOMP
!    !$omp end do
!#endif
!      deallocate(integrals)        
!         
!#ifdef USEOMP
!    !$omp end parallel
!#endif
!      deallocate(apwsize, losize, cmtfun, cmtfun0)
!      deallocate(evecfvu, evecfvo0)
!     
!      allocate(evecfvo20(n0, bc%n1))
!      allocate(evecfvu2(n, bc%n2))
!      evecfvu2(:, :) = evecfv(1:ngk(1, ikq), bc%il2:bc%iu2, 1)
!      evecfvo20(:, :) = evecfv0(1:ngk0(1, ik), bc%il1:bc%iu1, 1)
!     
!      ! Interstitial contribution
!      ! Matrix-matrix multiplication
!#ifdef USEOMP
!    !$omp parallel default(shared) private(igq, xihir, cpu00, cpu01, whichthread)
!        whichthread=omp_get_thread_num()
!#endif 
!        allocate(xihir(n0, n))
!#ifdef USEOMP
!    !$omp do
!#endif
!        do igq = 1, Gqgrid%ngk(1, iq)
!          ! Interstitial contribution
!          call ematqkgir_hack(iq, ik, igq, xihir, n0, n)
!
!          ! Interstitial contribution
!          call doublesummation_simple_cz(emat(:, :, igq), evecfvo20,&
!            & xihir, evecfvu2, zone, zone, .true.)
!
!        end do ! Igq
!#ifdef USEOMP
!    !$omp end do
!#endif
!        deallocate(xihir)
!#ifdef USEOMP
!    !$omp end parallel
!#endif
!
!    end subroutine ematqk_hack
!    !EOC
!
!    subroutine ematqkgmt_hack(iq, ik, igq, integrals, emat, bc)
!      use modinput, only: input
!      use mod_constants, only: zzero, zone, fourpi
!      use mod_atoms, only: natmtot, nspecies, natoms, idxas
!      use modxs, only: bcbs, apwmaxsize, lomaxsize, &
!                     & sfacgq,&
!                     & apwsize, losize, cmtfun, cmtfun0,&
!                     & cpumtaa
!#ifdef USEOMP
!    use omp_lib
!#endif
!
!      implicit none
!
!      ! Arguments
!      integer, intent(in) :: iq, ik, igq
!      type(bcbs), intent(in) :: bc
!      complex(8), intent(inout) :: emat(:,:)
!
!      complex(8) :: integrals(apwmaxsize+lomaxsize,apwmaxsize+lomaxsize,natmtot)
!
!      ! Local variables
!      character(*), parameter :: thisnam = 'ematqkgmt'
!      integer :: is, ia, ias
!      integer :: lmax1, lmax3, ikt, zmsize, whichthread
!      complex(8), allocatable :: zm(:,:)
!      complex(8) :: prefactor
!      real(8) :: cmt0, cmt1
!
!#ifdef USEOMP
!      whichthread=omp_get_thread_num()
!#else
!      whichthread=0
!#endif
!
!      ikt = ik
!      lmax1 = lmaxapw
!      lmax3 = lmax1
!      zmsize=apwmaxsize+lomaxsize
!
!      allocate(zm(1:bc%n2,zmsize))
!
!      emat(:, :) = zzero
!
!      ! Loop over species and atoms
!      do is = 1, nspecies
!        do ia = 1, natoms(is)
!
!          ias = idxas(ia, is)
!
!          !---------------------------!
!          !     apw-apw contribution  !
!          !---------------------------!
!          prefactor=fourpi*conjg(sfacgq(igq, ias, iq))
!
!          call zgemm('n', 'n', bc%n2, apwsize(is)+losize(is), apwsize(is)+losize(is),&
!            & zone, cmtfun(1,1,ias), bc%n2, integrals(1,1,ias), apwmaxsize+lomaxsize,&
!            & zzero, zm, bc%n2)
!          call zgemm('n', 't', bc%n1, bc%n2, apwsize(is)+losize(is),&
!            & prefactor, cmtfun0(1,1,ias), bc%n1, zm, bc%n2, zone, emat, bc%n1)
!
!          call timesec(cmt1)
!          if(whichthread.eq.0) then
!            cpumtaa = cpumtaa + cmt1 - cmt0
!          endif
!
!        ! End loop over species and atoms
!        end do ! ia
!      end do ! is
!
!      deallocate(zm)
!    end subroutine ematqkgmt_hack
!
!    subroutine ematgntsum_hack(iq, igq, integrals)
!      use modinput, only: input
!      use mod_misc, only: filext
!      use mod_constants, only: zzero, zil
!      use mod_muffin_tin, only: idxlm
!      use mod_atoms, only: natmtot, nspecies, natoms, idxas
!      use mod_APW_LO, only: lolmax, apwordmax, nlomax, apword,&
!                          & nlorb, lorbl
!      use m_findgntn0, only: l1map, l2map, l3map,&
!                           & m1map, m2map, m3map,&
!                           & l1shape, l2shape, l3shape,&
!                           & m1shape, m2shape, m3shape
!      use m_getunit
!
!      implicit none
!
!      ! Arguments
!      integer, intent(in) :: iq, igq
!      complex(8), intent(out) :: integrals(apwmaxsize+lomaxsize, apwmaxsize+lomaxsize, natmtot)
!
!      ! Local variables
!      integer :: is, ia, ias, iaug1, iaug2
!      integer :: l1, l2, l3, m2, lm2
!      integer :: ilo, ilo1, ilo2, io, io1, io2
!      integer :: lmax1, lmax2, lmax3, lmmax1, lmmax2, lmmax3
!      integer :: u1, u2, u3, u4
!      integer :: m1, m3, lm1, lm3, cl1, cm1, cl2, cm2, cl3, cm3
!      complex(8), dimension(:,:,:,:), allocatable :: intrgaa, intrgloa, intrglolo, intrgalo
!
!      ! Set lm related local variables
!      lmax1 = max(lmaxapw, lolmax)
!      lmax2 = lrayleigh
!      ! lmax1 and lmax3 should be the same
!      lmax3 = lmax1
!      lmmax1 = (lmax1+1) ** 2
!      lmmax2 = (lmax2+1) ** 2
!      lmmax3 = (lmax3+1) ** 2
!
!      ! Allocate arrays for radial integrals and Bessel functions
!      allocate(intrgaa(lmmax1, apwordmax, lmmax3, apwordmax))
!      allocate(intrgloa(-lolmax:lolmax, nlomax, lmmax3, apwordmax))
!      allocate(intrgalo(-lolmax:lolmax, nlomax, lmmax3, apwordmax))
!      allocate(intrglolo(-lolmax:lolmax, nlomax, -lolmax:lolmax, nlomax))
!
!      integrals=zzero
!
!      ! Begin loop over species
!      do is = 1, nspecies
!        ! Begin loop over atoms
!        do ia = 1, natoms(is)
!
!          intrgaa(:, :, :, :) = zzero
!          intrgloa(:, :, :, :) = zzero
!          intrgalo(:, :, :, :) = zzero
!          intrglolo(:, :, :, :) = zzero
!
!          ias = idxas(ia, is)
!
!          !---------------------------!
!          !     apw-apw integrals     !
!          !---------------------------!
!          do cl1 = 1, l1shape
!            l1 = l1map(cl1)
!            do cm1 = 1, m1shape(l1)
!              m1 = m1map(l1, cm1)
!              lm1 = idxlm(l1, m1)
!              do io1 = 1, apword(l1, is)
!
!                do cl2 = 1, l2shape(l1, m1)
!                  l3 = l2map(l1, m1, cl2)
!                  do cm2 = 1, m2shape(l1, m1, l3)
!                    m3 = m2map(l1, m1, l3, cm2)
!                    lm3 = idxlm(l3, m3)
!                    do io2 = 1, apword(l3, is)
!
!                      do cl3 = 1, l3shape(l1, m1, l3, m3)
!                        l2 = l3map(l1, m1, l3, m3, cl3)
!                        do cm3 = 1, m3shape(l1, m1, l3, m3, l2)
!                          m2 = m3map(l1, m1, l3, m3, l2, cm3)
!                          lm2 = idxlm(l2, m2)
!
!                          ! Sagmeister (8.26)
!                          intrgaa(lm1, io1, lm3, io2) =&
!                            & intrgaa(lm1, io1, lm3, io2)&
!                            &+ conjg(zil(l2)) * riaa(l1, io1, l3, io2, l2, ias, igq)&
!                            &* conjg(ylmgq(lm2, igq, iq)) * xsgnt(lm1, lm2, lm3)
!                        end do
!                      end do
!
!                    end do
!                  end do
!                end do
!
!              end do
!            end do
!          end do
!
!          iaug2=0
!          do l3=0, lmaxapw
!            do m3=-l3, l3
!              do io2=1, apword(l3, is)
!                iaug2=iaug2+1
!                lm3=idxlm(l3, m3)
!
!                iaug1=0
!                do l1=0, lmaxapw
!                  do m1=-l1, l1
!                     do io1=1, apword(l1, is)
!                        iaug1=iaug1+1
!                        lm1=idxlm(l1, m1)
!
!                        integrals(iaug2, iaug1, ias)=intrgaa(lm1, io1, lm3, io2)
!
!                      end do
!                    end do
!                  end do
!
!              end do
!            end do
!          end do
!
!          !-------------------------------------!
!          !     apw-local-orbital integrals     !
!          !-------------------------------------!
!          do cl1 = 1, l1shape
!            l1 = l1map(cl1)
!            do cm1 = 1, m1shape(l1)
!              m1 = m1map(l1, cm1)
!              lm1 = idxlm(l1, m1)
!              do io = 1, apword(l1, is)
!
!                do ilo = 1, nlorb(is)
!                  l3 = lorbl(ilo, is)
!                  do cm2 = 1, m2shape(l1, m1, l3)
!                    m3 = m2map(l1, m1, l3, cm2)
!                    lm3 = idxlm(l3, m3)
!
!                    do cl3 = 1, l3shape(l1, m1, l3, m3)
!                      l2 = l3map(l1, m1, l3, m3, cl3)
!                      do cm3 = 1, m3shape(l1, m1, l3, m3, l2)
!                        m2 = m3map(l1, m1, l3, m3, l2, cm3)
!                        lm2 = idxlm(l2, m2)
!
!                        intrgalo(m3, ilo, lm1, io) = intrgalo(m3, ilo, lm1, io)&
!                          &+ conjg(zil(l2)) * riloa(ilo, l1, io, l2, ias, igq)&
!                          &* conjg(ylmgq(lm2, igq, iq)) * xsgnt(lm1, lm2, lm3)
!                      end do
!                    end do
!
!                  end do
!                end do
!
!              end do
!            end do
!          end do
!
!          iaug2=0
!          do ilo = 1, nlorb(is)
!            l3 = lorbl(ilo, is)
!            do m3=-l3, l3
!              iaug2=iaug2+1
!              lm3=idxlm(l3, m3)
!             
!              iaug1=0
!              do l1=0, lmaxapw
!                do m1=-l1, l1
!                  do io=1, apword(l1, is)
!                    iaug1=iaug1+1
!                    lm1=idxlm(l1, m1)
!
!                    integrals(apwsize(is)+iaug2, iaug1, ias)=intrgalo(m3, ilo, lm1, io)
!
!                  end do
!                end do
!               end do
!
!            end do
!          end do
!
!          !-------------------------------------!
!          !     local-orbital-apw integrals     !
!          !-------------------------------------!
!          do ilo = 1, nlorb(is)
!            l1 = lorbl(ilo, is)
!            do cm1 = 1, m1shape(l1)
!              m1 = m1map(l1, cm1)
!              lm1 = idxlm(l1, m1)
!
!              do cl2 = 1, l2shape(l1, m1)
!                l3 = l2map(l1, m1, cl2)
!                do cm2 = 1, m2shape(l1, m1, l3)
!                  m3 = m2map(l1, m1, l3, cm2)
!                  lm3 = idxlm(l3, m3)
!                  do io = 1, apword(l3, is)
!
!                    do cl3 = 1, l3shape(l1, m1, l3, m3)
!                      l2 = l3map(l1, m1, l3, m3, cl3)
!                      do cm3 = 1, m3shape(l1, m1, l3, m3, l2)
!                        m2 = m3map(l1, m1, l3, m3, l2, cm3)
!                        lm2 = idxlm(l2, m2)
!
!                        intrgloa(m1, ilo, lm3, io) = intrgloa(m1, ilo, lm3, io)&
!                          &+ conjg(zil(l2)) * riloa(ilo, l3, io, l2, ias, igq)&
!                          &* conjg(ylmgq(lm2, igq, iq)) * xsgnt(lm1, lm2, lm3)
!                      end do
!                    end do
!
!                  end do
!                end do
!              end do
!
!            end do
!          end do
!
!          iaug2=0
!          do ilo = 1, nlorb(is)
!            l3 = lorbl(ilo, is)
!            do m3=-l3, l3
!              iaug2=iaug2+1
!              lm3=idxlm(l3, m3)
!
!              iaug1=0
!              do l1=0, lmaxapw
!                do m1=-l1, l1
!                  do io=1, apword(l1, is)
!                    iaug1=iaug1+1
!                    lm1=idxlm(l1, m1)
!                    integrals(iaug1, apwsize(is)+iaug2, ias)=intrgloa(m3, ilo, lm1, io)
!                  end do
!                end do
!               end do
!
!            end do
!          end do
!
!          !-----------------------------------------------!
!          !     local-orbital-local-orbital integrals     !
!          !-----------------------------------------------!
!          do ilo1 = 1, nlorb(is)
!            l1 = lorbl(ilo1, is)
!            do cm1 = 1, m1shape(l1)
!              m1 = m1map(l1, cm1)
!              lm1 = idxlm(l1, m1)
!
!              do ilo2 = 1, nlorb(is)
!                l3 = lorbl(ilo2, is)
!                do cm2 = 1, m2shape(l1, m1, l3)
!                  m3 = m2map(l1, m1, l3, cm2)
!                  lm3 = idxlm(l3, m3)
!
!                  do cl3 = 1, l3shape(l1, m1, l3, m3)
!                    l2 = l3map(l1, m1, l3, m3, cl3)
!                    do cm3 = 1, m3shape(l1, m1, l3, m3, l2)
!                      m2 = m3map(l1, m1, l3, m3, l2, cm3)
!                      lm2 = idxlm(l2, m2)
!
!                      intrglolo(m1, ilo1, m3, ilo2) = intrglolo(m1, ilo1, m3, ilo2)&
!                        &+ conjg(zil(l2)) * rilolo(ilo1, ilo2, l2, ias, igq)&
!                        &* conjg(ylmgq(lm2, igq, iq)) * xsgnt(lm1, lm2, lm3)
!                    end do
!                  end do
!
!                end do
!              end do
!
!            end do
!          end do
!
!          iaug2=0
!          do ilo2 = 1, nlorb(is)
!            l3 = lorbl(ilo2, is)
!            do m3=-l3, l3
!              iaug2=iaug2+1
!
!              iaug1=0
!              do ilo1 = 1, nlorb(is)
!                l1 = lorbl(ilo1, is)
!                do m1=-l1, l1
!                  iaug1=iaug1+1
!                  integrals(apwsize(is)+iaug2, apwsize(is)+iaug1, ias)=&
!                    & intrglolo(m1, ilo1, m3, ilo2)
!                end do
!              end do
!
!            end do
!          end do
!
!        ! End loops over atoms and species
!        end do
!      end do
!
!      ! Deallocate
!      if(allocated(intrgaa)) deallocate(intrgaa)
!      if(allocated(intrgloa)) deallocate(intrgloa)
!      if(allocated(intrgalo)) deallocate(intrgalo)
!      if(allocated(intrglolo)) deallocate(intrglolo)
!
!    end subroutine ematgntsum_hack
!
!    subroutine ematqkgir_hack(iq, ik, igq, xihir, n0, n)
!      use mod_Gvector, only: ivg, ivgig, cfunig
!      use mod_Gkvector, only: ngkmax, igkig, ngk
!      use mod_qpoint, only: vql
!      use mod_kpoint, only: vkl
!      use modxs, only: ikmapikq, ngkmax0, vkl0, igkig0,&
!                     & igqig, ngk0
!
!      implicit none
!
!      ! Arguments
!      integer, intent(in) :: iq, ik, igq, n0,n
!      complex(8), intent(out) :: xihir(n0,n) 
!
!      ! Local variables
!      character(*), parameter :: thisnam = 'ematqkgir'
!      integer :: ikq, ig, ig1, ig2, ig3, igk0, igk, iv(3), iv1(3), iv3(3), ivu(3)
!      integer, allocatable :: aigk0(:), aigk(:)
!
!      ! Grid index for k+q point
!      ikq = ikmapikq(ik, iq)
!
!      allocate(aigk0(ngkmax0), aigk(ngkmax))
!
!      ! Positive umklapp G-vector
!      ivu(:) = nint(vkl0(:, ik)+vql(:, iq)-vkl(:, ikq))
!
!      ! Precalculate for speedup
!      aigk0(:) = igkig0(:, 1, ik)
!      aigk(:) = igkig(:, 1, ikq)
!      ig3 = igqig(igq, iq)
!      iv3(:) = ivg(:, ig3)
!
!      do igk0 = 1, ngk0(1, ik)
!        ig1 = aigk0(igk0)
!        iv1(:) = ivg(:, ig1) + iv3(:)
!        do igk = 1, ngk(1, ikq)
!          ig2 = aigk(igk)
!          ! Umklapp of k+q vector included
!          iv(:) = iv1(:) - (ivg(:, ig2)-ivu(:))
!          ig = ivgig(iv(1), iv(2), iv(3))
!          xihir(igk0, igk) = cfunig(ig)
!        end do
!      end do
!
!      deallocate(aigk0, aigk)
!    end subroutine ematqkgir_hack
!
    !BOP
    ! !ROUTINE: ematrad
    ! !INTERFACE:
    subroutine ematrad_hack(iq)
    ! !USES:
      use mod_misc, only: filext
      use mod_muffin_tin, only: nrmtmax, nrmt
      use modinput, only: input
      use m_getunit
    ! !DESCRIPTION:
    ! This routine is used in the construction of the plane wave matrix elements.
    ! It calculates the involved radial integrals for the APW-APW, APW-LO and LO-LO
    ! combinations. (Compare the R quantities section 8.1 of Sagmeisters thesis.)
    !
    ! !REVISION HISTORY:
    ! Added to documentation scheme. (Aurich)
    !EOP
    !BOC

      implicit none

      ! Arguments
      integer, intent(in) :: iq

      ! Local variables
      integer :: is, ia, ias, nr, ir, igq
      integer :: l1, l2, l3,lio,liomax
      integer :: ilo, ilo1, ilo2, io, io1, io2
      real(8) :: t1
      integer :: lmax1, lmax2, lmax3
      integer :: u11, u22, u33
      ! Automatic arrays
      real(8) :: r2(nrmtmax), fr(nrmtmax), gr(nrmtmax), cf(3, nrmtmax)
      ! Allocatable arrays
      real(8), allocatable :: jl(:, :), jhelp(:)
      integer, allocatable :: lio2l(:), lio2io(:)

      ! Maximal l in basis
      lmax1 = max(lmaxapw, lolmax)
      ! Rayleigh expansions cutoff
      lmax2 = lrayleigh

      ! lmax1 and lmax3 should be the same!
      lmax3 = lmax1

      ! Allocate arrays for radial integrals and Bessel functions
      if(allocated(riaa)) deallocate(riaa)
      if(allocated(riloa)) deallocate(riloa)
      if(allocated(rilolo)) deallocate(rilolo)
      allocate(riaa(0:lmax1, apwordmax, 0:lmax3, apwordmax, 0:lmax2, natmtot, Gqgrid%ngk(1, iq)))
      allocate(riloa(nlomax, 0:lmax3, apwordmax, 0:lmax2, natmtot, Gqgrid%ngk(1, iq)))
      allocate(rilolo(nlomax, nlomax, 0:lmax2, natmtot, Gqgrid%ngk(1, iq)))

      ! Allocate temporary arrays
      allocate(jl(nrmtmax,0:lmax2))
      allocate(jhelp(0:lmax2))
      allocate(lio2l((lmax1+1)*apwordmax))
      allocate(lio2io((lmax1+1)*apwordmax))

      jl(:, :) = 0.d0
      jhelp(:) = 0.d0

      ! Zero arrays for radial integrals
      riaa(:, :, :, :, :, :, :) = 0.d0
      riloa(:, :, :, :, :, :) = 0.d0
      rilolo(:, :, :, :, :) = 0.d0

      ! Begin loop over G+q vectors
      do igq = 1, Gqgrid%ngk(1, iq)

        ! Begin loop over species
        do is = 1, nspecies

          nr = nrmt(is)
          do ir = 1, nr
            ! Calculate r^2
            r2(ir) = spr(ir, is) ** 2
            ! Calculate spherical Bessel functions of first kind j_l(|g+q|r_a)
            call sbessel(lmax2, Gqgrid%gkc(igq, 1, iq)*spr(ir, is), jhelp)
            jl(ir,:) = jhelp(:)
          end do

          lio=0
          do l1 = 0, lmax1
            do io1 = 1, apword(l1, is)
              lio=lio+1
              lio2l(lio)=l1
              lio2io(lio)=io1
            end do
          end do
          liomax=lio

          ! Begin loop over atoms
          do ia = 1, natoms(is)
            ias = idxas(ia, is)
            !----------------!
            !     apw-apw    !
            !----------------!
#ifdef USEOMP
    !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(lio,l1,io1,l3,io2,ir,t1,fr,cf,gr)
    !$OMP DO
#endif
            do lio=1,liomax
              l1=lio2l(lio)
              io1=lio2io(lio)
              do l3 = 0, lmax3
                do io2 = 1, apword(l3, is)
                  do l2 = 0, lmax2
                    do ir = 1, nr
                      t1 = apwfr(ir, 1, io1, l1, ias) * apwfr(ir, 1, io2, l3, ias) * r2(ir)
                      fr(ir) = t1 * jl(ir,l2)
                    end do
                    call fderiv(-1, nr, spr(:, is), fr, gr, cf)
                    riaa(l3, io2, l1, io1, l2, ias, igq) = gr(nr)
                  end do ! l2
                end do ! io2
              end do ! l3
            end do ! l1
#ifdef USEOMP
    !$OMP END DO
    !$OMP END PARALLEL
#endif
            !----------------------------!
            !     local-orbital-apw      !
            !----------------------------!
            do ilo = 1, nlorb(is)
              l1 = lorbl(ilo, is)
              do l3 = 0, lmax3
                do io = 1, apword(l3, is)
                  do l2 = 0, lmax2
                    do ir = 1, nr
                      t1 = lofr(ir, 1, ilo, ias) * apwfr(ir, 1, io, l3, ias) * r2(ir)
                      fr(ir) = t1 * jl(ir,l2)
                    end do
                    call fderiv(-1, nr, spr(:, is), fr, gr, cf)
                    riloa(ilo, l3, io, l2, ias, igq) = gr(nr)
                  end do ! l2
                end do ! io
              end do ! l3
            end do ! ilo

            !------------------------------------!
            !     local-orbital-local-orbital    !
            !------------------------------------!
            do ilo1 = 1, nlorb(is)
              l1 = lorbl(ilo1, is)
              do ilo2 = 1, nlorb(is)
                l3 = lorbl(ilo2, is)
                do l2 = 0, lmax2
                  do ir = 1, nr
                    t1 = lofr(ir, 1, ilo1, ias) * lofr(ir, 1, ilo2, ias) * r2(ir)
                    fr(ir) = t1 * jl(ir,l2)
                  end do
                  call fderiv(-1, nr, spr(:, is), fr, gr, cf)
                  rilolo(ilo1, ilo2, l2, ias, igq) = gr(nr)
                end do ! l2
              end do ! ilo2
            end do ! ilo1

          ! End loops over atoms and species
          end do
        end do

      ! End loop over G+q vectors
      end do

      ! Deallocate
      deallocate(jl, jhelp)
      deallocate(lio2l,lio2io)

    end subroutine ematrad_hack

!    !EOC
!
!    subroutine ematqalloc_hack
!      use mod_eigensystem, only: nmatmax
!      use mod_eigenvalue_occupancy, only: nstfv, nstsv
!      use mod_spin, only: nspnfv
!      use mod_APW_LO, only: apwordmax
!      use mod_atoms, only: natmtot
!      use modxs, only: evecfv, evecfv0, evalsv0, nmatmax0,&
!                     & apwcmt, lmmaxapwwf, apwcmt0,&
!
!      implicit none
!
!      ! Allocate eigenvalue and eigenvector arrays
!      if(allocated(evecfv)) deallocate(evecfv)
!      allocate(evecfv(nmatmax, nstfv, nspnfv))
!      if(allocated(evecfv0)) deallocate(evecfv0)
!      allocate(evecfv0(nmatmax0, nstfv, nspnfv))
!      
!      ! Allocate contracted coefficients array for apw-part
!      if(allocated(apwcmt)) deallocate(apwcmt)
!      allocate(apwcmt(nstfv, apwordmax, lmmaxapwwf, natmtot))
!      if(allocated(apwcmt0)) deallocate(apwcmt0)
!      allocate(apwcmt0(nstfv, apwordmax, lmmaxapwwf, natmtot))
!
!    end subroutine ematqalloc_hack

end module m_ematqk_hack
