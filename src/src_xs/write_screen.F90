! Copyright (C) 2008-2010 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module mod_write_screen

  contains
  
  !BOP
  ! !ROUTINE: scrcoulint
  ! !INTERFACE:
  subroutine write_screen 
    ! !USES:
      use mod_misc, only: filext
      use modinput, only: input
      use modmpi
      use constants, only: zzero, zone, fourpi
      use mod_APW_LO, only: lolmax
      use mod_qpoint, only: iqmap, vql, vqc, nqpt, ivq, wqpt
      use mod_lattice, only: omega
      use mod_symmetry, only: maxsymcrys
      use modxs, only: xsgnt, unitout,&
                    & ngqmax,&
                    & nqptr, qpari, qparf, ivqr,&
                    & ngq, ppari, pparf, iqmapr,&
                    & vqlr, vqcr, wqptr, ngqr,&
                    & bcbs, ematraddir, eps0dirname,&
                    & filext0, usefilext0, iqmt0, iqmt1,&
                    & iqmtgamma,&
                    & bsedl, bsedu, bsedd, bsed
      use m_xsgauntgen
      use m_findgntn0
      use m_writevars
      use m_genfilname
      use m_getunit
      use m_ematqk
      use m_putgetbsemat
      use modbse
      use mod_xsgrids
      use mod_Gkvector, only: gkmax
      use xhdf5, only: xhdf5_type
      use mod_hdf5, only: fhdf5, hdf5_exist_group, hdf5_create_group, hdf5_write 
      use os_utils, only: join_paths
    ! !DESCRIPTION:
    !   Calculates the screened Coulomb interaction $W_{\mathbf{GG}'}(\mathbf{q})$
    !   on the full $\mathbf{q}$-grid and writes it to hdf5 file
    !
    ! !REVISION HISTORY:
    !   Extracted from scroulint.F90 (Vorwerk)
    !   Forked from scrcoulint.F90 and adapted for non-TDA BSE and finite Q. (Aurich)
    !EOP
    !BOC      
    
      implicit none
    
    
      ! Local variables
      character(*), parameter :: thisname = 'scrcoulint'
    
      ! ik,jk block of W matrix (final product)
      complex(8), allocatable :: sccli(:,:)
    
      ! Diagonal of W (needed for MB TDDFT kernels)
      complex(8), allocatable :: bsedt(:, :)
    
      ! Auxilliary arrays for the construction of sccli
      complex(8), allocatable :: sccli_t1(:, :), sccli_t2(:, :, :, :)
      complex(8), allocatable :: zm(:,:)
      ! Plane wave arrays
      complex(8), allocatable :: muu(:, :, :), cmuu(:, :)
      complex(8), allocatable :: moo(:, :, :), cmoo(:, :)
      complex(8), allocatable :: mou(:, :, :), cmou(:, :)
      complex(8), allocatable :: muo(:, :, :), cmuo(:, :)
      ! The fourier coefficients of the screend coulomb potential
      complex(8), allocatable :: wfc(:, :)
      ! Auxilliary arrays/vars for the construction of wfc
      complex(8), allocatable :: scieffg(:, :, :)
      complex(8), allocatable :: phf(:, :)
      ! Symmerty maps creation 
      integer(4) :: sc(maxsymcrys), ivgsc(3, maxsymcrys)
      integer(4), allocatable :: igqmap(:)
      integer(4) :: jsym, jsymi
      integer(4) :: nsc, ivgsym(3)
      logical :: tphf
      ! Mappings of jk ik combinations to q points
      integer(4) :: ikkp, iknr, jknr, ikpnr, ikmnr, jkpnr, jkmnr, ik, jk
      integer(4) :: iqrnr, iqr, iq
      real(8) :: vqr(3), vq(3)
      integer(4) :: numgq
      logical :: tq0
      ! Number of occupied/unoccupied states at ik and jk
      integer(4) :: ino, inu, jno, jnu
      ! Number of transitions at ik and jk
      integer(4) :: inou, jnou
      ! Number of (o/u)_i (o/u)_j combinations
      integer(4) :: noo, nuu, nou, nuo
      ! State loop indices
      integer(4) :: io, jo, iu, ju
      ! Combinded loop indices 
      integer(4) :: jaoff, iaoff, ia, ja
      ! Maximal l used in the APWs and LOs
      !   Influences quality of plane wave matrix elements
      integer(4) :: maxl_apwlo
      ! Maximal l used in the Reghley expansion of exponential
      !   Influences quality of plane wave matrix elements
      integer(4) :: maxl_e
      ! Maximal l used in the APWs and LOs in the groundstate calculations
      !   Influences quality of eigencoefficients.
      integer(4) :: maxl_mat
      ! Aux.
      integer(4) :: j1, j2
      complex(8) :: pref, zt1
      real(8) :: t1
      ! Timing vars
      real(8) :: tscc1, tscc0
    
      ! Auxilliary strings
      character(256) :: syscommand, fileext_scr_read, fileext_ematrad_write
      ! HDF5 variables
      character(:), allocatable :: gname, group
      character(4) :: ciq
      ! External functions
      logical, external :: tqgamma
    
      real(8) :: vqoff(3)
      real(8), parameter :: epslat = 1.0d-8
    
      logical :: fsameq, fsamekp, fsamekm

      type(xhdf5_type) :: h5 
    
      !---------------!
      !   main part   !
      !---------------!
    
      !write(*,*) "Hello, this is scrcoulint at rank:", rank
    
      ! General setup
      call init0()
      ! k-point setup
      call init1()
      ! Save variables of the unshifinit0
      ! q-point and qmt-point setup
      !   Init 2 sets up (task 440):
      !   * A list of momentum transfer vectors form the q-point list (modxs::vqmtl)
      !   * The reduced unshifted q-grid (modxs::vqlr etc)
      !   * The non-reduced unshifted q-grid (mod_qpoint::vql etc)
      !   * Offset of the k+q grid derived from k offset an q point (modxs::qvkloff)
      !     which is just equal to the k-offset, since it is the unshifted q grid
      !   * non-reduced mapping between ik,q and ik' grids (modxs::ikmapikq)
      !   * G+q quantities for the non-reduced unshifted q-grid (modxs)
      !   * The square root of the Coulomb potential for the non-reduced q points
      !   * Reads STATE.OUT
      !   * Generates radial functions (mod_APW_LO)
      call init2
    
      ! Making folder for the radial integals pertaining to the plane wave matrix elements
      ematraddir = 'EMATRAD'
    
      syscommand = 'test ! -e '//trim(adjustl(ematraddir))&
        & //' && mkdir '//trim(adjustl(ematraddir))//' &> /dev/null'
      call system(trim(adjustl(syscommand)))
    
      ! Generate gaunt coefficients used in the construction of 
      ! the plane wave matrix elements in ematqk.
    
      ! gs%lmaxapw defaults to 8, BUT xs%lmaxapw defaults to 10 
      ! in xsinit gs%lmaxapw is overwritten with the xs%lmaxapw value.
      ! --> xs%lmaxapw = 10
      maxl_apwlo = max(input%groundstate%lmaxapw, lolmax)
    
      ! xs%lmaxapwwf defaults to -1, in which case it is set to gs%lmaxmat by init2
      ! which in trun was previously overwritten by xs%lmaxmat by xsinit. 
      ! xs%lmatmax defaults to 5.
      ! --> xs%lmaxapwwf = xs%lmaxmat = 5
      maxl_mat = max(input%xs%lmaxapwwf, lolmax)
    
      ! xs%lmaxemat defaults to 3
      maxl_e = input%xs%lmaxemat
    
      ! Allocates xsgnt array of shape 
      ! ((maxl_apwlo+1)**2, (maxl_e+1)**2, (maxl_apwlo+1)**2)
      ! and fills it with the corresponding gaunt coefficients
      ! Note: In the generation of the gaunt coefficients l1 and l3 correspond
      !       to the APW/LOs while l2 corresponds to the exponetial.
      call xsgauntgen(maxl_apwlo, maxl_e, maxl_apwlo)
      ! Find indices for non-zero gaunt coefficients in xsgnt,
      ! and creates index maps, e.g. given l1,m1,l2,m2 -> non zero l3,m3
      ! Up l1 up to maxl_mat, l2 up to maxl_mat, l3 up to maxl_e
      ! Note: In the maps and following plane wave matrix elements l1 and l2 correspond
      !       to the APW/LOs while l3 corresponds to the exponetial.
      call findgntn0(maxl_mat, maxl_mat, maxl_e, xsgnt)
    
      if (input%xs%BSE%outputlevelnumber == 1) then
        write(unitout, '(a)') 'Info(' // thisname // '):&
          & Gaunt coefficients generated within lmax values:'
        write(unitout, '(a, i8)') "lmax1 = lmaxapw =", input%groundstate%lmaxapw
        write(unitout, '(a, i8)') "lmax2 = lmaxemat=", input%xs%lmaxemat
        write(unitout, '(a, i8)') "lmax3 = lmaxapw =", input%groundstate%lmaxapw
        call flushifc(unitout)
      end if
    
      ! Read Fermi energy from file EFERMI
      ! Use EFERMI_QMT001.OUT (corresponding to the xs groundstate run for the unshifted k grid)
      call genfilname(iqmt=iqmtgamma, setfilext=.true.)
      call readfermi
    
      ! Set ist* variables and ksgap in modxs using findocclims
      ! This also reads in 
      ! mod_eigevalue_occupancy:evalsv, mod_eigevalue_occupancy:occsv 
      ! modxs:evalsv0, modxs:occsv0
      call printline(unitout, '-')
      write(unitout, '(a)') 'Info(' // thisname // '):&
        & Inspecting occupations...'
      call flushifc(unitout)
    
      !call setranges_modxs(iqmt)
    
      ! Select relevant transitions for the construction
      ! of the BSE hamiltonian
      ! Also sets nkkp_bse, nk_bse 
      call printline(unitout, '-')
      write(unitout, '(a)') 'Info(' // thisname // '):&
        & Selecting transitions...'
      call flushifc(unitout)
    
      !call select_transitions(iqmt, serial=.false.)
    
    
      call printline(unitout, '-')
      write(unitout, '("Info(",a,"): Size of file ",a," will be about ", f12.6, " GB" )')&
        & trim(thisname), trim(scclifbasename),&
        & int(nou_bse_max,8)**2*int(nkkp_bse,8)*16.0d0/1024.0d0**3
      call flushifc(unitout)
    
      !------------------------------------!
      ! GENERATE FOURIER COEFFICIENTS OF W !     
      ! (and radial integrals for emat)    !
      !------------------------------------!
    
      call xsgrids_init(totalqlmt(1:3, iqmtgamma), gkmax)
    
      !--------------------------------------------------------------------------------!
      ! Setup q points
      !--------------------------------------------------------------------------------!
      ! Get q grid offset
      ! This is zero
      vqoff = q%qset%vkloff
      fsameq=.true.
    
      ! Make reduced q-grid with possible offset
      ! This sets up also G+q quantities and the square root of the Coulomb potential
      ! but the second call below will override these
      call init1offs(k_kqmtp%kset%vkloff)
      vql = 0.0d0
      call init2offs(vqoff, input%xs%reduceq)
      ! Copy results q-ponit results form mod_qpoint into modxs variables
      nqptr = nqpt
      ivqr = ivq
      iqmapr = iqmap
      vqlr = vql
      vqcr = vqc
      wqptr = wqpt
      ! Make non-reduced q-grid with possible offset (written into mod_qpoint variables)
      ! This sets up also G+q quantities and the square root of the Coulomb potential
      call init2offs(vqoff, .false.)
      write(unitout, '(a, i6)') 'Info(' // thisname // '):&
        & Number of reduced q-points: ', nqptr
      write(unitout, '(a, i6)') 'Info(' // thisname // '):&
        & Number of non-reduced q-points: ', nqpt
      write(unitout,*)
      call flushifc(unitout)
      ! Make also ngqr
      if(allocated(ngqr)) deallocate(ngqr)
      allocate(ngqr(nqptr))
      ngqr=0
      do iq= 1, nqpt
        iqr = iqmapr(ivq(1,iq), ivq(2,iq), ivq(3,iq))
        ngqr(iqr) = ngq(iq)
      end do
      !--------------------------------------------------------------------------------!
    
      !--------------------------------------------------------------------------------!
      ! Make W(G,G',q) Fourier coefficients and radial integrals for the PW elements
      !--------------------------------------------------------------------------------!
      ! Allocate local arrays for screened coulomb interaction and
      ! W(G,G',qr)
      allocate(scieffg(ngqmax, ngqmax, nqptr))
      scieffg(:, :, :) = zzero
      ! Phases for transformations form reduced q points to non-reduced ones.
      allocate(phf(ngqmax, ngqmax))
    
      ! Parallelize over reduced q-point set
      call genparidxran('q', nqptr)
    
      ! Set file to read the static, non-broadened RPA screening from 
      eps0dirname = 'EPS0'
        ! q-grid (zero offset)
      call genfilname(fileext=fileext_scr_read)
      call genfilname(iqmt=iqmtgamma, fileext=fileext_ematrad_write)
    
      write(unitout, '("Info(scrcoulint):&
        & Calculating W(G1,G2,qr) fourier coefficients")')
      call timesec(tscc0)
    
      do iqr = qpari, qparf ! Reduced q
    
        ! Locate reduced q-point in non-reduced set
        iqrnr = iqmap(ivqr(1,iqr), ivqr(2,iqr), ivqr(3,iqr))
    
        ! Get number of G+q vectors for current q
        numgq = ngq(iqrnr)
    
        ! Calculate effective screened coulomb interaction
        ! by inverting the symmetrized IP dielectric matrix for a given q and
        ! 0 frequency and then multiplying
        ! it with v^{1/2} from both sides.
        filext = fileext_scr_read
        call genscclieff(iqr, iqrnr, ngqmax, numgq, scieffg(:,:,iqr))
    
        ! Generate radial integrals for matrix elements of plane wave
        ! and save them to disk.
        filext = fileext_ematrad_write
        call putematrad(iqr, iqrnr)
      end do
    
    
      ! Set file extesion for later read EMATRAD in getematrad
      ! (some ranks may not participate in the qr loop above)
      filext = fileext_ematrad_write
    
      ! Communicate array-parts wrt. q-points
      call mpi_allgatherv_ifc(set=nqptr, rlen=ngqmax*ngqmax,&
        & zbuf=scieffg, inplace=.true., comm=mpiglobal)
      ! write W(G,G,q) to file if necessary


      call barrier(callername=trim(thisname))
    
      if (input%xs%BSE%writepotential.and. mpiglobal%rank == 0) then
        call h5%initialize(fhdf5, mpiglobal%comm)
        call h5%initialize_group('/', 'screenedpotential')
        gname = 'screenedpotential'

        ! loop over all reduced q-vectors
        do iqr=1, nqptr
          write(ciq, '(I4.4)') iqr
          call h5%initialize_group(gname, ciq)
          group = join_paths(gname, ciq)
          call h5%write(group, "wqq", scieffg(:, :, iqr), [1, 1], shape(scieffg(:, :, iqr)))
        end do
      end if

    end subroutine write_screen 
    !EOC
    
end module mod_write_screen   