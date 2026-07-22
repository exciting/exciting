!BOP
!
!!ROUTINE: qdepwtet
!
!!INTERFACE:
!
subroutine qdepwtet(iq, iomstart, iomend, ndim)
!
!!DESCRIPTION:
!
! This subroutine calculates the weights for q dependent BZ integration
! using LIBBZINT
!
!
!!USES:
    use constants, only: zzero
    use mod_atoms, only: idxas
    use mod_bands, only: nomax, numin, evalfv, nstdf
    use mod_corestate, only: evalcr
    use mod_eigenvalue_occupancy, only: efermi
    use mod_lattice, only: bvec, binv
    use mod_symmetry, only: symlat, lsplsymc, nsymcrys, symlatc, & 
                            find_equivalent_wavevectors, get_equivalent_qpairs
    use modinput, only: input
    use modgw,   only : fnm_tet, freq, ncmax, kset, kqset, &
    &                   ncg, corind, fdebug, time_bzinit
        use precision, only: i32, dp

    implicit none

    integer(i32), intent(in) :: iq
    integer(i32), intent(in) :: iomstart, iomend
    integer(i32), intent(in) :: ndim
    

!!LOCAL VARIABLES:
    integer(i32) :: ik, ikp, ib, ic, icg
    integer(i32) :: ia, is, ias, ie1, ie2
    integer(i32) :: iom
    integer(i32) :: fflg, sgw
    integer(i32) :: mini, mend

    real(dp) :: emaxb ! maximum energy of the second band
    real(dp) :: edif, edsq, omsq
    real(dp) :: sfact

    real(dp), allocatable :: eval(:,:), eval_pair(:,:)
    real(dp), allocatable :: cwpar(:,:,:)
    real(dp), allocatable :: cwparsurf(:,:,:)

    logical :: is_core, is_realfreq

    real(dp) :: tstart, tend

    ! Equivalent q-points fnm 
    complex(dp), allocatable :: fnm_equiv(:,:)

    ! Symmetry related indexes
    integer(i32), allocatable :: iqeq_list(:), isymeq_list(:)

    ! Index for iq equivalence
    integer(i32) :: isym, iqeq, nsym, neq, neqpairs
    integer(i32), allocatable :: point_pairs(:,:)

!EOP
!BOC
    call timesec(tstart)
    mini = lbound( fnm_tet, 2 )
    mend = ubound( fnm_tet, 2 )

    !---------------------------------------------------------------------
    ! Initialization
    !---------------------------------------------------------------------
    nsym = merge(nsymcrys, 1_i32, input%gw%enforceCrystalSymmetryTetrahedron)
    allocate(fnm_equiv(kqset%nkpt,nsym), source=zzero)

    ! Find equivalent q-points
    ! Notice that a point can be found more than once
    call find_equivalent_wavevectors( 3, kqset%vql(:, iq), kqset%vql(:,:), kqset%nkpt, &
      symlat(:, :, lsplsymc(1:nsym)), nsym, iqeq_list, isymeq_list)
    neq = size(iqeq_list)

    ! real or imaginary frequencies
    select case (freq%fconv)
      case('nofreq')
        fflg = 1
      case('refreq')
        fflg = 2
      case('imfreq')
        fflg = 3
    end select
    sgw = 5-2*fflg

    ! Are we dealing with real frequencies
    is_realfreq = (fflg == 2)

    !====================
    ! valence-valence
    !====================
    allocate(eval(nstdf,kqset%nkpt))
    do ik = 1, kqset%nkpt
      ikp = kset%ik2ikp(ik)
      eval(1:nstdf,ik) = evalfv(1:nstdf,ikp)
    end do

    allocate(cwpar(2,2,kqset%nkpt))
    if (is_realfreq) allocate(cwparsurf(2,2,kqset%nkpt))

    allocate(eval_pair(2,kqset%nkpt))
    allocate(point_pairs(2,nsym))

    !$omp parallel do collapse(3) default(none) schedule(dynamic) &
    !$omp shared(iq, iomstart, iomend, mini, mend, ndim, nomax, corind, idxas, eval, evalcr) &
    !$omp shared(efermi, sgw, freq, neq, iqeq_list, kqset, fflg, is_realfreq, fnm_tet, nsym) &
    !$omp shared(lsplsymc, symlatc, bvec, binv, input) &
    !$omp private(iom, ie2, ie1, is_core, icg, is, ia, ic, ias, emaxb, omsq, isym, iqeq, ik, edif, edsq, neqpairs) &
    !$omp firstprivate(eval_pair, point_pairs, cwpar, cwparsurf, fnm_equiv) 
    do iom = iomstart, iomend
      do ie2 = mini, mend
        do ie1 = 1, ndim

          eval_pair(2,:) = eval(ie2,:)

          is_core = ( ie1 > nomax )

          ! Valence band
          if (.not. is_core) then
            eval_pair(1,:) = eval(ie1,:)
          else
            icg = ie1 - nomax
            is  = corind(icg,1)
            ia  = corind(icg,2)
            ic  = corind(icg,6)
            ias = idxas(ia,is)
            eval_pair(1,1:kqset%nkpt) = evalcr(ic,ias)
            emaxb = maxval(eval_pair(2,:))
            if (emaxb <= efermi) cycle 
            omsq = sgw*freq%freqs(iom)*freq%freqs(iom)
          end if
          
          do isym = 1, neq
            
            iqeq = iqeq_list(isym)

            call tetcw(kqset%nkpt, kqset%ntet, 2, kqset%wtet, eval_pair, &
                       kqset%tnodes, kqset%linkq(:,iqeq), kqset%kqid(:,iqeq), &
                       kqset%tvol, efermi, freq%freqs(iom), merge(1, fflg, is_core), &
                       cwpar)

            if (is_realfreq) call tetcw(kqset%nkpt, kqset%ntet, 2, kqset%wtet, eval_pair, &
                                    kqset%tnodes, kqset%linkq(:,iqeq), kqset%kqid(:,iqeq), &
                                    kqset%tvol, efermi, freq%freqs(iom), 4, cwparsurf)

            if (is_core) then
              do ik = 1, kqset%nkpt
                edif = eval(ie2,ik)-evalcr(ic,ias)
                edsq = edif*edif
                cwpar(1,2,ik) = 2.0_dp*cwpar(1,2,ik)*edif/(omsq-edsq) ! <-- why 2?
              end do
            else 
              do ik = 1, kqset%nkpt
                if (eval(ie2,kqset%kqid(ik,iqeq)) > 900.0_dp) then
                  cwpar(1,2,ik) = 0.0_dp
                  if (is_realfreq) cwparsurf(1,2,ik) = 0.0_dp
                end if
              end do
            end if

            if (is_realfreq) then
              fnm_equiv(1:kqset%nkpt,isym) = cmplx(cwpar(1,2,1:kqset%nkpt),cwparsurf(1,2,1:kqset%nkpt))
            else 
              fnm_equiv(1:kqset%nkpt,isym) = cmplx(cwpar(1,2,1:kqset%nkpt),0.0)
            end if

          end do ! isym

          ! The symmetry equivalent pairs are not equivalent
          ! Thus, we average over symmetry equivalent pairs with their appropiate weight
          ! While this approach cannot be physically justified, and indeed a much better
          ! solution will be to modify LIBBZINT to be compliant with the crystal symmetry,
          ! averaging is a commonly used solution in other fields as transport.
          if (nsym /= 1) then
            do ik = 1, kqset%nkpt
              ! Computing the symmetry equivalent pairs (q,k), notice that the 
              ! call will provide all even repeated pairs, that is with their appropriate weight
              ! within the equivalence group.  
              call get_equivalent_qpairs(iq, ik, kqset%vkl(:,:), point_pairs, nsym, lsplsymc, symlatc, bvec, binv, &
                                         input%gw%ngridq, input%gw%vqloff)
              ! Few choices for the mesh can break crystal symmetry; most of times this is unwanted
              ! Nevertheless, it has few useful cases, so we check for this case.
              ! For the weight, we only consider
              ! operations mapping point-pairs within the mesh
              neqpairs = nsym - count(point_pairs(1,:) == -1)
              do isym = 1, nsym
                if (point_pairs(1,isym) == -1) cycle ! If the symmetry operation do not map the points to valid ones, cycle
                fnm_tet(ie1,ie2,iom,ik) = fnm_tet(ie1,ie2,iom,ik) + fnm_equiv(point_pairs(2,isym),isym)
              end do
              ! Apply the approapiate weight to the summation over the equivalent pairs
              fnm_tet(ie1,ie2,iom,ik) = fnm_tet(ie1,ie2,iom,ik) / neqpairs
            end do
          else
            fnm_tet(ie1,ie2,iom,:) = fnm_equiv(:,1)
          end if

        end do !ie1
      end do !ie2
    end do !iom
    !$omp end parallel do

    ! spin degeneracy: I'm not sure about this prefactor for core states
    sfact = 2.0_dp
    fnm_tet = sfact*fnm_tet

    !-------------------------
    ! Debugging info
    !-------------------------
    if (input%gw%debug) then
      write(fdebug,*)'------------------------------------------------------'
      write(fdebug,*)'       convolution weights for iq =',iq
      write(fdebug,*)'------------------------------------------------------'
      write(fdebug,*)
      iom = 1
      do ik = 1, kqset%nkpt
      do ic = 1, ndim
        do ib = mini, mend
          write(fdebug,1) ic, ib, iom, ik, fnm_tet(ic,ib,iom,ik)
        end do
        end do
      end do
      1 format('  ic =',i4,' ib =',i4,' iom =',i4,' ik =',i4,' fnm =',2g16.8)
    end if ! debug

    ! timing info
    call timesec(tend)
    time_bzinit =  time_bzinit+tend-tstart

    return

end subroutine
!EOC
