
subroutine task_dos()

    use modinput
    use modmain
    use modgw
    implicit none

    integer(4) :: ik, ist, nstsv0
    integer(4) :: iw, n, nsk(3)
    integer(4) :: nsmdos, nwdos, ngrdos
    real(8) :: winddos(2)
    real(8) :: dw, egap
    real(8), allocatable :: e(:,:)
    real(8), allocatable :: f(:,:)
    real(8), allocatable :: w(:)
    real(8), allocatable :: g(:)
    real(8), allocatable :: bc(:,:,:,:,:)

    integer :: is, ia, ias
    integer :: lmax, lmmax, l, m, lm
    integer :: nstqp
    character(80) :: fname

    !-----------------
    ! Initialization
    !-----------------

    ! Note: DOS calculations are performed employing Groundstate Options
    ! - ngridk
    ! - smearing etc.

    call rereadinput() ! overwrite some GW parameters
    call init0()
    call init1()
    ! write(*,*) 'nkpt=', nkpt
    ! write(*,*) input%groundstate%ngridk

    ! read KS data
    do ik = 1, nkpt
      call getevalsv(vkl(:,ik), evalsv(:,ik))
    end do

    ! shift KS energies
    call readfermi

    ! read QP energies from file and perform Fourier interpolation (if required)
    call getevalqp(nkpt,vkl,evalsv)
    
    ! GW number of states
    nbgw = min(nstsv,nbgw)
    nstqp = nbgw-ibgw+1
    allocate(e(ibgw:nbgw,nkpt))
    e(ibgw:nbgw,:) = evalsv(ibgw:nbgw,:)

    ! diagonal of spin density matrix for weight
    allocate(f(ibgw:nbgw,nkpt))
    f(:,:) = 1.d0

    !-----------------------------      
    ! DOS parameters
    !-----------------------------
    if (.not.associated(input%properties)) &
    &  input%properties => getstructproperties(emptynode)
    if (.not.associated(input%properties%dos)) &
    &  input%properties%dos => getstructdos(emptynode)

    nsmdos = input%properties%dos%nsmdos
    nwdos = input%properties%dos%nwdos
    ngrdos = input%properties%dos%ngrdos
    winddos(:) = input%properties%dos%winddos(:)    
    
    ! generate energy grid
    allocate(w(nwdos))
    dw = (winddos(2)-winddos(1))/dble(nwdos)
    do iw = 1, nwdos
      w(iw) = dw*dble(iw-1)+winddos(1)
    end do
      
    ! number of subdivisions used for interpolation
    nsk(:) = max(ngrdos/input%groundstate%ngridk(:),1)
    
    ! BZ integration
    allocate(g(nwdos)) ! DOS
    call brzint(nsmdos, input%groundstate%ngridk, nsk, ikmap, &
    &           nwdos, winddos, nstqp, nstqp, e, f, g)
    g(:) = occmax*g(:)

    ! output file
    open(50,file='TDOS-QP.OUT',action='WRITE',form='FORMATTED')
    do iw = 1, nwdos
      write(50,'(2G18.10)') w(iw), g(iw)
    end do
    close(50)

    if (input%properties%dos%lmirep) then
      lmax = 4
      lmmax = (lmax+1)**2
      allocate(bc(lmmax,nspinor,natmtot,nstsv,nkpt))
      call calc_band_character(lmax,lmmax,bc)
      do is = 1, nspecies
      do ia = 1, natoms(is)
        ias = idxas(ia,is)
        write(fname, '("PDOS-QP_S", I2.2, "_A", I4.4, ".OUT")') is, ia
        open(50, File=trim(fname), Action='WRITE', Form='FORMATTED')
        do l = 0, lmax
        do m = - l, l
          lm = idxlm(l,m)
          do ik = 1, nkpt
            do ist = ibgw, nbgw
              f(ist,ik) = bc(lm,1,ias,ist,ik)
            end do
          end do
          call brzint(nsmdos, input%groundstate%ngridk, nsk, ikmap, &
          &           nwdos, winddos, nstqp, nstqp, e, f, g)
          do iw = 1, nwdos
            write(50, '(2G18.10)') w(iw), occmax*g(iw)
          end do
          write(50,*)
        end do
        end do
        close(50)
      end do
      end do
      deallocate(bc)

    end if
    
    deallocate(e,w,f,g)
 
    return

contains

    !-----------------------------------------------------------------
    subroutine calc_band_character(lmax,lmmax,bc)
      use modmain
      implicit none
      integer, intent(in)  :: lmax
      integer, intent(in)  :: lmmax
      real(8), intent(out) :: bc(lmmax,nspinor,natmtot,nstsv,nkpt)
      ! local
      integer :: ik, ispn, jspn, is, ia, ist, ias, lm
      real(8),    allocatable :: elm(:,:)
      complex(8), allocatable :: ulm(:,:,:)
      complex(8), allocatable :: a(:,:)
      complex(8), allocatable :: dmat(:,:,:,:,:)
      complex(8), allocatable :: apwalm(:,:,:,:,:)
      complex(8), allocatable :: evecfv(:,:,:)
      complex(8), allocatable :: evecsv(:,:)

      allocate(elm(lmmax,natmtot))
      allocate(ulm(lmmax,lmmax,natmtot))
      allocate(a(lmmax,lmmax))

      allocate(dmat(lmmax,lmmax,nspinor,nspinor,nstsv))
      allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
      allocate(evecfv(nmatmax,nstfv,nspnfv))
      allocate(evecsv(nstsv,nstsv))

      ! read density and potentials from file
      call readstate
      ! read Fermi energy from file
      call readfermi
      ! find the new linearisation energies
      call linengy
      ! generate the APW radial functions
      call genapwfr
      ! generate the local-orbital radial functions
      call genlofr
      
      ! generate unitary matrices which convert the (l,m) basis into the irreducible
      ! representation basis of the symmetry group at each atomic site
      call genlmirep(lmax,lmmax,elm,ulm)
      
      ! loop over k-points
      do ik = 1, nkpt

        call getevecfv(vkl(1,ik), vgkl(:,:,:,ik), evecfv)
        call getevecsv(vkl(1,ik), evecsv)
        
        ! find the matching coefficients
        do ispn = 1, nspnfv
          call match(ngk(ispn,ik), gkc(:,ispn,ik), &
          &          tpgkc(:,:,ispn,ik), sfacgk(:,:,ispn,ik), &
          &          apwalm(:,:,:,:,ispn))
        end do

        do is = 1, nspecies
          do ia = 1, natoms (is)
            ias = idxas (ia, is)

            ! generate the density matrix
            call gendmat(.False., .False., 0, lmax, is, ia, &
            &            ngk(:,ik), apwalm, evecfv, evecsv, lmmax, dmat)

            ! convert (l,m) part to an irreducible representation if required
            do ist = 1, nstsv
              do ispn = 1, nspinor
                do jspn = 1, nspinor
                  call zgemm('N', 'N', lmmax, lmmax, lmmax, &
                  &          zone, ulm(:,:,ias), lmmax, &
                  &          dmat(:,:,ispn,jspn,ist), lmmax, &
                  &          zzero, a, lmmax)
                  call zgemm('N', 'C', lmmax, lmmax, lmmax, &
                  &          zone, a, lmmax, ulm(:,:,ias), lmmax, &
                  &          zzero, dmat(:,:,ispn,jspn,ist), lmmax)
                end do
              end do
            end do

            ! determine the band characters from the density matrix
            do ist = 1, nstsv
              do ispn = 1, nspinor
                do lm = 1, lmmax
                  bc(lm,ispn,ias,ist,ik) = dble(dmat(lm,lm,ispn,ispn,ist))
                end do
              end do
            end do

          end do ! ia
        end do ! is
         
      end do ! ik

      deallocate(elm)
      deallocate(ulm)
      deallocate(a)
      deallocate(dmat)
      deallocate(apwalm)
      deallocate(evecfv)
      deallocate(evecsv)

      return
    end subroutine

end subroutine
