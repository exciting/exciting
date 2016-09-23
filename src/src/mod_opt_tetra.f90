
MODULE mod_opt_tetra
  implicit none
  
  PRIVATE
  
  integer :: ntetra                        ! total number of tetrahedra
  real(8), allocatable, save :: tetra(:,:) ! k-point index of tetrahedra corners                           
  real(8), allocatable, save :: wlsm(:,:)  ! weights for the optimized tetrahedrom method

  PUBLIC opt_tetra_init, opt_tetra_efermi, opt_tetra_occ, opt_tetra_dos

CONTAINS
!
!-------------------------------------------------------------------------------------
subroutine opt_tetra_init(tetra_type, bvec, ngridk, nkp, ikmap)
  !
  ! This routine set the corners and additional points for each tetrahedra
  !
  implicit none
  integer, intent(in) :: tetra_type  ! = 1 for Linear tetrahedron method
                                     ! = 2 for Optimized tetrahedron method
  real(8), intent(in) :: bvec(3,3)   ! Reciplocal lattice vectors [2pi/a]
  integer, intent(in) :: ngridk(3)   ! Grid size
  integer, intent(in) :: nkp         ! number of k-points
  integer, intent(in) :: ikmap(0:ngridk(1)-1,0:ngridk(2)-1,0:ngridk(3)-1)
  ! local variables
  real(8), parameter :: eps = 1.d-6
  integer :: isym, i1, i2, i3, itet, itettot, ii 
  integer :: ik, jk, nk1, nk2, nk3, iv(3)
  integer :: ivvec(3,20,6), divvec(4,4), ivvec0(4)
  integer :: ikv(3)
  real(8) :: l(4), bvec2(3,3), bvec3(3,4)
  !
  nk1 = ngridk(1)
  nk2 = ngridk(2)
  nk3 = ngridk(3)
  !
  ! Take the shortest diagonal line as the "shaft" of tetrahedral devision
  !
  bvec2(1:3,1) = bvec(1:3,1) / dble(nk1)
  bvec2(1:3,2) = bvec(1:3,2) / dble(nk2)
  bvec2(1:3,3) = bvec(1:3,3) / dble(nk3)
  !
  bvec3(1:3,1) = -bvec2(1:3,1) + bvec2(1:3,2) + bvec2(1:3,3)
  bvec3(1:3,2) =  bvec2(1:3,1) - bvec2(1:3,2) + bvec2(1:3,3)
  bvec3(1:3,3) =  bvec2(1:3,1) + bvec2(1:3,2) - bvec2(1:3,3)
  bvec3(1:3,4) =  bvec2(1:3,1) + bvec2(1:3,2) + bvec2(1:3,3)
  !
  do ii = 1, 4
     l(ii) = DOT_PRODUCT(bvec3(1:3, ii), bvec3(1:3, ii))
  end do
  !
  ii = MINLOC(l(1:4),1)
  !
  ivvec0(1:4) = (/ 0, 0, 0, 0 /)
  !
  divvec(1:4,1) = (/ 1, 0, 0, 0 /)
  divvec(1:4,2) = (/ 0, 1, 0, 0 /)
  divvec(1:4,3) = (/ 0, 0, 1, 0 /)
  divvec(1:4,4) = (/ 0, 0, 0, 1 /)
  !
  ivvec0(ii) = 1
  divvec(ii,ii) = -1
  !
  ! Devide a subcell into 6 tetrahedra
  !
  itet = 0
  do i1 = 1, 3
     do i2 = 1, 3
        if(i2 == i1) cycle
        do i3 = 1, 3
           if(i3 == i1 .or. i3 == i2) cycle
           itet = itet + 1
           ivvec(1:3,1,itet) = ivvec0(1:3)
           ivvec(1:3,2,itet) = ivvec(1:3,1,itet) + divvec(1:3,i1)
           ivvec(1:3,3,itet) = ivvec(1:3,2,itet) + divvec(1:3,i2)
           ivvec(1:3,4,itet) = ivvec(1:3,3,itet) + divvec(1:3,i3)
           !
        end do
     end do
  end do
  !
  ! Additional points surrounding the tetrahedron
  !
  ivvec(1:3, 5,1:6) = 2 * ivvec(1:3,1,1:6) - ivvec(1:3,2,1:6)
  ivvec(1:3, 6,1:6) = 2 * ivvec(1:3,2,1:6) - ivvec(1:3,3,1:6)
  ivvec(1:3, 7,1:6) = 2 * ivvec(1:3,3,1:6) - ivvec(1:3,4,1:6)
  ivvec(1:3, 8,1:6) = 2 * ivvec(1:3,4,1:6) - ivvec(1:3,1,1:6)
  !
  ivvec(1:3, 9,1:6) = 2 * ivvec(1:3,1,1:6) - ivvec(1:3,3,1:6)
  ivvec(1:3,10,1:6) = 2 * ivvec(1:3,2,1:6) - ivvec(1:3,4,1:6)
  ivvec(1:3,11,1:6) = 2 * ivvec(1:3,3,1:6) - ivvec(1:3,1,1:6)
  ivvec(1:3,12,1:6) = 2 * ivvec(1:3,4,1:6) - ivvec(1:3,2,1:6)
  !
  ivvec(1:3,13,1:6) = 2 * ivvec(1:3,1,1:6) - ivvec(1:3,4,1:6)
  ivvec(1:3,14,1:6) = 2 * ivvec(1:3,2,1:6) - ivvec(1:3,1,1:6)
  ivvec(1:3,15,1:6) = 2 * ivvec(1:3,3,1:6) - ivvec(1:3,2,1:6)
  ivvec(1:3,16,1:6) = 2 * ivvec(1:3,4,1:6) - ivvec(1:3,3,1:6)
  !
  ivvec(1:3,17,1:6) =  ivvec(1:3,4,1:6) - ivvec(1:3,1,1:6) + ivvec(1:3,2,1:6)
  ivvec(1:3,18,1:6) =  ivvec(1:3,1,1:6) - ivvec(1:3,2,1:6) + ivvec(1:3,3,1:6)
  ivvec(1:3,19,1:6) =  ivvec(1:3,2,1:6) - ivvec(1:3,3,1:6) + ivvec(1:3,4,1:6)
  ivvec(1:3,20,1:6) =  ivvec(1:3,3,1:6) - ivvec(1:3,4,1:6) + ivvec(1:3,1,1:6)
  !
  ! Set the weight for the each tetrahedron method
  !
  if (allocated(wlsm)) deallocate(wlsm)
  allocate(wlsm(4,20))
  !
  if (tetra_type == 1) then
     !
     !write(6,*) "(opt_tetra) Linear tetrahedron method is used."
     !
     wlsm(1:4,1:20) = 0.0d0
     !
     wlsm(1,1) = 1.0d0
     wlsm(2,2) = 1.0d0
     wlsm(3,3) = 1.0d0
     wlsm(4,4) = 1.0d0
     !
  else if(tetra_type == 2) then
     !
     !write(6,*) "(opt_tetra) Optimized tetrahedron method is used."
     !
     wlsm(1, 1: 4) = dble((/1440,    0,   30,    0/))
     wlsm(2, 1: 4) = dble((/   0, 1440,    0,   30/))
     wlsm(3, 1: 4) = dble((/  30,    0, 1440,    0/))
     wlsm(4, 1: 4) = dble((/   0,   30,    0, 1440/))
     !                                              
     wlsm(1, 5: 8) = dble((/ -38,    7,   17,  -28/))
     wlsm(2, 5: 8) = dble((/ -28,  -38,    7,   17/))
     wlsm(3, 5: 8) = dble((/  17,  -28,  -38,    7/))
     wlsm(4, 5: 8) = dble((/   7,   17,  -28,  -38/))
     !                                              
     wlsm(1, 9:12) = dble((/ -56,    9,  -46,    9/))
     wlsm(2, 9:12) = dble((/   9,  -56,    9,  -46/))
     wlsm(3, 9:12) = dble((/ -46,    9,  -56,    9/))
     wlsm(4, 9:12) = dble((/   9,  -46,    9,  -56/))
     !                                              
     wlsm(1,13:16) = dble((/ -38,  -28,   17,    7/))
     wlsm(2,13:16) = dble((/   7,  -38,  -28,   17/))
     wlsm(3,13:16) = dble((/  17,    7,  -38,  -28/))
     wlsm(4,13:16) = dble((/ -28,   17,    7,  -38/))
     !                                              
     wlsm(1,17:20) = dble((/ -18,  -18,   12,  -18/))
     wlsm(2,17:20) = dble((/ -18,  -18,  -18,   12/))
     wlsm(3,17:20) = dble((/  12,  -18,  -18,  -18/))
     wlsm(4,17:20) = dble((/ -18,   12,  -18,  -18/))
     !
     wlsm(1:4,1:20) = wlsm(1:4,1:20) / dble(1260)
     !
  else
     !
     write(6,*) '(opt_tetra_init) tetra_type is wrong !', tetra_type
     !
  end if
  !
  ! construct tetrahedra
  !
  ntetra = 6*nk1*nk2*nk3
  !
  if (allocated(tetra)) deallocate(tetra)
  allocate(tetra(20,ntetra))
  !
  itettot = 0
  do i3 = 0, nk3-1
    do i2 = 0, nk2-1
      do i1 = 0, nk1-1
        do itet = 1, 6
          itettot = itettot + 1
          do ii = 1, 20
            ikv(1:3) = (/i1,i2,i3/)
            ikv(1:3) = ikv(1:3) + ivvec(1:3,ii,itet)
            ikv(1:3) = modulo(ikv(1:3),(/nk1,nk2,nk3/))
            tetra(ii,itettot) = ikmap(ikv(1),ikv(2),ikv(3))
          end do ! ii
        end do
      end do
    end do
  end do
  !
end subroutine opt_tetra_init
!
!----------------------------------------------------------------------------
subroutine opt_tetra_efermi(nelec,nkp,nbnd,ebnd,ef,occ)
  !----------------------------------------------------------------------------
  !
  ! calculate fermi energy by using bisection method
  !
  ! input
  real(8), intent(in) :: nelec           ! number of electrons
  integer, intent(in) :: nkp             ! # of k-points in irr-BZ
  integer, intent(in) :: nbnd            ! # of bands
  real(8), intent(in) :: ebnd(nbnd,nkp)  ! Kohn-Sham energies
  ! output
  real(8), intent(out) :: ef             ! the fermi energy
  real(8), intent(out) :: occ(nbnd,nkp)  ! occupation numbers of each k
  ! local variables
  integer :: ik, iter, maxiter = 500
  real(8) :: elw, eup, sumk
  !
  ! find bounds for the fermi energy
  !
  elw = minval(ebnd(1:nbnd,1:nkp))
  eup = maxval(ebnd(1:nbnd,1:nkp))
  !
  ! bisection method
  !
  do iter = 1, maxiter
     !
     ef = (eup+elw)*0.5d0
     !
     ! calc. # of electrons 
     !
     call opt_tetra_occ(nkp, nbnd, ebnd, ef, occ)
     sumk = sum(occ(1:nbnd,1:nkp))
     !write(*,*) iter, ef, sumk
     !
     ! convergence check
     !
     if (abs(sumk - nelec) < 1.d-6) then
        exit
     else if(sumk < nelec) then
        elw = ef
     else
        eup = ef
     end if
     !
  end do ! iter
  if (iter >= maxiter) write(6,*) "(opt_tetra_weights) Not converged", iter
  !
end subroutine opt_tetra_efermi
!
!--------------------------------------------------------------------------------------
subroutine opt_tetra_occ(nkp, nbnd, ebnd, ef, occ)
  !------------------------------------------------------------------------------------
  !
  ! calculate occupation with given fermi energy
  !
  integer, intent(in)    :: nkp, nbnd
  real(8), intent(in)    :: ef, ebnd(nbnd,nkp)
  real(8), intent(out)   :: occ(nbnd,nkp)
  !
  integer :: ik, it, ib, i, ii
  real(8) :: e(4), c(3), a(4,4)
  real(8) :: wocc(4)
  integer :: idx(4)
  real(8) :: rar(4)
  !
  occ(1:nbnd,1:nkp) = 0.d0
  !
  do it = 1, ntetra
     !
     do ib = 1, nbnd
        !
        e(:) = 0.d0
        do ii = 1, 20
           ik = tetra(ii, it)
           do i = 1, 4
             e(i) = e(i) + wlsm(i,ii) * ebnd(ib,ik)
           end do
        end do
        !
        call sortidx(4, e, idx)
        rar(:) = e(:)
        do i = 1, 4
          e(i) = rar(idx(i))
        end do
        !
        do ii = 1, 4
          a(ii,1:4) = (ef-e(1:4)) / (e(ii)-e(1:4))
        end do
        !
        if( e(1) <= ef .and. ef < e(2) ) then
           !
           c(1) = a(2,1) * a(3,1) * a(4,1) * 0.25d0
           wocc(1) = c(1) * (1.d0 + a(1,2) + a(1,3) + a(1,4))
           wocc(2) = c(1) * a(2,1)
           wocc(3) = c(1) * a(3,1)
           wocc(4) = c(1) * a(4,1)
           !
        else if( e(2) <= ef .and. ef < e(3)) then
           !
           c(1) = a(4,1) * a(3,1) * 0.25d0
           c(2) = a(4,1) * a(3,2) * a(1,3) * 0.25d0
           c(3) = a(4,2) * a(3,2) * a(1,4) * 0.25d0
           !
           wocc(1) = c(1) + (c(1) + c(2)) * a(1,3) + (c(1) + c(2) + c(3)) * a(1,4)
           wocc(2) = c(1) + c(2) + c(3) + (c(2) + c(3)) * a(2,3) + c(3) * a(2,4)
           wocc(3) = (c(1) + c(2)) * a(3,1) + (c(2) + c(3)) * a(3,2)
           wocc(4) = (c(1) + c(2) + c(3)) * a(4,1) + c(3) * a(4,2)
           !
        else if( e(3) <= ef .and. ef < e(4)) then
           !
           c(1) = a(1,4) * a(2,4) * a(3,4)
           wocc(1) = 1.d0 - c(1) * a(1,4)
           wocc(2) = 1.d0 - c(1) * a(2,4)
           wocc(3) = 1.d0 - c(1) * a(3,4)
           wocc(4) = 1.d0 - c(1) * (1.d0 + a(4,1) + a(4,2) + a(4,3))
           wocc(1:4) = wocc(1:4) * 0.25d0
           !
        else if (e(4) <= ef) then
           !
           wocc(1:4) = 0.25d0
           !
        else
           !
           wocc(1:4) = 0.d0
           !
        end if
        !
        do ii = 1, 20
           ik = tetra(ii,it)
           occ(ib,ik) = occ(ib,ik) + dot_product(wlsm(1:4,ii), wocc(idx(1:4)))
        end do
        !
     end do ! ib
     !
  end do ! it
  !
  occ(:,:) = occ(:,:) / dble(ntetra)
  !
  return
end subroutine opt_tetra_occ
!
!--------------------------------------------------------------------------------------
subroutine opt_tetra_dos(nkp, nbnd, ebnd, ef, dfde)
  !------------------------------------------------------------------------------------
  !
  ! calculate occupation with given fermi energy
  !
  integer, intent(in)    :: nkp, nbnd
  real(8), intent(in)    :: ef, ebnd(nbnd,nkp)
  real(8), intent(out)   :: dfde(nbnd,nkp)
  !
  integer :: ik, it, ib, i, ii
  real(8) :: e(4), c(3), a(4,4)
  real(8) :: wdos(4)
  integer :: idx(4)
  real(8) :: rar(4)
  !
  dfde(1:nbnd,1:nkp) = 0.d0
  !
  do it = 1, ntetra
     !
     do ib = 1, nbnd
        !
        e(1:4) = 0.d0
        do ii = 1, 20
           ik = tetra(ii, it)
           e(1:4) = e(1:4) + wlsm(1:4,ii) * ebnd(ib,ik)
        end do
        !
        call sortidx(4, e, idx)
        rar(:) = e(:)
        do i = 1, 4
          e(i) = rar(idx(i))
        end do
        !
        do ii = 1, 4
          a(ii,1:4) = (ef-e(1:4)) / (e(ii)-e(1:4))
        end do
        !
        if( e(1) <= ef .and. ef < e(2) ) then
           !
           c(1) = a(2,1) * a(3,1) * a(4,1) / (ef - e(1))
           wdos(1) = a(1,2) + a(1,3) + a(1,4)
           wdos(2:4) = a(2:4,1)
           wdos(1:4) = wdos(1:4) * c(1)
           !
        else if( e(2) <= ef .and. ef < e(3)) then
           !
           c(1) = a(2,3) * a(3,1) + a(3,2) * a(2,4)
           wdos(1) = a(1,4) * c(1) + a(1,3) * a(3,1) * a(2,3)
           wdos(2) = a(2,3) * c(1) + a(2,4) * a(2,4) * a(3,2)
           wdos(3) = a(3,2) * c(1) + a(3,1) * a(3,1) * a(2,3)
           wdos(4) = a(4,1) * c(1) + a(4,2) * a(2,4) * a(3,2)
           c(1)  = 1d0 / (e(4) - e(1))
           wdos(1:4) = wdos(1:4) * c(1)
           !
        else if( e(3) <= ef .and. ef < e(4)) then
           !
           c(1) = a(1,4) * a(2,4) * a(3,4) / (e(4) - ef)
           wdos(1:3)  = a(1:3,4)
           wdos(4)  = a(4,1) + a(4,2) + a(4,3)
           wdos(1:4) = wdos(1:4) * c(1)
           !
        else if (e(4) <= ef) then
           !
           wdos(1:4) = 0.d0
           !
        else
           !
           wdos(1:4) = 0.d0
           !
        end if
        !
        wdos(1:4) = wdos(1:4) / dble(ntetra)
        !
        do ii = 1, 20
           ik = tetra(ii,it)
           dfde(ib,ik) = dfde(ib,ik) + dot_product(wlsm(1:4,ii), wdos(idx(1:4)))
        end do
        !
     end do ! ib
     !
  end do ! it
  !
end subroutine opt_tetra_dos
!
END MODULE

