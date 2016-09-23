
real(8) function idos_exciting(nbd,nik,eband,ntet,tetc,wtet,vt,e)
    use order
    implicit none     
    ! input parameters
    integer(4), intent(in) :: nbd            ! Number of bands
    integer(4), intent(in) :: nik            ! Number of irreducible k-points
    real(8),    intent(in) :: eband(nbd,nik) ! Band energies
    integer(4), intent(in) :: ntet           ! Number of tetrahedra
    integer(4), intent(in) :: tetc(4,*)      ! id. numbers of the corners of the tetrahedra
    integer(4), intent(in) :: wtet(*)        ! weight of each tetrahedron
    real(8), intent(in)    :: vt             ! the volume of the tetrahedra
    real(8), intent(in)    :: e              ! energy
    ! local variables
    integer(4) :: itet, i, ib
    real(8)    :: idos, idt
    real(8)    :: ee(4)
    
    real(8), external :: intdos1t

    idos = 0.d0
    do itet = 1, ntet
      do ib = 1, nbd 
        do i = 1, 4
          ee(i) = eband(ib,tetc(i,itet))
        end do ! i 
        call sort(4,ee)
        idt = intdos1t(ee,e,vt)
        idos = idos+wtet(itet)*idt
      end do ! ib 
    end do ! itet
    
    idos_exciting = idos
    return
end function
          
        
      
