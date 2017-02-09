subroutine xasfinit
	Use modmain
	Use modinput
	Use modxs
	Use modxas
	Implicit none
    ! Deallocation of global arrays
	deallocate(ucore)
	deallocate(ecore)
	deallocate(mj2ml)
	deallocate(preml)

end subroutine xasfinit
