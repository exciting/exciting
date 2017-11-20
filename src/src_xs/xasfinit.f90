subroutine xasfinit
  Use modxas, only: ucore, ecore, mj, spj
  Implicit none
    ! Deallocation of global arrays
  deallocate(ucore)
  deallocate(ecore)
  deallocate(spj)
  deallocate(mj)
end subroutine xasfinit
