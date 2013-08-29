subroutine warning(message)
#ifdef MPI
    use modmpi
#endif
    use m_getunit
        
    implicit none
    character :: message*(*)
    integer   :: unit
#ifdef MPI
    if (rank==0) then
#endif
      call getunit(unit)
      open(unit,File='WARNINGS.OUT',Action='WRITE',Position='APPEND')
      write(unit,*) trim(message)
      close(unit)
#ifdef MPI
    end if
#endif

end subroutine
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
subroutine delete_warnings()
#ifdef MPI
    use modmpi
#endif
    use m_getunit

    implicit none
    integer   :: unit
    logical :: exist
#ifdef MPI    
    if (rank==0) then
#endif
      inquire (File='WARNINGS.OUT', Exist=exist)
      if (exist) Then
        call getunit(unit)
        open(unit, File='WARNINGS.OUT')
        close(unit, Status='DELETE')
      end if
#ifdef MPI
    end if
#endif
end subroutine
