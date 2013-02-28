subroutine warning(message)
#ifdef MPI
    use modmpi
#endif
    
    implicit none
    character :: message*(*)
#ifdef MPI
    if (rank==0) then
#endif
      open(100,File='WARNINGS.OUT',Action='WRITE',Position='APPEND')
      write(100,*) trim(message)
      close(100)
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

    implicit none
    logical :: exist
#ifdef MPI    
    if (rank==0) then
#endif
      inquire (File='WARNINGS.OUT', Exist=exist)
      if (exist) Then
        open(100, File='WARNINGS.OUT')
        close(100, Status='DELETE')
      end if
#ifdef MPI
    end if
#endif
end subroutine
