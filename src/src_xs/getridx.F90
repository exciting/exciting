
subroutine getridx(procs,set,i,ir)
  ! Get index ir relative to partition determined by overall index i.
  use modmpi, only: firstofset,procofindex
  implicit none
  ! arguments
  integer :: procs
  integer, intent(in) :: set,i
  integer, intent(out) :: ir
  ! executable statements begin
  ir=i-firstofset(procofindex(i,set),set)+1
end subroutine getridx
