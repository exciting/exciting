
module m_writesumrls
  implicit none
contains

  subroutine writesumrls(iq,s,fn)
    use modmain
    use modxs
    use m_getunit
    use m_tdwriteh
    implicit none
    ! arguments
    integer, intent(in) :: iq
    real(8), intent(in) :: s(3)
    character(*), intent(in) :: fn
    ! local variables
    character(*), parameter :: thisnam = 'writesumrls'
    integer :: n1(1),n,iw

    call getunit(unit1)
    open(unit1,file=trim(fn),action='write')
    ! write parameters as header to file
    call tdwriteh(unit1,iq)
    ! zeroth frequency moment sumrule
    write(unit1,'(a,g18.10,a,g18.10,a)') 'zeroth frequency moment sumrule &
         &(num. val. el.):', s(1), '(', chgval, ')'
    ! first frequency moment sumrule
    write(unit1,'(a,g18.10,a,g18.10,a)') 'first frequency moment sumrule  &
         &(num. val. el.):', s(2), '(', chgval, ')'
    ! one over frequency sumrule
    write(unit1,'(a,g18.10,a,g18.10,a)') 'pi half sumrule                 &
         &(target)       :', s(3), '(', pi/2.d0, ')'
    ! close files
    close(unit1)

  end subroutine writesumrls

end module m_writesumrls
