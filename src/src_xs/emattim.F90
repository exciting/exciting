
! Copyright (C) 2006-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module m_emattim

contains

  subroutine emattim(iq,ik,filnam,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,&
       t13,t14,t15,t16,t17)
    use m_getunit
    implicit none
    ! arguments
    integer, intent(in) :: iq,ik
    character(*), intent(in) :: filnam
    real(8), intent(in) :: t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12
    real(8), intent(in) :: t13,t14,t15,t16,t17!,t18,t19,t20,t21
    ! local variables
    integer :: un
    real(8) :: tm
    call getunit(un)
    open(un,file=trim(filnam),action='write',form='formatted', &
         position='append')
    tm=t9+t10+t11+t12
    write(un,*)
    write(un,'("Timings (CPU seconds) for q-point/k-point : ",2i6)') iq,ik
    write(un,'("  initialisation                          : ",f14.4)') t1
    write(un,'("  reading APW coefficients                : ",f14.4)') t2
    write(un,'("  loop over G(+q) vectors                 : ",f14.4)') t3
    write(un,'("    summation wrt Gaunt coefficients      : ",f14.4)') t6
    write(un,'("    muffin-tin contribution               : ",f14.4)') t7
    write(un,'("      APW-APW                             : ",f14.4)') t14
    write(un,'("      APW-lo                              : ",f14.4)') t15
    write(un,'("      lo -APW                             : ",f14.4)') t16
    write(un,'("      lo -lo                              : ",f14.4)') t17
    write(un,'("    interstitial contribution             : ",f14.4)') t8
    write(un,'("    matrix multiplications                : ",f14.4)') tm
    write(un,'("      APW-lo                              : ",f14.4)') t9
    write(un,'("      lo -APW                             : ",f14.4)') t10
    write(un,'("      lo -lo                              : ",f14.4)') t11
    write(un,'("      interstitial                        : ",f14.4)') t12
    write(un,'("    debugging                             : ",f14.4)') t13
    write(un,'("  writing matrix elements to file         : ",f14.4)') t4
    write(un,'("  total                                   : ",f14.4)') t5
    close(un)
  end subroutine emattim

end module m_emattim
