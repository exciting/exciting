!
!
!
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.
!
!
Subroutine timesec (ts)
#ifdef MPI
      Use modmpi
#endif
      Implicit None
! arguments
      Real (8), Intent (Out) :: ts
! local variables
      Integer(8) :: count, count_rate
    
#ifdef MPI
      ts = MPI_wtime()
#else
      Call system_clock (count=count, count_rate=count_rate)
      ts = dble (count) / dble (count_rate)
#endif
      Return
End Subroutine
