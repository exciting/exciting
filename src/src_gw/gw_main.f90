!BOP
!
! !ROUTINE: gw_main
!
! !INTERFACE:
subroutine gw_main

! !DESCRIPTION:
!
! This is the main program unit. It calculates the corrections to the DFT
! eigenvalues within the GW aproximation...
! 
! !USES:
    use modmain
    use modgw
    use modmpi
!
! !DEFINED PARAMETERS:

    implicit none

!
! !LOCAL VARIABLES:
!
    real(8) :: tstart, tend, t(2)
    character(124)::sbuffer
! !INTRINSIC ROUTINES:

    intrinsic cpu_time

! !EXTERNAL ROUTINES: 


! !REVISION HISTORY:
!       
! Created Nov. 2003 by RGA
! Last Modified: Augus 27th. 2004 by RGA
! Revisited: April 2011 by DIN
!
! !TO DO:
!
!EOP
!BOC

      if (input%gw%taskname.eq.'skip') return

      debug=input%gw%debug
      if(debug)open(55,file='debug.info',action='write')

!---------------------------------------
!     Initializations
!---------------------------------------
      call cpu_time(tstart)
      if (tstart.lt.0.0d0) write(6,*)'warning, tstart < 0'

!     General output file
      fgw=700

#ifdef MPI
      write(sbuffer,*)rank
      open(fgw,file='GWINFO'//trim(adjustl(sbuffer))//'.OUT',action='write',access='Append')
#endif
#ifndef MPI
      open(fgw,file='GWINFO.OUT',action='write',access='Append')
#endif
      call boxmsg(fgw,'=','Main GW output file')
!
!     testid = 16: Calculate bandstructure
!            
      select case(input%gw%taskname)
        case('BAND','band')
          testid=16
          call task_band
          goto 1000
      end select

!     initialise global variables
      call init0
      call init1

!     Parse input data
      call cpu_time(t(1))
      call readingw
      call cpu_time(t(2))
      call write_cputime(fgw,t(2)-t(1),'READINGW')
      
!     Initialize GW global parameters
      call cpu_time(t(1))
      call init_gw
      call barrier
      call cpu_time(t(2))
      call write_cputime(fgw,t(2)-t(1),'INITGW')

      select case(testid)
!
!       testid = 1: Calculate LAPW basis functions for plotting
!
        case (1) 
          call plotlapw
          goto 1000
!
!       testid = 2: Calculate LAPW eigenvectors for plotting
!              
        case(2)
          call plotevec
          goto 1000
            
      end select
 
!     Mixed basis initialization

      call cpu_time(t(1))
      call init_mb
     
      call cpu_time(t(2))
      
      call write_cputime(fgw,t(2)-t(1),'INIMB')
      
      call boxmsg(fgw,'-',"Mixed WF info")
      write(fgw,*)'Max num of IPW for APW:', ngkmax
      write(fgw,*)'Max num of IPW for Mixbasis:', ngqmax
      write(fgw,*)'Max. nr. of mixed functions per atom:', lmixmax
      write(fgw,*)'Mixed basis set size:', lmixmax+ngqmax
     
!     Read the eigenenergies from the  "EIGVAL.OUT" file
      call cpu_time(t(1))      
 
      call readevaldft
     
      call cpu_time(t(2))
      call write_cputime(fgw,t(2)-t(1),'READEVALDFT')
     
      call barrier()
      select case(testid)
!!
!!        testid = 0 or none of below:  Run the GW calculation
!!            
          case (0)
              call gwcycle
!!
!!        testid = 3: Calculate LAPW eigenvectors products for plotting
!!
          case (3) 
              call testprodfun
!
!         testid = 4: Calculate eigenvectors products compared with mix basis
!                    expansion for plotting
          case (4) 
              call testmixfun
!!            
!!        testid = 5 : Integrate eigenvector products directly and as a
!!                      sum of the Minm matrix elements
          case (5) 
              call testmixcomp
!
!         testid = 6 : Test the bare coulomb matrix for various q-points
!            
          case (6)
              call testcoulomb
!!
!!        testid = 7: Calculate the macroscopic dielectric function
!!           
          case (7) 
              call task_emac
!!
!!        testid = 9 Only analityc continuation.
!!            
          case (9)
              call testacont
!!
!!        testid = 10: Calculate the eigenvalues the dielectric function and its inverse
!!            
          case (10) 
              call task_epsev
!            
!!
!!        testid = 11: Calculate the eigenvalues of the screened coulomb
!!                      potential
!!            
!          case (11) 
!            call task_wev
!            
!!
!!        testid = 12: Plot the selfenergy
!!            
          case (12) 
              call testselfeplot
!            
!!
!!        testid = 13: Calculate LDA exchange-correlation matrix elements
!!            
          case (13) 
              call calcvxcnn
!            
!
!         testid = 14: Run the GW calculation exchange only
!            
          case (14)
              call excycle
!
!         testid = 18: Check the rotational matrix for MB functions
!
          case (18)
             call testmbrotmat
!   
      end select 

 1000 call cpu_time(tend)
      if (tend.lt.0.0d0) write(6,*)'warning, tend < 0'
      call linmsg(fgw,'=','')      
      call write_cputime(fgw,tend-tstart,'GW TOTAL CPU TIME')

      if(debug)close(55)
      close(fgw)

#ifdef MPI
      call  cat_logfiles('GWINFO')
#endif 
      return
end subroutine gw_main
!!EOC
