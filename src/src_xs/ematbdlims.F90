
subroutine ematbdlims(typ)
  use modmain
  use modxs
  implicit none
  ! arguments
  integer, intent(in) :: typ
  ! set limits for band combinatios
  select case(typ)
  case(-1)
     istlo1=1
     isthi1=nstsv
     istlo2=1
     isthi2=nstsv
     nst1=nstsv
     nst2=nstsv
  case(0)
     istlo1=1
     isthi1=nstval
     istlo2=nstval+1
     isthi2=nstsv
     nst1=nstval
     nst2=nstcon
  case(1)
     istlo1=1
     isthi1=nstval
     istlo2=1
     isthi2=nstval
     nst1=nstval
     nst2=nstval
  case(2)
     istlo1=nstval+1
     isthi1=nstsv
     istlo2=nstval+1
     isthi1=nstsv
     nst1=nstcon
     nst2=nstcon
  end select
end subroutine ematbdlims
