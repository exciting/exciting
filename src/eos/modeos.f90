module modeos
! crystal name
character(256) cname
! number of atoms
integer natoms
! EOS type
integer etype
! number of volume points to plot
integer nvplt
! volume plot range
real(8) vplt1,vplt2
! number of energy data points to fit
integer nevpt
! volume and energy data point sets
real(8), allocatable :: vpt(:)
real(8), allocatable :: ept(:)
! maximum number of parameters for an EOS
integer, parameter :: maxparam=100
! number of parameters
integer nparam
! EOS name
character(256) ename(2)
! optimized parameter set
real(8) popt(maxparam)
! parameter names
character(256) pname(maxparam)
! constants
real(8) pi
real(8) aumass_si,aucharge_si,auenergy_si,autime_si,aulength_si
real(8) eps0_si,planck_si,kb_si,kb_au,mu_au,aupress_gpa
real(8) avogadro
! code version
integer version(3)
data version /1,2,0/
end module
