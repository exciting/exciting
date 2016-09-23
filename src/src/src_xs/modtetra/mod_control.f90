! <sag>
module control
  implicit none
  ! interface parameter
  character(32), save :: tetraifc
  ! default to WIEN2k style (corresponding to original version)
  data tetraifc / 'wien2k' /
  ! level of debug output
  integer, save :: tetradbglv
  ! high default value to mimic original version
  data tetradbglv / 1000 /
  ! handling of pointers (problems with Portland compiler)
  integer, save :: pointerhandling
  ! default is 0 (original version); 1...explicit target assignment in
  ! routine "tetcw" (other routine to be followed)
  data pointerhandling / 0 /
  ! resonance type: 1...resonant term; 2...anti-resonant term; 0...both terms
  integer, save :: restype
  ! initialize with "0" according to original version
  data restype / 0 /
  ! switch wether k+q or k-q should be calculated
  logical, save :: kplusq
  ! initialize with ".false." according to original version
  data kplusq / .false. /
end module control
