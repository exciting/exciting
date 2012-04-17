module mod_addons

! number of atom classes (non-equivalent atoms)
integer natmcls
! number of atoms in each class
integer, allocatable :: natoms_in_class(:)
! i-th class -> ias mapping
integer, allocatable :: ic2ias(:)
! ias -> ic mapping
integer, allocatable :: ias2ic(:)
! ias -> is mapping
integer, allocatable :: ias2is(:)
! ias -> ia mapping
integer, allocatable :: ias2ia(:)

end module

