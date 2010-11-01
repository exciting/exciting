
module modinput
use inputdom
implicit none
type origin_type
 real(8)::coord(3)
end type
type point_type
 real(8)::coord(3)
 character(512)::label
end type

type point_type_array
type(point_type),pointer::point
 end type
    type plot1d_type
  type(path_type),pointer::path
end type
type path_type
 integer::steps
 character(512)::outfileprefix
  type(point_type_array),pointer::pointarray(:)
end type
type plot2d_type
  type(parallelogram_type),pointer::parallelogram
end type
type parallelogram_type
 integer::grid(2)
 character(512)::outfileprefix
  type(origin_type),pointer::origin
  type(point_type_array),pointer::pointarray(:)
end type
type plot3d_type
  type(box_type),pointer::box
end type
type box_type
 integer::grid(3)
 character(512)::outfileprefix
  type(origin_type),pointer::origin
  type(point_type_array),pointer::pointarray(:)
end type
type kstlist_type
 integer,pointer::pointstatepair(:,:)
end type
type energywindow_type
 integer::points
 real(8)::intv(2)
end type
type input_type
 character(1024)::xsltpath
 character(1024)::scratchpath
 character(512)::title
  type(structure_type),pointer::structure
  type(groundstate_type),pointer::groundstate
  type(structureoptimization_type),pointer::structureoptimization
  type(properties_type),pointer::properties
  type(phonons_type),pointer::phonons
  type(xs_type),pointer::xs
 character(512)::keywords
end type
type structure_type
 character(1024)::speciespath
 logical::molecule
 real(8)::vacuum
 real(8)::epslat
 logical::autormt
 logical::primcell
 logical::tshift
  type(symmetries_type),pointer::symmetries
  type(crystal_type),pointer::crystal
  type(species_type_array),pointer::speciesarray(:)
end type
type symmetries_type
 character(512)::HermannMauguinSymbol
 character(512)::HallSymbol
 character(512)::SchoenfliesSymbol
 character(512)::spaceGroupNumber
  type(lattice_type),pointer::lattice
  type(WyckoffPositions_type),pointer::WyckoffPositions
end type
type lattice_type
 real(8)::a
 real(8)::b
 real(8)::c
 real(8)::ab
 real(8)::ac
 real(8)::bc
 integer::ncell(3)
 real(8)::scale
 real(8)::stretch(3)
end type
type WyckoffPositions_type
  type(wspecies_type_array),pointer::wspeciesarray(:)
end type
type wspecies_type
 character(512)::speciesfile
  type(wpos_type_array),pointer::wposarray(:)
end type

type wspecies_type_array
type(wspecies_type),pointer::wspecies
 end type
    type wpos_type
 real(8)::coord(3)
end type

type wpos_type_array
type(wpos_type),pointer::wpos
 end type
    type crystal_type
 real(8)::scale
 real(8)::stretch(3)
 real(8),pointer::basevect(:,:)
end type
type species_type
 character(1024)::speciesfile
 character(512)::chemicalSymbol
 integer::atomicNumber
 real(8)::rmt
 character(1024)::href
  type(atom_type_array),pointer::atomarray(:)
  type(LDAplusU_type),pointer::LDAplusU
end type

type species_type_array
type(species_type),pointer::species
 end type
    type atom_type
 real(8)::coord(3)
 real(8)::bfcmt(3)
 real(8)::mommtfix(3)
end type

type atom_type_array
type(atom_type),pointer::atom
 end type
    type LDAplusU_type
 integer::l
 real(8)::U
 real(8)::J
end type
type groundstate_type
 character(512)::do
 integer::ngridk(3)
 real(8)::rgkmax
 real(8)::epspot
 real(8)::epsengy
 real(8)::epsforce
 real(8)::rmtapm(2)
 real(8)::swidth
 character(512)::stype
 integer::stypenumber
 character(512)::findlinentype
 integer::findlinentypenumber
 logical::fermilinengy
 integer::isgkmax
 real(8)::gmaxvr
 integer::nempty
 logical::nosym
 logical::symmorph
 logical::frozencore
 logical::autokpt
 real(8)::radkpt
 integer::nktot
 logical::reducek
 logical::tfibs
 logical::tforce
 integer::lmaxapw
 integer::maxscl
 real(8)::chgexs
 real(8)::deband
 real(8)::epsband
 real(8)::dlinengyfermi
 real(8)::epschg
 real(8)::epsocc
 character(512)::mixer
 integer::mixernumber
 real(8)::beta0
 real(8)::betainc
 real(8)::betadec
 integer::lradstep
 integer::nprad
 character(512)::xctype
 integer::xctypenumber
 character(512)::ldapu
 integer::ldapunumber
 integer::lmaxvr
 real(8)::fracinr
 integer::lmaxinr
 integer::lmaxmat
 real(8)::vkloff(3)
 integer::npsden
 real(8)::cfdamp
 logical::nosource
 logical::tevecsv
 integer::nwrite
 logical::ptnucl
  type(spin_type),pointer::spin
  type(HartreeFock_type),pointer::HartreeFock
  type(solver_type),pointer::solver
  type(OEP_type),pointer::OEP
  type(RDMFT_type),pointer::RDMFT
  type(output_type),pointer::output
  type(libxc_type),pointer::libxc
end type
type spin_type
 real(8)::momfix(3)
 real(8)::bfieldc(3)
 logical::spinorb
 logical::spinsprl
 real(8)::vqlss(3)
 real(8)::taufsm
 real(8)::reducebf
 character(512)::fixspin
 integer::fixspinnumber
end type
type HartreeFock_type
 real(8)::epsengy
end type
type solver_type
 character(512)::type
 integer::typenumber
 logical::packedmatrixstorage
 real(8)::epsarpack
 real(8)::evaltol
end type
type OEP_type
 integer::maxitoep
 real(8)::tauoep(3)
end type
type RDMFT_type
 integer::rdmxctype
 integer::rdmmaxscl
 integer::maxitn
 integer::maxitc
 real(8)::taurdmn
 real(8)::taurdmc
 real(8)::rdmalpha
 real(8)::rdmtemp
end type
type output_type
 character(512)::state
 integer::statenumber
end type
type libxc_type
 character(512)::exchange
 integer::exchangenumber
 character(512)::correlation
 integer::correlationnumber
 character(512)::xc
 integer::xcnumber
end type
type structureoptimization_type
 real(8)::epsforce
 real(8)::tau0atm
 logical::resume
end type
type properties_type
  type(bandstructure_type),pointer::bandstructure
  type(STM_type),pointer::STM
  type(wfplot_type),pointer::wfplot
  type(dos_type),pointer::dos
  type(LSJ_type),pointer::LSJ
  type(masstensor_type),pointer::masstensor
  type(chargedensityplot_type),pointer::chargedensityplot
  type(exccplot_type),pointer::exccplot
  type(elfplot_type),pointer::elfplot
  type(mvecfield_type),pointer::mvecfield
  type(xcmvecfield_type),pointer::xcmvecfield
  type(electricfield_type),pointer::electricfield
  type(gradmvecfield_type),pointer::gradmvecfield
  type(fermisurfaceplot_type),pointer::fermisurfaceplot
  type(EFG_type),pointer::EFG
  type(mossbauer_type),pointer::mossbauer
  type(momentummatrix_type),pointer::momentummatrix
  type(dielectric_type),pointer::dielectric
  type(moke_type),pointer::moke
  type(expiqr_type),pointer::expiqr
  type(elnes_type),pointer::elnes
  type(eliashberg_type),pointer::eliashberg
end type
type bandstructure_type
 real(8)::scissor
 logical::character
  type(plot1d_type),pointer::plot1d
end type
type STM_type
  type(plot2d_type),pointer::plot2d
end type
type wfplot_type
  type(kstlist_type),pointer::kstlist
  type(plot1d_type),pointer::plot1d
  type(plot2d_type),pointer::plot2d
  type(plot3d_type),pointer::plot3d
end type
type dos_type
 real(8)::sqados(3)
 logical::lmirep
 integer::nwdos
 integer::ngrdos
 integer::nsmdos
 real(8)::winddos(2)
 real(8)::scissor
end type
type LSJ_type
  type(kstlist_type),pointer::kstlist
end type
type masstensor_type
 real(8)::deltaem
 integer::ndspem
 real(8)::vklem(3)
end type
type chargedensityplot_type
  type(plot1d_type),pointer::plot1d
  type(plot2d_type),pointer::plot2d
  type(plot3d_type),pointer::plot3d
end type
type exccplot_type
  type(plot1d_type),pointer::plot1d
  type(plot2d_type),pointer::plot2d
  type(plot3d_type),pointer::plot3d
end type
type elfplot_type
  type(plot1d_type),pointer::plot1d
  type(plot2d_type),pointer::plot2d
  type(plot3d_type),pointer::plot3d
end type
type mvecfield_type
  type(plot2d_type),pointer::plot2d
  type(plot3d_type),pointer::plot3d
end type
type xcmvecfield_type
  type(plot2d_type),pointer::plot2d
  type(plot3d_type),pointer::plot3d
end type
type electricfield_type
  type(plot2d_type),pointer::plot2d
  type(plot3d_type),pointer::plot3d
end type
type gradmvecfield_type
  type(plot1d_type),pointer::plot1d
  type(plot2d_type),pointer::plot2d
  type(plot3d_type),pointer::plot3d
end type
type fermisurfaceplot_type
 integer::nstfsp
 logical::separate
end type

type EFG_type
logical::exists
 end type
    
type mossbauer_type
logical::exists
 end type
    type momentummatrix_type
 logical::fastpmat
end type
type dielectric_type
 real(8)::scissor
 logical::intraband
 logical::usegdft
 integer,pointer::optcomp(:,:)
end type

type moke_type
logical::exists
 end type
    
type expiqr_type
logical::exists
 end type
    type elnes_type
 real(8)::vecql(3)
end type
type eliashberg_type
 real(8)::mustar
end type
type phonons_type
 character(512)::do
 integer::donumber
 integer::ngridq(3)
 logical::reduceq
 real(8)::deltaph
  type(qpointset_type),pointer::qpointset
  type(phonondos_type),pointer::phonondos
  type(phonondispplot_type),pointer::phonondispplot
  type(reformatdynmat_type),pointer::reformatdynmat
  type(interpolate_type),pointer::interpolate
  type(parts_type),pointer::parts
end type
type phonondos_type
 integer::nwdos
 integer::ngrdos
 integer::nsmdos
end type
type phonondispplot_type
  type(plot1d_type),pointer::plot1d
end type

type reformatdynmat_type
logical::exists
 end type
    type interpolate_type
 integer::ngridq(3)
 real(8)::vqloff(3)
 logical::writeeigenvectors
end type
type xs_type
 integer::emattype
 logical::dfoffdiag
 integer::lmaxapwwf
 integer::lmaxemat
 real(8)::emaxdf
 real(8)::broad
 real(8)::epsdfde
 logical::tevout
 character(512)::xstype
 integer::xstypenumber
 logical::fastpmat
 logical::fastemat
 logical::tappinfo
 integer::dbglev
 character(512)::gqmaxtype
 integer::gqmaxtypenumber
 real(8)::gqmax
 logical::nosym
 integer::ngridk(3)
 real(8)::vkloff(3)
 logical::reducek
 integer::ngridq(3)
 logical::reduceq
 real(8)::rgkmax
 real(8)::swidth
 integer::lmaxapw
 integer::lmaxmat
 integer::nempty
 real(8)::scissor
  type(tddft_type),pointer::tddft
  type(screening_type),pointer::screening
  type(BSE_type),pointer::BSE
  type(transitions_type),pointer::transitions
  type(qpointset_type),pointer::qpointset
  type(tetra_type),pointer::tetra
  type(energywindow_type),pointer::energywindow
  type(plan_type),pointer::plan
end type
type tddft_type
 logical::intraband
 logical::torddf
 logical::tordfxc
 logical::aresdf
 logical::aresfxc
 real(8)::fxcbsesplit
 logical::acont
 integer::nwacont
 logical::lindhard
 logical::kerndiag
 integer::lmaxalda
 real(8)::alphalrc
 real(8)::alphalrcdyn
 real(8)::betalrcdyn
 integer::mdfqtype
 character(512)::fxctype
 integer::fxctypenumber
 character(512)::do
 integer::donumber
end type
type screening_type
 character(512)::do
 integer::donumber
 logical::nosym
 integer::ngridk(3)
 logical::reducek
 real(8)::vkloff(3)
 real(8)::rgkmax
 integer::nempty
 character(512)::screentype
 integer::screentypenumber
end type
type BSE_type
 logical::nosym
 logical::reducek
 real(8)::vkloff(3)
 real(8)::rgkmax
 integer::scrherm
 logical::fbzq
 character(512)::sciavtype
 integer::sciavtypenumber
 logical::sciavbd
 logical::sciavqhd
 logical::sciavqwg
 logical::sciavqbd
 logical::bsedirsing
 integer::lmaxdielt
 integer::nleblaik
 integer::nexcitmax
 integer::nstlbsemat(4)
 integer::nstlbse(4)
 logical::aresbse
 character(512)::bsetype
 integer::bsetypenumber
end type
type transitions_type
  type(individual_type),pointer::individual
  type(ranges_type),pointer::ranges
  type(lists_type),pointer::lists
end type
type individual_type
  type(trans_type_array),pointer::transarray(:)
end type
type trans_type
 character(512)::action
 integer::kpointnumber
 integer::initial
 integer::final
end type

type trans_type_array
type(trans_type),pointer::trans
 end type
    type ranges_type
  type(range_type_array),pointer::rangearray(:)
end type
type range_type
 character(512)::action
 character(512)::statestype
 integer::kpointnumber
 integer::start
 integer::stop
end type

type range_type_array
type(range_type),pointer::range
 end type
    type lists_type
  type(istate_type_array),pointer::istatearray(:)
end type
type istate_type
 character(512)::action
 character(512)::statestype
 integer::kpointnumber
 integer::state
end type

type istate_type_array
type(istate_type),pointer::istate
 end type
    type tetra_type
 logical::tetraocc
 logical::tetradf
 logical::kordexc
 logical::cw1k
 integer::qweights
end type
type plan_type
  type(doonly_type_array),pointer::doonlyarray(:)
end type
type doonly_type
 character(512)::task
 integer::tasknumber
end type

type doonly_type_array
type(doonly_type),pointer::doonly
 end type
    type qpointset_type
 real(8),pointer::qpoint(:,:)
end type
type parts_type
  type(dopart_type_array),pointer::dopartarray(:)
end type
type dopart_type
 character(512)::id
end type

type dopart_type_array
type(dopart_type),pointer::dopart
 end type
    
   type(input_type)::input
contains

function getstructorigin(thisnode)

implicit none
type(Node),pointer::thisnode
type(origin_type),pointer::getstructorigin
type(Node),pointer::np


integer::len=1,i=0
allocate(getstructorigin)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at origin"
#endif
      
nullify(np)  
np=>getAttributeNode(thisnode,"coord")
if(associated(np)) then
       call extractDataAttribute(thisnode,"coord",getstructorigin%coord)
       call removeAttribute(thisnode,"coord")  
endif

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructpoint(thisnode)

implicit none
type(Node),pointer::thisnode
type(point_type),pointer::getstructpoint
type(Node),pointer::np


integer::len=1,i=0
allocate(getstructpoint)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at point"
#endif
      
nullify(np)  
np=>getAttributeNode(thisnode,"coord")
if(associated(np)) then
       call extractDataAttribute(thisnode,"coord",getstructpoint%coord)
       call removeAttribute(thisnode,"coord")  
        else
        write(*,*)"Parser ERROR: The element 'point' requires the attribute 'coord' to be defined."
        write(*,*)"stopped"
        stop
        
endif

nullify(np)  
np=>getAttributeNode(thisnode,"label")
getstructpoint%label= ""
if(associated(np)) then
       call extractDataAttribute(thisnode,"label",getstructpoint%label)
       call removeAttribute(thisnode,"label")  
endif

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructplot1d(thisnode)

implicit none
type(Node),pointer::thisnode
type(plot1d_type),pointer::getstructplot1d

integer::len=1,i=0
allocate(getstructplot1d)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at plot1d"
#endif
      
            len= countChildEmentsWithName(thisnode,"path")
getstructplot1d%path=>null()
Do i=0,len-1
getstructplot1d%path=>getstructpath(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"path"),0)) ) 
enddo

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructpath(thisnode)

implicit none
type(Node),pointer::thisnode
type(path_type),pointer::getstructpath
type(Node),pointer::np


integer::len=1,i=0
allocate(getstructpath)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at path"
#endif
      
nullify(np)  
np=>getAttributeNode(thisnode,"steps")
if(associated(np)) then
       call extractDataAttribute(thisnode,"steps",getstructpath%steps)
       call removeAttribute(thisnode,"steps")  
        else
        write(*,*)"Parser ERROR: The element 'path' requires the attribute 'steps' to be defined."
        write(*,*)"stopped"
        stop
        
endif

nullify(np)  
np=>getAttributeNode(thisnode,"outfileprefix")
if(associated(np)) then
       call extractDataAttribute(thisnode,"outfileprefix",getstructpath%outfileprefix)
       call removeAttribute(thisnode,"outfileprefix")  
endif

            len= countChildEmentsWithName(thisnode,"point")
     
allocate(getstructpath%pointarray(len))
Do i=0,len-1
getstructpath%pointarray(i+1)%point=>getstructpoint(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"point"),0)) ) 
enddo

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructplot2d(thisnode)

implicit none
type(Node),pointer::thisnode
type(plot2d_type),pointer::getstructplot2d

integer::len=1,i=0
allocate(getstructplot2d)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at plot2d"
#endif
      
            len= countChildEmentsWithName(thisnode,"parallelogram")
getstructplot2d%parallelogram=>null()
Do i=0,len-1
getstructplot2d%parallelogram=>getstructparallelogram(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"parallelogram"),0)) ) 
enddo

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructparallelogram(thisnode)

implicit none
type(Node),pointer::thisnode
type(parallelogram_type),pointer::getstructparallelogram
type(Node),pointer::np


integer::len=1,i=0
allocate(getstructparallelogram)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at parallelogram"
#endif
      
nullify(np)  
np=>getAttributeNode(thisnode,"grid")
if(associated(np)) then
       call extractDataAttribute(thisnode,"grid",getstructparallelogram%grid)
       call removeAttribute(thisnode,"grid")  
        else
        write(*,*)"Parser ERROR: The element 'parallelogram' requires the attribute 'grid' to be defined."
        write(*,*)"stopped"
        stop
        
endif

nullify(np)  
np=>getAttributeNode(thisnode,"outfileprefix")
if(associated(np)) then
       call extractDataAttribute(thisnode,"outfileprefix",getstructparallelogram%outfileprefix)
       call removeAttribute(thisnode,"outfileprefix")  
endif

            len= countChildEmentsWithName(thisnode,"origin")
getstructparallelogram%origin=>null()
Do i=0,len-1
getstructparallelogram%origin=>getstructorigin(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"origin"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"point")
     
allocate(getstructparallelogram%pointarray(len))
Do i=0,len-1
getstructparallelogram%pointarray(i+1)%point=>getstructpoint(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"point"),0)) ) 
enddo

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructplot3d(thisnode)

implicit none
type(Node),pointer::thisnode
type(plot3d_type),pointer::getstructplot3d

integer::len=1,i=0
allocate(getstructplot3d)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at plot3d"
#endif
      
            len= countChildEmentsWithName(thisnode,"box")
getstructplot3d%box=>null()
Do i=0,len-1
getstructplot3d%box=>getstructbox(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"box"),0)) ) 
enddo

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructbox(thisnode)

implicit none
type(Node),pointer::thisnode
type(box_type),pointer::getstructbox
type(Node),pointer::np


integer::len=1,i=0
allocate(getstructbox)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at box"
#endif
      
nullify(np)  
np=>getAttributeNode(thisnode,"grid")
if(associated(np)) then
       call extractDataAttribute(thisnode,"grid",getstructbox%grid)
       call removeAttribute(thisnode,"grid")  
        else
        write(*,*)"Parser ERROR: The element 'box' requires the attribute 'grid' to be defined."
        write(*,*)"stopped"
        stop
        
endif

nullify(np)  
np=>getAttributeNode(thisnode,"outfileprefix")
if(associated(np)) then
       call extractDataAttribute(thisnode,"outfileprefix",getstructbox%outfileprefix)
       call removeAttribute(thisnode,"outfileprefix")  
endif

            len= countChildEmentsWithName(thisnode,"origin")
getstructbox%origin=>null()
Do i=0,len-1
getstructbox%origin=>getstructorigin(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"origin"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"point")
     
allocate(getstructbox%pointarray(len))
Do i=0,len-1
getstructbox%pointarray(i+1)%point=>getstructpoint(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"point"),0)) ) 
enddo

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructkstlist(thisnode)

implicit none
type(Node),pointer::thisnode
type(kstlist_type),pointer::getstructkstlist

integer::len=1,i=0
allocate(getstructkstlist)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at kstlist"
#endif
      
      len= countChildEmentsWithName (thisnode,"pointstatepair")           
allocate(getstructkstlist%pointstatepair(2,len))
Do i=1,len

getstructkstlist%pointstatepair(:,i)=getvalueofpointstatepair(&
      removechild(thisnode,item(getElementsByTagname(thisnode,&
      "pointstatepair"),0)))
end do

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructenergywindow(thisnode)

implicit none
type(Node),pointer::thisnode
type(energywindow_type),pointer::getstructenergywindow
type(Node),pointer::np


integer::len=1,i=0
allocate(getstructenergywindow)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at energywindow"
#endif
      
nullify(np)  
np=>getAttributeNode(thisnode,"points")
getstructenergywindow%points=500
if(associated(np)) then
       call extractDataAttribute(thisnode,"points",getstructenergywindow%points)
       call removeAttribute(thisnode,"points")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"intv")
getstructenergywindow%intv=(/-0.5d0,0.5d0/)
if(associated(np)) then
       call extractDataAttribute(thisnode,"intv",getstructenergywindow%intv)
       call removeAttribute(thisnode,"intv")  
endif

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructinput(thisnode)

implicit none
type(Node),pointer::thisnode
type(input_type),pointer::getstructinput
type(Node),pointer::np


integer::len=1,i=0
allocate(getstructinput)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at input"
#endif
      
nullify(np)  
np=>getAttributeNode(thisnode,"xsltpath")
getstructinput%xsltpath= "http://xml.exciting-code.org"
if(associated(np)) then
       call extractDataAttribute(thisnode,"xsltpath",getstructinput%xsltpath)
       call removeAttribute(thisnode,"xsltpath")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"scratchpath")
if(associated(np)) then
       call extractDataAttribute(thisnode,"scratchpath",getstructinput%scratchpath)
       call removeAttribute(thisnode,"scratchpath")  
endif

            len= countChildEmentsWithName(thisnode,"structure")
getstructinput%structure=>null()
Do i=0,len-1
getstructinput%structure=>getstructstructure(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"structure"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"groundstate")
getstructinput%groundstate=>null()
Do i=0,len-1
getstructinput%groundstate=>getstructgroundstate(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"groundstate"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"structureoptimization")
getstructinput%structureoptimization=>null()
Do i=0,len-1
getstructinput%structureoptimization=>getstructstructureoptimization(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"structureoptimization"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"properties")
getstructinput%properties=>null()
Do i=0,len-1
getstructinput%properties=>getstructproperties(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"properties"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"phonons")
getstructinput%phonons=>null()
Do i=0,len-1
getstructinput%phonons=>getstructphonons(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"phonons"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"xs")
getstructinput%xs=>null()
Do i=0,len-1
getstructinput%xs=>getstructxs(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"xs"),0)) ) 
enddo

      len= countChildEmentsWithName (thisnode,"title")
Do i=1,len

getstructinput%title=getvalueoftitle(&
      removechild(thisnode,item(getElementsByTagname(thisnode,&
      "title"),0)))
end do

      len= countChildEmentsWithName (thisnode,"keywords")
Do i=1,len

getstructinput%keywords=getvalueofkeywords(&
      removechild(thisnode,item(getElementsByTagname(thisnode,&
      "keywords"),0)))
end do

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructstructure(thisnode)

implicit none
type(Node),pointer::thisnode
type(structure_type),pointer::getstructstructure
type(Node),pointer::np


integer::len=1,i=0
allocate(getstructstructure)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at structure"
#endif
      
nullify(np)  
np=>getAttributeNode(thisnode,"speciespath")
if(associated(np)) then
       call extractDataAttribute(thisnode,"speciespath",getstructstructure%speciespath)
       call removeAttribute(thisnode,"speciespath")  
        else
        write(*,*)"Parser ERROR: The element 'structure' requires the attribute 'speciespath' to be defined."
        write(*,*)"stopped"
        stop
        
endif

nullify(np)  
np=>getAttributeNode(thisnode,"molecule")
getstructstructure%molecule= .false.
if(associated(np)) then
       call extractDataAttribute(thisnode,"molecule",getstructstructure%molecule)
       call removeAttribute(thisnode,"molecule")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"vacuum")
getstructstructure%vacuum=10.0d0
if(associated(np)) then
       call extractDataAttribute(thisnode,"vacuum",getstructstructure%vacuum)
       call removeAttribute(thisnode,"vacuum")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"epslat")
getstructstructure%epslat=1.0d-6
if(associated(np)) then
       call extractDataAttribute(thisnode,"epslat",getstructstructure%epslat)
       call removeAttribute(thisnode,"epslat")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"autormt")
getstructstructure%autormt= .false.
if(associated(np)) then
       call extractDataAttribute(thisnode,"autormt",getstructstructure%autormt)
       call removeAttribute(thisnode,"autormt")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"primcell")
getstructstructure%primcell= .false.
if(associated(np)) then
       call extractDataAttribute(thisnode,"primcell",getstructstructure%primcell)
       call removeAttribute(thisnode,"primcell")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"tshift")
getstructstructure%tshift= .true.
if(associated(np)) then
       call extractDataAttribute(thisnode,"tshift",getstructstructure%tshift)
       call removeAttribute(thisnode,"tshift")  
endif

            len= countChildEmentsWithName(thisnode,"symmetries")
getstructstructure%symmetries=>null()
Do i=0,len-1
getstructstructure%symmetries=>getstructsymmetries(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"symmetries"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"crystal")
getstructstructure%crystal=>null()
Do i=0,len-1
getstructstructure%crystal=>getstructcrystal(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"crystal"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"species")
     
allocate(getstructstructure%speciesarray(len))
Do i=0,len-1
getstructstructure%speciesarray(i+1)%species=>getstructspecies(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"species"),0)) ) 
enddo

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructsymmetries(thisnode)

implicit none
type(Node),pointer::thisnode
type(symmetries_type),pointer::getstructsymmetries
type(Node),pointer::np


integer::len=1,i=0
allocate(getstructsymmetries)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at symmetries"
#endif
      
nullify(np)  
np=>getAttributeNode(thisnode,"HermannMauguinSymbol")
if(associated(np)) then
       call extractDataAttribute(thisnode,"HermannMauguinSymbol",getstructsymmetries%HermannMauguinSymbol)
       call removeAttribute(thisnode,"HermannMauguinSymbol")  
        else
        write(*,*)"Parser ERROR: The element 'symmetries' requires the attribute 'HermannMauguinSymbol' to be defined."
        write(*,*)"stopped"
        stop
        
endif

nullify(np)  
np=>getAttributeNode(thisnode,"HallSymbol")
if(associated(np)) then
       call extractDataAttribute(thisnode,"HallSymbol",getstructsymmetries%HallSymbol)
       call removeAttribute(thisnode,"HallSymbol")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"SchoenfliesSymbol")
if(associated(np)) then
       call extractDataAttribute(thisnode,"SchoenfliesSymbol",getstructsymmetries%SchoenfliesSymbol)
       call removeAttribute(thisnode,"SchoenfliesSymbol")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"spaceGroupNumber")
if(associated(np)) then
       call extractDataAttribute(thisnode,"spaceGroupNumber",getstructsymmetries%spaceGroupNumber)
       call removeAttribute(thisnode,"spaceGroupNumber")  
endif

            len= countChildEmentsWithName(thisnode,"lattice")
getstructsymmetries%lattice=>null()
Do i=0,len-1
getstructsymmetries%lattice=>getstructlattice(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"lattice"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"WyckoffPositions")
getstructsymmetries%WyckoffPositions=>null()
Do i=0,len-1
getstructsymmetries%WyckoffPositions=>getstructWyckoffPositions(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"WyckoffPositions"),0)) ) 
enddo

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructlattice(thisnode)

implicit none
type(Node),pointer::thisnode
type(lattice_type),pointer::getstructlattice
type(Node),pointer::np


integer::len=1,i=0
allocate(getstructlattice)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at lattice"
#endif
      
nullify(np)  
np=>getAttributeNode(thisnode,"a")
if(associated(np)) then
       call extractDataAttribute(thisnode,"a",getstructlattice%a)
       call removeAttribute(thisnode,"a")  
        else
        write(*,*)"Parser ERROR: The element 'lattice' requires the attribute 'a' to be defined."
        write(*,*)"stopped"
        stop
        
endif

nullify(np)  
np=>getAttributeNode(thisnode,"b")
if(associated(np)) then
       call extractDataAttribute(thisnode,"b",getstructlattice%b)
       call removeAttribute(thisnode,"b")  
        else
        write(*,*)"Parser ERROR: The element 'lattice' requires the attribute 'b' to be defined."
        write(*,*)"stopped"
        stop
        
endif

nullify(np)  
np=>getAttributeNode(thisnode,"c")
if(associated(np)) then
       call extractDataAttribute(thisnode,"c",getstructlattice%c)
       call removeAttribute(thisnode,"c")  
        else
        write(*,*)"Parser ERROR: The element 'lattice' requires the attribute 'c' to be defined."
        write(*,*)"stopped"
        stop
        
endif

nullify(np)  
np=>getAttributeNode(thisnode,"ab")
if(associated(np)) then
       call extractDataAttribute(thisnode,"ab",getstructlattice%ab)
       call removeAttribute(thisnode,"ab")  
        else
        write(*,*)"Parser ERROR: The element 'lattice' requires the attribute 'ab' to be defined."
        write(*,*)"stopped"
        stop
        
endif

nullify(np)  
np=>getAttributeNode(thisnode,"ac")
if(associated(np)) then
       call extractDataAttribute(thisnode,"ac",getstructlattice%ac)
       call removeAttribute(thisnode,"ac")  
        else
        write(*,*)"Parser ERROR: The element 'lattice' requires the attribute 'ac' to be defined."
        write(*,*)"stopped"
        stop
        
endif

nullify(np)  
np=>getAttributeNode(thisnode,"bc")
if(associated(np)) then
       call extractDataAttribute(thisnode,"bc",getstructlattice%bc)
       call removeAttribute(thisnode,"bc")  
        else
        write(*,*)"Parser ERROR: The element 'lattice' requires the attribute 'bc' to be defined."
        write(*,*)"stopped"
        stop
        
endif

nullify(np)  
np=>getAttributeNode(thisnode,"ncell")
getstructlattice%ncell=(/1,1,1/)
if(associated(np)) then
       call extractDataAttribute(thisnode,"ncell",getstructlattice%ncell)
       call removeAttribute(thisnode,"ncell")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"scale")
getstructlattice%scale=1
if(associated(np)) then
       call extractDataAttribute(thisnode,"scale",getstructlattice%scale)
       call removeAttribute(thisnode,"scale")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"stretch")
getstructlattice%stretch=(/1.0d0,1.0d0,1.0d0/)
if(associated(np)) then
       call extractDataAttribute(thisnode,"stretch",getstructlattice%stretch)
       call removeAttribute(thisnode,"stretch")  
endif

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructWyckoffPositions(thisnode)

implicit none
type(Node),pointer::thisnode
type(WyckoffPositions_type),pointer::getstructWyckoffPositions

integer::len=1,i=0
allocate(getstructWyckoffPositions)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at WyckoffPositions"
#endif
      
            len= countChildEmentsWithName(thisnode,"wspecies")
     
allocate(getstructWyckoffPositions%wspeciesarray(len))
Do i=0,len-1
getstructWyckoffPositions%wspeciesarray(i+1)%wspecies=>getstructwspecies(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"wspecies"),0)) ) 
enddo

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructwspecies(thisnode)

implicit none
type(Node),pointer::thisnode
type(wspecies_type),pointer::getstructwspecies
type(Node),pointer::np


integer::len=1,i=0
allocate(getstructwspecies)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at wspecies"
#endif
      
nullify(np)  
np=>getAttributeNode(thisnode,"speciesfile")
if(associated(np)) then
       call extractDataAttribute(thisnode,"speciesfile",getstructwspecies%speciesfile)
       call removeAttribute(thisnode,"speciesfile")  
endif

            len= countChildEmentsWithName(thisnode,"wpos")
     
allocate(getstructwspecies%wposarray(len))
Do i=0,len-1
getstructwspecies%wposarray(i+1)%wpos=>getstructwpos(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"wpos"),0)) ) 
enddo

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructwpos(thisnode)

implicit none
type(Node),pointer::thisnode
type(wpos_type),pointer::getstructwpos
type(Node),pointer::np


integer::len=1,i=0
allocate(getstructwpos)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at wpos"
#endif
      
nullify(np)  
np=>getAttributeNode(thisnode,"coord")
if(associated(np)) then
       call extractDataAttribute(thisnode,"coord",getstructwpos%coord)
       call removeAttribute(thisnode,"coord")  
endif

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructcrystal(thisnode)

implicit none
type(Node),pointer::thisnode
type(crystal_type),pointer::getstructcrystal
type(Node),pointer::np


integer::len=1,i=0
allocate(getstructcrystal)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at crystal"
#endif
      
nullify(np)  
np=>getAttributeNode(thisnode,"scale")
getstructcrystal%scale=1.0d0
if(associated(np)) then
       call extractDataAttribute(thisnode,"scale",getstructcrystal%scale)
       call removeAttribute(thisnode,"scale")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"stretch")
getstructcrystal%stretch=(/1.0d0,1.0d0,1.0d0/)
if(associated(np)) then
       call extractDataAttribute(thisnode,"stretch",getstructcrystal%stretch)
       call removeAttribute(thisnode,"stretch")  
endif

      len= countChildEmentsWithName (thisnode,"basevect")           
allocate(getstructcrystal%basevect(3,len))
Do i=1,len

getstructcrystal%basevect(:,i)=getvalueofbasevect(&
      removechild(thisnode,item(getElementsByTagname(thisnode,&
      "basevect"),0)))
end do

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructspecies(thisnode)

implicit none
type(Node),pointer::thisnode
type(species_type),pointer::getstructspecies
type(Node),pointer::np


integer::len=1,i=0
allocate(getstructspecies)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at species"
#endif
      
nullify(np)  
np=>getAttributeNode(thisnode,"speciesfile")
if(associated(np)) then
       call extractDataAttribute(thisnode,"speciesfile",getstructspecies%speciesfile)
       call removeAttribute(thisnode,"speciesfile")  
        else
        write(*,*)"Parser ERROR: The element 'species' requires the attribute 'speciesfile' to be defined."
        write(*,*)"stopped"
        stop
        
endif

nullify(np)  
np=>getAttributeNode(thisnode,"chemicalSymbol")
getstructspecies%chemicalSymbol= ""
if(associated(np)) then
       call extractDataAttribute(thisnode,"chemicalSymbol",getstructspecies%chemicalSymbol)
       call removeAttribute(thisnode,"chemicalSymbol")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"atomicNumber")
if(associated(np)) then
       call extractDataAttribute(thisnode,"atomicNumber",getstructspecies%atomicNumber)
       call removeAttribute(thisnode,"atomicNumber")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"rmt")
getstructspecies%rmt=-1.0d0
if(associated(np)) then
       call extractDataAttribute(thisnode,"rmt",getstructspecies%rmt)
       call removeAttribute(thisnode,"rmt")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"href")
getstructspecies%href= ""
if(associated(np)) then
       call extractDataAttribute(thisnode,"href",getstructspecies%href)
       call removeAttribute(thisnode,"href")  
endif

            len= countChildEmentsWithName(thisnode,"atom")
     
allocate(getstructspecies%atomarray(len))
Do i=0,len-1
getstructspecies%atomarray(i+1)%atom=>getstructatom(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"atom"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"LDAplusU")
getstructspecies%LDAplusU=>null()
Do i=0,len-1
getstructspecies%LDAplusU=>getstructLDAplusU(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"LDAplusU"),0)) ) 
enddo

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructatom(thisnode)

implicit none
type(Node),pointer::thisnode
type(atom_type),pointer::getstructatom
type(Node),pointer::np


integer::len=1,i=0
allocate(getstructatom)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at atom"
#endif
      
nullify(np)  
np=>getAttributeNode(thisnode,"coord")
if(associated(np)) then
       call extractDataAttribute(thisnode,"coord",getstructatom%coord)
       call removeAttribute(thisnode,"coord")  
        else
        write(*,*)"Parser ERROR: The element 'atom' requires the attribute 'coord' to be defined."
        write(*,*)"stopped"
        stop
        
endif

nullify(np)  
np=>getAttributeNode(thisnode,"bfcmt")
getstructatom%bfcmt=(/0.0d0,0.0d0,0.0d0/)
if(associated(np)) then
       call extractDataAttribute(thisnode,"bfcmt",getstructatom%bfcmt)
       call removeAttribute(thisnode,"bfcmt")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"mommtfix")
getstructatom%mommtfix=(/0.0d0,0.0d0,0.0d0/)
if(associated(np)) then
       call extractDataAttribute(thisnode,"mommtfix",getstructatom%mommtfix)
       call removeAttribute(thisnode,"mommtfix")  
endif

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructLDAplusU(thisnode)

implicit none
type(Node),pointer::thisnode
type(LDAplusU_type),pointer::getstructLDAplusU
type(Node),pointer::np


integer::len=1,i=0
allocate(getstructLDAplusU)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at LDAplusU"
#endif
      
nullify(np)  
np=>getAttributeNode(thisnode,"l")
getstructLDAplusU%l=-1
if(associated(np)) then
       call extractDataAttribute(thisnode,"l",getstructLDAplusU%l)
       call removeAttribute(thisnode,"l")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"U")
getstructLDAplusU%U=0.0d0
if(associated(np)) then
       call extractDataAttribute(thisnode,"U",getstructLDAplusU%U)
       call removeAttribute(thisnode,"U")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"J")
getstructLDAplusU%J=0.0d0
if(associated(np)) then
       call extractDataAttribute(thisnode,"J",getstructLDAplusU%J)
       call removeAttribute(thisnode,"J")  
endif

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructgroundstate(thisnode)

implicit none
type(Node),pointer::thisnode
type(groundstate_type),pointer::getstructgroundstate
type(Node),pointer::np


integer::len=1,i=0
allocate(getstructgroundstate)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at groundstate"
#endif
      
nullify(np)  
np=>getAttributeNode(thisnode,"do")
getstructgroundstate%do= "fromscratch"
if(associated(np)) then
       call extractDataAttribute(thisnode,"do",getstructgroundstate%do)
       call removeAttribute(thisnode,"do")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"ngridk")
if(associated(np)) then
       call extractDataAttribute(thisnode,"ngridk",getstructgroundstate%ngridk)
       call removeAttribute(thisnode,"ngridk")  
        else
        write(*,*)"Parser ERROR: The element 'groundstate' requires the attribute 'ngridk' to be defined."
        write(*,*)"stopped"
        stop
        
endif

nullify(np)  
np=>getAttributeNode(thisnode,"rgkmax")
getstructgroundstate%rgkmax=7.0d0
if(associated(np)) then
       call extractDataAttribute(thisnode,"rgkmax",getstructgroundstate%rgkmax)
       call removeAttribute(thisnode,"rgkmax")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"epspot")
getstructgroundstate%epspot=1.0d-6
if(associated(np)) then
       call extractDataAttribute(thisnode,"epspot",getstructgroundstate%epspot)
       call removeAttribute(thisnode,"epspot")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"epsengy")
getstructgroundstate%epsengy=1.0d-4
if(associated(np)) then
       call extractDataAttribute(thisnode,"epsengy",getstructgroundstate%epsengy)
       call removeAttribute(thisnode,"epsengy")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"epsforce")
getstructgroundstate%epsforce=5.0d-5
if(associated(np)) then
       call extractDataAttribute(thisnode,"epsforce",getstructgroundstate%epsforce)
       call removeAttribute(thisnode,"epsforce")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"rmtapm")
getstructgroundstate%rmtapm=(/0.25d0,0.95d0/)
if(associated(np)) then
       call extractDataAttribute(thisnode,"rmtapm",getstructgroundstate%rmtapm)
       call removeAttribute(thisnode,"rmtapm")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"swidth")
getstructgroundstate%swidth=0.001d0
if(associated(np)) then
       call extractDataAttribute(thisnode,"swidth",getstructgroundstate%swidth)
       call removeAttribute(thisnode,"swidth")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"stype")
getstructgroundstate%stype= "Gaussian"
if(associated(np)) then
       call extractDataAttribute(thisnode,"stype",getstructgroundstate%stype)
       call removeAttribute(thisnode,"stype")  
endif
getstructgroundstate%stypenumber=stringtonumbergroundstatestype(getstructgroundstate%stype)

nullify(np)  
np=>getAttributeNode(thisnode,"findlinentype")
getstructgroundstate%findlinentype= "advanced"
if(associated(np)) then
       call extractDataAttribute(thisnode,"findlinentype",getstructgroundstate%findlinentype)
       call removeAttribute(thisnode,"findlinentype")  
endif
getstructgroundstate%findlinentypenumber=stringtonumbergroundstatefindlinentype(getstructgroundstate%findlinentype)

nullify(np)  
np=>getAttributeNode(thisnode,"fermilinengy")
getstructgroundstate%fermilinengy= .false.
if(associated(np)) then
       call extractDataAttribute(thisnode,"fermilinengy",getstructgroundstate%fermilinengy)
       call removeAttribute(thisnode,"fermilinengy")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"isgkmax")
getstructgroundstate%isgkmax=-1
if(associated(np)) then
       call extractDataAttribute(thisnode,"isgkmax",getstructgroundstate%isgkmax)
       call removeAttribute(thisnode,"isgkmax")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"gmaxvr")
getstructgroundstate%gmaxvr=12.0d0
if(associated(np)) then
       call extractDataAttribute(thisnode,"gmaxvr",getstructgroundstate%gmaxvr)
       call removeAttribute(thisnode,"gmaxvr")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"nempty")
getstructgroundstate%nempty=5
if(associated(np)) then
       call extractDataAttribute(thisnode,"nempty",getstructgroundstate%nempty)
       call removeAttribute(thisnode,"nempty")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"nosym")
getstructgroundstate%nosym= .false.
if(associated(np)) then
       call extractDataAttribute(thisnode,"nosym",getstructgroundstate%nosym)
       call removeAttribute(thisnode,"nosym")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"symmorph")
getstructgroundstate%symmorph= .false.
if(associated(np)) then
       call extractDataAttribute(thisnode,"symmorph",getstructgroundstate%symmorph)
       call removeAttribute(thisnode,"symmorph")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"frozencore")
getstructgroundstate%frozencore= .false.
if(associated(np)) then
       call extractDataAttribute(thisnode,"frozencore",getstructgroundstate%frozencore)
       call removeAttribute(thisnode,"frozencore")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"autokpt")
getstructgroundstate%autokpt= .false.
if(associated(np)) then
       call extractDataAttribute(thisnode,"autokpt",getstructgroundstate%autokpt)
       call removeAttribute(thisnode,"autokpt")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"radkpt")
getstructgroundstate%radkpt=40.0d0
if(associated(np)) then
       call extractDataAttribute(thisnode,"radkpt",getstructgroundstate%radkpt)
       call removeAttribute(thisnode,"radkpt")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"nktot")
getstructgroundstate%nktot=0
if(associated(np)) then
       call extractDataAttribute(thisnode,"nktot",getstructgroundstate%nktot)
       call removeAttribute(thisnode,"nktot")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"reducek")
getstructgroundstate%reducek= .true.
if(associated(np)) then
       call extractDataAttribute(thisnode,"reducek",getstructgroundstate%reducek)
       call removeAttribute(thisnode,"reducek")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"tfibs")
getstructgroundstate%tfibs= .true.
if(associated(np)) then
       call extractDataAttribute(thisnode,"tfibs",getstructgroundstate%tfibs)
       call removeAttribute(thisnode,"tfibs")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"tforce")
getstructgroundstate%tforce= .false.
if(associated(np)) then
       call extractDataAttribute(thisnode,"tforce",getstructgroundstate%tforce)
       call removeAttribute(thisnode,"tforce")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"lmaxapw")
getstructgroundstate%lmaxapw=10
if(associated(np)) then
       call extractDataAttribute(thisnode,"lmaxapw",getstructgroundstate%lmaxapw)
       call removeAttribute(thisnode,"lmaxapw")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"maxscl")
getstructgroundstate%maxscl=200
if(associated(np)) then
       call extractDataAttribute(thisnode,"maxscl",getstructgroundstate%maxscl)
       call removeAttribute(thisnode,"maxscl")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"chgexs")
getstructgroundstate%chgexs=0.0d0
if(associated(np)) then
       call extractDataAttribute(thisnode,"chgexs",getstructgroundstate%chgexs)
       call removeAttribute(thisnode,"chgexs")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"deband")
getstructgroundstate%deband=0.0025d0
if(associated(np)) then
       call extractDataAttribute(thisnode,"deband",getstructgroundstate%deband)
       call removeAttribute(thisnode,"deband")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"epsband")
getstructgroundstate%epsband=1.0d-6
if(associated(np)) then
       call extractDataAttribute(thisnode,"epsband",getstructgroundstate%epsband)
       call removeAttribute(thisnode,"epsband")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"dlinengyfermi")
getstructgroundstate%dlinengyfermi=-0.1d0
if(associated(np)) then
       call extractDataAttribute(thisnode,"dlinengyfermi",getstructgroundstate%dlinengyfermi)
       call removeAttribute(thisnode,"dlinengyfermi")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"epschg")
getstructgroundstate%epschg=1.0d-3
if(associated(np)) then
       call extractDataAttribute(thisnode,"epschg",getstructgroundstate%epschg)
       call removeAttribute(thisnode,"epschg")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"epsocc")
getstructgroundstate%epsocc=1.0d-8
if(associated(np)) then
       call extractDataAttribute(thisnode,"epsocc",getstructgroundstate%epsocc)
       call removeAttribute(thisnode,"epsocc")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"mixer")
getstructgroundstate%mixer= "msec"
if(associated(np)) then
       call extractDataAttribute(thisnode,"mixer",getstructgroundstate%mixer)
       call removeAttribute(thisnode,"mixer")  
endif
getstructgroundstate%mixernumber=stringtonumbergroundstatemixer(getstructgroundstate%mixer)

nullify(np)  
np=>getAttributeNode(thisnode,"beta0")
getstructgroundstate%beta0=0.4d0
if(associated(np)) then
       call extractDataAttribute(thisnode,"beta0",getstructgroundstate%beta0)
       call removeAttribute(thisnode,"beta0")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"betainc")
getstructgroundstate%betainc=1.1d0
if(associated(np)) then
       call extractDataAttribute(thisnode,"betainc",getstructgroundstate%betainc)
       call removeAttribute(thisnode,"betainc")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"betadec")
getstructgroundstate%betadec=0.6d0
if(associated(np)) then
       call extractDataAttribute(thisnode,"betadec",getstructgroundstate%betadec)
       call removeAttribute(thisnode,"betadec")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"lradstep")
getstructgroundstate%lradstep=4
if(associated(np)) then
       call extractDataAttribute(thisnode,"lradstep",getstructgroundstate%lradstep)
       call removeAttribute(thisnode,"lradstep")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"nprad")
getstructgroundstate%nprad=4
if(associated(np)) then
       call extractDataAttribute(thisnode,"nprad",getstructgroundstate%nprad)
       call removeAttribute(thisnode,"nprad")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"xctype")
getstructgroundstate%xctype= "LSDAPerdew-Wang"
if(associated(np)) then
       call extractDataAttribute(thisnode,"xctype",getstructgroundstate%xctype)
       call removeAttribute(thisnode,"xctype")  
endif
getstructgroundstate%xctypenumber=stringtonumbergroundstatexctype(getstructgroundstate%xctype)

nullify(np)  
np=>getAttributeNode(thisnode,"ldapu")
getstructgroundstate%ldapu= "none"
if(associated(np)) then
       call extractDataAttribute(thisnode,"ldapu",getstructgroundstate%ldapu)
       call removeAttribute(thisnode,"ldapu")  
endif
getstructgroundstate%ldapunumber=stringtonumbergroundstateldapu(getstructgroundstate%ldapu)

nullify(np)  
np=>getAttributeNode(thisnode,"lmaxvr")
getstructgroundstate%lmaxvr=6
if(associated(np)) then
       call extractDataAttribute(thisnode,"lmaxvr",getstructgroundstate%lmaxvr)
       call removeAttribute(thisnode,"lmaxvr")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"fracinr")
getstructgroundstate%fracinr=0.25d0
if(associated(np)) then
       call extractDataAttribute(thisnode,"fracinr",getstructgroundstate%fracinr)
       call removeAttribute(thisnode,"fracinr")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"lmaxinr")
getstructgroundstate%lmaxinr=2
if(associated(np)) then
       call extractDataAttribute(thisnode,"lmaxinr",getstructgroundstate%lmaxinr)
       call removeAttribute(thisnode,"lmaxinr")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"lmaxmat")
getstructgroundstate%lmaxmat=5
if(associated(np)) then
       call extractDataAttribute(thisnode,"lmaxmat",getstructgroundstate%lmaxmat)
       call removeAttribute(thisnode,"lmaxmat")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"vkloff")
getstructgroundstate%vkloff=(/0.0d0,0.0d0,0.0d0/)
if(associated(np)) then
       call extractDataAttribute(thisnode,"vkloff",getstructgroundstate%vkloff)
       call removeAttribute(thisnode,"vkloff")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"npsden")
getstructgroundstate%npsden=9
if(associated(np)) then
       call extractDataAttribute(thisnode,"npsden",getstructgroundstate%npsden)
       call removeAttribute(thisnode,"npsden")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"cfdamp")
getstructgroundstate%cfdamp=0.0d0
if(associated(np)) then
       call extractDataAttribute(thisnode,"cfdamp",getstructgroundstate%cfdamp)
       call removeAttribute(thisnode,"cfdamp")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"nosource")
getstructgroundstate%nosource= .false.
if(associated(np)) then
       call extractDataAttribute(thisnode,"nosource",getstructgroundstate%nosource)
       call removeAttribute(thisnode,"nosource")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"tevecsv")
getstructgroundstate%tevecsv= .false.
if(associated(np)) then
       call extractDataAttribute(thisnode,"tevecsv",getstructgroundstate%tevecsv)
       call removeAttribute(thisnode,"tevecsv")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"nwrite")
getstructgroundstate%nwrite=0
if(associated(np)) then
       call extractDataAttribute(thisnode,"nwrite",getstructgroundstate%nwrite)
       call removeAttribute(thisnode,"nwrite")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"ptnucl")
getstructgroundstate%ptnucl= .true.
if(associated(np)) then
       call extractDataAttribute(thisnode,"ptnucl",getstructgroundstate%ptnucl)
       call removeAttribute(thisnode,"ptnucl")  
endif

            len= countChildEmentsWithName(thisnode,"spin")
getstructgroundstate%spin=>null()
Do i=0,len-1
getstructgroundstate%spin=>getstructspin(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"spin"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"HartreeFock")
getstructgroundstate%HartreeFock=>null()
Do i=0,len-1
getstructgroundstate%HartreeFock=>getstructHartreeFock(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"HartreeFock"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"solver")
getstructgroundstate%solver=>null()
Do i=0,len-1
getstructgroundstate%solver=>getstructsolver(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"solver"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"OEP")
getstructgroundstate%OEP=>null()
Do i=0,len-1
getstructgroundstate%OEP=>getstructOEP(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"OEP"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"RDMFT")
getstructgroundstate%RDMFT=>null()
Do i=0,len-1
getstructgroundstate%RDMFT=>getstructRDMFT(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"RDMFT"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"output")
getstructgroundstate%output=>null()
Do i=0,len-1
getstructgroundstate%output=>getstructoutput(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"output"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"libxc")
getstructgroundstate%libxc=>null()
Do i=0,len-1
getstructgroundstate%libxc=>getstructlibxc(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"libxc"),0)) ) 
enddo

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructspin(thisnode)

implicit none
type(Node),pointer::thisnode
type(spin_type),pointer::getstructspin
type(Node),pointer::np


integer::len=1,i=0
allocate(getstructspin)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at spin"
#endif
      
nullify(np)  
np=>getAttributeNode(thisnode,"momfix")
getstructspin%momfix=(/0.0d0,0.0d0,0.0d0/)
if(associated(np)) then
       call extractDataAttribute(thisnode,"momfix",getstructspin%momfix)
       call removeAttribute(thisnode,"momfix")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"bfieldc")
getstructspin%bfieldc=(/0.0d0,0.0d0,0.0d0/)
if(associated(np)) then
       call extractDataAttribute(thisnode,"bfieldc",getstructspin%bfieldc)
       call removeAttribute(thisnode,"bfieldc")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"spinorb")
if(associated(np)) then
       call extractDataAttribute(thisnode,"spinorb",getstructspin%spinorb)
       call removeAttribute(thisnode,"spinorb")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"spinsprl")
getstructspin%spinsprl= .false.
if(associated(np)) then
       call extractDataAttribute(thisnode,"spinsprl",getstructspin%spinsprl)
       call removeAttribute(thisnode,"spinsprl")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"vqlss")
getstructspin%vqlss=(/0.0d0,0.0d0,0.0d0/)
if(associated(np)) then
       call extractDataAttribute(thisnode,"vqlss",getstructspin%vqlss)
       call removeAttribute(thisnode,"vqlss")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"taufsm")
getstructspin%taufsm=0.01d0
if(associated(np)) then
       call extractDataAttribute(thisnode,"taufsm",getstructspin%taufsm)
       call removeAttribute(thisnode,"taufsm")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"reducebf")
getstructspin%reducebf=1.0d0
if(associated(np)) then
       call extractDataAttribute(thisnode,"reducebf",getstructspin%reducebf)
       call removeAttribute(thisnode,"reducebf")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"fixspin")
getstructspin%fixspin= "none"
if(associated(np)) then
       call extractDataAttribute(thisnode,"fixspin",getstructspin%fixspin)
       call removeAttribute(thisnode,"fixspin")  
endif
getstructspin%fixspinnumber=stringtonumberspinfixspin(getstructspin%fixspin)

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructHartreeFock(thisnode)

implicit none
type(Node),pointer::thisnode
type(HartreeFock_type),pointer::getstructHartreeFock
type(Node),pointer::np


integer::len=1,i=0
allocate(getstructHartreeFock)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at HartreeFock"
#endif
      
nullify(np)  
np=>getAttributeNode(thisnode,"epsengy")
getstructHartreeFock%epsengy=1.0d-4
if(associated(np)) then
       call extractDataAttribute(thisnode,"epsengy",getstructHartreeFock%epsengy)
       call removeAttribute(thisnode,"epsengy")  
endif

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructsolver(thisnode)

implicit none
type(Node),pointer::thisnode
type(solver_type),pointer::getstructsolver
type(Node),pointer::np


integer::len=1,i=0
allocate(getstructsolver)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at solver"
#endif
      
nullify(np)  
np=>getAttributeNode(thisnode,"type")
getstructsolver%type= "Lapack"
if(associated(np)) then
       call extractDataAttribute(thisnode,"type",getstructsolver%type)
       call removeAttribute(thisnode,"type")  
endif
getstructsolver%typenumber=stringtonumbersolvertype(getstructsolver%type)

nullify(np)  
np=>getAttributeNode(thisnode,"packedmatrixstorage")
getstructsolver%packedmatrixstorage= .true.
if(associated(np)) then
       call extractDataAttribute(thisnode,"packedmatrixstorage",getstructsolver%packedmatrixstorage)
       call removeAttribute(thisnode,"packedmatrixstorage")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"epsarpack")
getstructsolver%epsarpack=1.0d-8
if(associated(np)) then
       call extractDataAttribute(thisnode,"epsarpack",getstructsolver%epsarpack)
       call removeAttribute(thisnode,"epsarpack")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"evaltol")
getstructsolver%evaltol=1.0d-8
if(associated(np)) then
       call extractDataAttribute(thisnode,"evaltol",getstructsolver%evaltol)
       call removeAttribute(thisnode,"evaltol")  
endif

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructOEP(thisnode)

implicit none
type(Node),pointer::thisnode
type(OEP_type),pointer::getstructOEP
type(Node),pointer::np


integer::len=1,i=0
allocate(getstructOEP)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at OEP"
#endif
      
nullify(np)  
np=>getAttributeNode(thisnode,"maxitoep")
getstructOEP%maxitoep=120
if(associated(np)) then
       call extractDataAttribute(thisnode,"maxitoep",getstructOEP%maxitoep)
       call removeAttribute(thisnode,"maxitoep")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"tauoep")
getstructOEP%tauoep=(/1.0d0,0.2d0,1.5d0/)
if(associated(np)) then
       call extractDataAttribute(thisnode,"tauoep",getstructOEP%tauoep)
       call removeAttribute(thisnode,"tauoep")  
endif

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructRDMFT(thisnode)

implicit none
type(Node),pointer::thisnode
type(RDMFT_type),pointer::getstructRDMFT
type(Node),pointer::np


integer::len=1,i=0
allocate(getstructRDMFT)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at RDMFT"
#endif
      
nullify(np)  
np=>getAttributeNode(thisnode,"rdmxctype")
getstructRDMFT%rdmxctype=2
if(associated(np)) then
       call extractDataAttribute(thisnode,"rdmxctype",getstructRDMFT%rdmxctype)
       call removeAttribute(thisnode,"rdmxctype")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"rdmmaxscl")
getstructRDMFT%rdmmaxscl=1
if(associated(np)) then
       call extractDataAttribute(thisnode,"rdmmaxscl",getstructRDMFT%rdmmaxscl)
       call removeAttribute(thisnode,"rdmmaxscl")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"maxitn")
getstructRDMFT%maxitn=250
if(associated(np)) then
       call extractDataAttribute(thisnode,"maxitn",getstructRDMFT%maxitn)
       call removeAttribute(thisnode,"maxitn")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"maxitc")
getstructRDMFT%maxitc=10
if(associated(np)) then
       call extractDataAttribute(thisnode,"maxitc",getstructRDMFT%maxitc)
       call removeAttribute(thisnode,"maxitc")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"taurdmn")
getstructRDMFT%taurdmn=1.0d0
if(associated(np)) then
       call extractDataAttribute(thisnode,"taurdmn",getstructRDMFT%taurdmn)
       call removeAttribute(thisnode,"taurdmn")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"taurdmc")
getstructRDMFT%taurdmc=0.5d0
if(associated(np)) then
       call extractDataAttribute(thisnode,"taurdmc",getstructRDMFT%taurdmc)
       call removeAttribute(thisnode,"taurdmc")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"rdmalpha")
getstructRDMFT%rdmalpha=0.7d0
if(associated(np)) then
       call extractDataAttribute(thisnode,"rdmalpha",getstructRDMFT%rdmalpha)
       call removeAttribute(thisnode,"rdmalpha")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"rdmtemp")
getstructRDMFT%rdmtemp=0.0d0
if(associated(np)) then
       call extractDataAttribute(thisnode,"rdmtemp",getstructRDMFT%rdmtemp)
       call removeAttribute(thisnode,"rdmtemp")  
endif

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructoutput(thisnode)

implicit none
type(Node),pointer::thisnode
type(output_type),pointer::getstructoutput
type(Node),pointer::np


integer::len=1,i=0
allocate(getstructoutput)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at output"
#endif
      
nullify(np)  
np=>getAttributeNode(thisnode,"state")
getstructoutput%state= "binary"
if(associated(np)) then
       call extractDataAttribute(thisnode,"state",getstructoutput%state)
       call removeAttribute(thisnode,"state")  
endif
getstructoutput%statenumber=stringtonumberoutputstate(getstructoutput%state)

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructlibxc(thisnode)

implicit none
type(Node),pointer::thisnode
type(libxc_type),pointer::getstructlibxc
type(Node),pointer::np


integer::len=1,i=0
allocate(getstructlibxc)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at libxc"
#endif
      
nullify(np)  
np=>getAttributeNode(thisnode,"exchange")
getstructlibxc%exchange= "XC_GGA_X_PBE"
if(associated(np)) then
       call extractDataAttribute(thisnode,"exchange",getstructlibxc%exchange)
       call removeAttribute(thisnode,"exchange")  
endif
getstructlibxc%exchangenumber=stringtonumberlibxcexchange(getstructlibxc%exchange)

nullify(np)  
np=>getAttributeNode(thisnode,"correlation")
getstructlibxc%correlation= "XC_GGA_C_PBE"
if(associated(np)) then
       call extractDataAttribute(thisnode,"correlation",getstructlibxc%correlation)
       call removeAttribute(thisnode,"correlation")  
endif
getstructlibxc%correlationnumber=stringtonumberlibxccorrelation(getstructlibxc%correlation)

nullify(np)  
np=>getAttributeNode(thisnode,"xc")
getstructlibxc%xc= "none"
if(associated(np)) then
       call extractDataAttribute(thisnode,"xc",getstructlibxc%xc)
       call removeAttribute(thisnode,"xc")  
endif
getstructlibxc%xcnumber=stringtonumberlibxcxc(getstructlibxc%xc)

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructstructureoptimization(thisnode)

implicit none
type(Node),pointer::thisnode
type(structureoptimization_type),pointer::getstructstructureoptimization
type(Node),pointer::np


integer::len=1,i=0
allocate(getstructstructureoptimization)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at structureoptimization"
#endif
      
nullify(np)  
np=>getAttributeNode(thisnode,"epsforce")
getstructstructureoptimization%epsforce=5.0d-5
if(associated(np)) then
       call extractDataAttribute(thisnode,"epsforce",getstructstructureoptimization%epsforce)
       call removeAttribute(thisnode,"epsforce")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"tau0atm")
getstructstructureoptimization%tau0atm=0.2d0
if(associated(np)) then
       call extractDataAttribute(thisnode,"tau0atm",getstructstructureoptimization%tau0atm)
       call removeAttribute(thisnode,"tau0atm")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"resume")
getstructstructureoptimization%resume= .false.
if(associated(np)) then
       call extractDataAttribute(thisnode,"resume",getstructstructureoptimization%resume)
       call removeAttribute(thisnode,"resume")  
endif

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructproperties(thisnode)

implicit none
type(Node),pointer::thisnode
type(properties_type),pointer::getstructproperties

integer::len=1,i=0
allocate(getstructproperties)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at properties"
#endif
      
            len= countChildEmentsWithName(thisnode,"bandstructure")
getstructproperties%bandstructure=>null()
Do i=0,len-1
getstructproperties%bandstructure=>getstructbandstructure(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"bandstructure"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"STM")
getstructproperties%STM=>null()
Do i=0,len-1
getstructproperties%STM=>getstructSTM(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"STM"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"wfplot")
getstructproperties%wfplot=>null()
Do i=0,len-1
getstructproperties%wfplot=>getstructwfplot(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"wfplot"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"dos")
getstructproperties%dos=>null()
Do i=0,len-1
getstructproperties%dos=>getstructdos(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"dos"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"LSJ")
getstructproperties%LSJ=>null()
Do i=0,len-1
getstructproperties%LSJ=>getstructLSJ(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"LSJ"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"masstensor")
getstructproperties%masstensor=>null()
Do i=0,len-1
getstructproperties%masstensor=>getstructmasstensor(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"masstensor"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"chargedensityplot")
getstructproperties%chargedensityplot=>null()
Do i=0,len-1
getstructproperties%chargedensityplot=>getstructchargedensityplot(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"chargedensityplot"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"exccplot")
getstructproperties%exccplot=>null()
Do i=0,len-1
getstructproperties%exccplot=>getstructexccplot(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"exccplot"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"elfplot")
getstructproperties%elfplot=>null()
Do i=0,len-1
getstructproperties%elfplot=>getstructelfplot(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"elfplot"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"mvecfield")
getstructproperties%mvecfield=>null()
Do i=0,len-1
getstructproperties%mvecfield=>getstructmvecfield(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"mvecfield"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"xcmvecfield")
getstructproperties%xcmvecfield=>null()
Do i=0,len-1
getstructproperties%xcmvecfield=>getstructxcmvecfield(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"xcmvecfield"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"electricfield")
getstructproperties%electricfield=>null()
Do i=0,len-1
getstructproperties%electricfield=>getstructelectricfield(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"electricfield"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"gradmvecfield")
getstructproperties%gradmvecfield=>null()
Do i=0,len-1
getstructproperties%gradmvecfield=>getstructgradmvecfield(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"gradmvecfield"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"fermisurfaceplot")
getstructproperties%fermisurfaceplot=>null()
Do i=0,len-1
getstructproperties%fermisurfaceplot=>getstructfermisurfaceplot(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"fermisurfaceplot"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"EFG")
getstructproperties%EFG=>null()
Do i=0,len-1
getstructproperties%EFG=>getstructEFG(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"EFG"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"mossbauer")
getstructproperties%mossbauer=>null()
Do i=0,len-1
getstructproperties%mossbauer=>getstructmossbauer(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"mossbauer"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"momentummatrix")
getstructproperties%momentummatrix=>null()
Do i=0,len-1
getstructproperties%momentummatrix=>getstructmomentummatrix(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"momentummatrix"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"dielectric")
getstructproperties%dielectric=>null()
Do i=0,len-1
getstructproperties%dielectric=>getstructdielectric(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"dielectric"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"moke")
getstructproperties%moke=>null()
Do i=0,len-1
getstructproperties%moke=>getstructmoke(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"moke"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"expiqr")
getstructproperties%expiqr=>null()
Do i=0,len-1
getstructproperties%expiqr=>getstructexpiqr(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"expiqr"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"elnes")
getstructproperties%elnes=>null()
Do i=0,len-1
getstructproperties%elnes=>getstructelnes(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"elnes"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"eliashberg")
getstructproperties%eliashberg=>null()
Do i=0,len-1
getstructproperties%eliashberg=>getstructeliashberg(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"eliashberg"),0)) ) 
enddo

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructbandstructure(thisnode)

implicit none
type(Node),pointer::thisnode
type(bandstructure_type),pointer::getstructbandstructure
type(Node),pointer::np


integer::len=1,i=0
allocate(getstructbandstructure)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at bandstructure"
#endif
      
nullify(np)  
np=>getAttributeNode(thisnode,"scissor")
getstructbandstructure%scissor=0.0d0
if(associated(np)) then
       call extractDataAttribute(thisnode,"scissor",getstructbandstructure%scissor)
       call removeAttribute(thisnode,"scissor")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"character")
getstructbandstructure%character= .false.
if(associated(np)) then
       call extractDataAttribute(thisnode,"character",getstructbandstructure%character)
       call removeAttribute(thisnode,"character")  
endif

            len= countChildEmentsWithName(thisnode,"plot1d")
getstructbandstructure%plot1d=>null()
Do i=0,len-1
getstructbandstructure%plot1d=>getstructplot1d(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"plot1d"),0)) ) 
enddo

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructSTM(thisnode)

implicit none
type(Node),pointer::thisnode
type(STM_type),pointer::getstructSTM

integer::len=1,i=0
allocate(getstructSTM)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at STM"
#endif
      
            len= countChildEmentsWithName(thisnode,"plot2d")
getstructSTM%plot2d=>null()
Do i=0,len-1
getstructSTM%plot2d=>getstructplot2d(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"plot2d"),0)) ) 
enddo

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructwfplot(thisnode)

implicit none
type(Node),pointer::thisnode
type(wfplot_type),pointer::getstructwfplot

integer::len=1,i=0
allocate(getstructwfplot)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at wfplot"
#endif
      
            len= countChildEmentsWithName(thisnode,"kstlist")
getstructwfplot%kstlist=>null()
Do i=0,len-1
getstructwfplot%kstlist=>getstructkstlist(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"kstlist"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"plot1d")
getstructwfplot%plot1d=>null()
Do i=0,len-1
getstructwfplot%plot1d=>getstructplot1d(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"plot1d"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"plot2d")
getstructwfplot%plot2d=>null()
Do i=0,len-1
getstructwfplot%plot2d=>getstructplot2d(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"plot2d"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"plot3d")
getstructwfplot%plot3d=>null()
Do i=0,len-1
getstructwfplot%plot3d=>getstructplot3d(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"plot3d"),0)) ) 
enddo

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructdos(thisnode)

implicit none
type(Node),pointer::thisnode
type(dos_type),pointer::getstructdos
type(Node),pointer::np


integer::len=1,i=0
allocate(getstructdos)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at dos"
#endif
      
nullify(np)  
np=>getAttributeNode(thisnode,"sqados")
getstructdos%sqados=(/0.0d0,0.0d0,1.0d0/)
if(associated(np)) then
       call extractDataAttribute(thisnode,"sqados",getstructdos%sqados)
       call removeAttribute(thisnode,"sqados")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"lmirep")
getstructdos%lmirep= .false.
if(associated(np)) then
       call extractDataAttribute(thisnode,"lmirep",getstructdos%lmirep)
       call removeAttribute(thisnode,"lmirep")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"nwdos")
getstructdos%nwdos=500
if(associated(np)) then
       call extractDataAttribute(thisnode,"nwdos",getstructdos%nwdos)
       call removeAttribute(thisnode,"nwdos")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"ngrdos")
getstructdos%ngrdos=100
if(associated(np)) then
       call extractDataAttribute(thisnode,"ngrdos",getstructdos%ngrdos)
       call removeAttribute(thisnode,"ngrdos")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"nsmdos")
getstructdos%nsmdos=0
if(associated(np)) then
       call extractDataAttribute(thisnode,"nsmdos",getstructdos%nsmdos)
       call removeAttribute(thisnode,"nsmdos")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"winddos")
getstructdos%winddos=(/-0.5d0,0.5d0/)
if(associated(np)) then
       call extractDataAttribute(thisnode,"winddos",getstructdos%winddos)
       call removeAttribute(thisnode,"winddos")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"scissor")
getstructdos%scissor=0.0d0
if(associated(np)) then
       call extractDataAttribute(thisnode,"scissor",getstructdos%scissor)
       call removeAttribute(thisnode,"scissor")  
endif

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructLSJ(thisnode)

implicit none
type(Node),pointer::thisnode
type(LSJ_type),pointer::getstructLSJ

integer::len=1,i=0
allocate(getstructLSJ)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at LSJ"
#endif
      
            len= countChildEmentsWithName(thisnode,"kstlist")
getstructLSJ%kstlist=>null()
Do i=0,len-1
getstructLSJ%kstlist=>getstructkstlist(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"kstlist"),0)) ) 
enddo

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructmasstensor(thisnode)

implicit none
type(Node),pointer::thisnode
type(masstensor_type),pointer::getstructmasstensor
type(Node),pointer::np


integer::len=1,i=0
allocate(getstructmasstensor)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at masstensor"
#endif
      
nullify(np)  
np=>getAttributeNode(thisnode,"deltaem")
getstructmasstensor%deltaem=0.025d0
if(associated(np)) then
       call extractDataAttribute(thisnode,"deltaem",getstructmasstensor%deltaem)
       call removeAttribute(thisnode,"deltaem")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"ndspem")
getstructmasstensor%ndspem=1
if(associated(np)) then
       call extractDataAttribute(thisnode,"ndspem",getstructmasstensor%ndspem)
       call removeAttribute(thisnode,"ndspem")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"vklem")
getstructmasstensor%vklem=(/0.0d0,0.0d0,0.0d0/)
if(associated(np)) then
       call extractDataAttribute(thisnode,"vklem",getstructmasstensor%vklem)
       call removeAttribute(thisnode,"vklem")  
endif

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructchargedensityplot(thisnode)

implicit none
type(Node),pointer::thisnode
type(chargedensityplot_type),pointer::getstructchargedensityplot

integer::len=1,i=0
allocate(getstructchargedensityplot)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at chargedensityplot"
#endif
      
            len= countChildEmentsWithName(thisnode,"plot1d")
getstructchargedensityplot%plot1d=>null()
Do i=0,len-1
getstructchargedensityplot%plot1d=>getstructplot1d(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"plot1d"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"plot2d")
getstructchargedensityplot%plot2d=>null()
Do i=0,len-1
getstructchargedensityplot%plot2d=>getstructplot2d(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"plot2d"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"plot3d")
getstructchargedensityplot%plot3d=>null()
Do i=0,len-1
getstructchargedensityplot%plot3d=>getstructplot3d(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"plot3d"),0)) ) 
enddo

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructexccplot(thisnode)

implicit none
type(Node),pointer::thisnode
type(exccplot_type),pointer::getstructexccplot

integer::len=1,i=0
allocate(getstructexccplot)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at exccplot"
#endif
      
            len= countChildEmentsWithName(thisnode,"plot1d")
getstructexccplot%plot1d=>null()
Do i=0,len-1
getstructexccplot%plot1d=>getstructplot1d(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"plot1d"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"plot2d")
getstructexccplot%plot2d=>null()
Do i=0,len-1
getstructexccplot%plot2d=>getstructplot2d(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"plot2d"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"plot3d")
getstructexccplot%plot3d=>null()
Do i=0,len-1
getstructexccplot%plot3d=>getstructplot3d(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"plot3d"),0)) ) 
enddo

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructelfplot(thisnode)

implicit none
type(Node),pointer::thisnode
type(elfplot_type),pointer::getstructelfplot

integer::len=1,i=0
allocate(getstructelfplot)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at elfplot"
#endif
      
            len= countChildEmentsWithName(thisnode,"plot1d")
getstructelfplot%plot1d=>null()
Do i=0,len-1
getstructelfplot%plot1d=>getstructplot1d(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"plot1d"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"plot2d")
getstructelfplot%plot2d=>null()
Do i=0,len-1
getstructelfplot%plot2d=>getstructplot2d(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"plot2d"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"plot3d")
getstructelfplot%plot3d=>null()
Do i=0,len-1
getstructelfplot%plot3d=>getstructplot3d(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"plot3d"),0)) ) 
enddo

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructmvecfield(thisnode)

implicit none
type(Node),pointer::thisnode
type(mvecfield_type),pointer::getstructmvecfield

integer::len=1,i=0
allocate(getstructmvecfield)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at mvecfield"
#endif
      
            len= countChildEmentsWithName(thisnode,"plot2d")
getstructmvecfield%plot2d=>null()
Do i=0,len-1
getstructmvecfield%plot2d=>getstructplot2d(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"plot2d"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"plot3d")
getstructmvecfield%plot3d=>null()
Do i=0,len-1
getstructmvecfield%plot3d=>getstructplot3d(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"plot3d"),0)) ) 
enddo

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructxcmvecfield(thisnode)

implicit none
type(Node),pointer::thisnode
type(xcmvecfield_type),pointer::getstructxcmvecfield

integer::len=1,i=0
allocate(getstructxcmvecfield)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at xcmvecfield"
#endif
      
            len= countChildEmentsWithName(thisnode,"plot2d")
getstructxcmvecfield%plot2d=>null()
Do i=0,len-1
getstructxcmvecfield%plot2d=>getstructplot2d(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"plot2d"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"plot3d")
getstructxcmvecfield%plot3d=>null()
Do i=0,len-1
getstructxcmvecfield%plot3d=>getstructplot3d(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"plot3d"),0)) ) 
enddo

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructelectricfield(thisnode)

implicit none
type(Node),pointer::thisnode
type(electricfield_type),pointer::getstructelectricfield

integer::len=1,i=0
allocate(getstructelectricfield)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at electricfield"
#endif
      
            len= countChildEmentsWithName(thisnode,"plot2d")
getstructelectricfield%plot2d=>null()
Do i=0,len-1
getstructelectricfield%plot2d=>getstructplot2d(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"plot2d"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"plot3d")
getstructelectricfield%plot3d=>null()
Do i=0,len-1
getstructelectricfield%plot3d=>getstructplot3d(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"plot3d"),0)) ) 
enddo

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructgradmvecfield(thisnode)

implicit none
type(Node),pointer::thisnode
type(gradmvecfield_type),pointer::getstructgradmvecfield

integer::len=1,i=0
allocate(getstructgradmvecfield)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at gradmvecfield"
#endif
      
            len= countChildEmentsWithName(thisnode,"plot1d")
getstructgradmvecfield%plot1d=>null()
Do i=0,len-1
getstructgradmvecfield%plot1d=>getstructplot1d(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"plot1d"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"plot2d")
getstructgradmvecfield%plot2d=>null()
Do i=0,len-1
getstructgradmvecfield%plot2d=>getstructplot2d(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"plot2d"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"plot3d")
getstructgradmvecfield%plot3d=>null()
Do i=0,len-1
getstructgradmvecfield%plot3d=>getstructplot3d(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"plot3d"),0)) ) 
enddo

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructfermisurfaceplot(thisnode)

implicit none
type(Node),pointer::thisnode
type(fermisurfaceplot_type),pointer::getstructfermisurfaceplot
type(Node),pointer::np


integer::len=1,i=0
allocate(getstructfermisurfaceplot)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at fermisurfaceplot"
#endif
      
nullify(np)  
np=>getAttributeNode(thisnode,"nstfsp")
getstructfermisurfaceplot%nstfsp=6
if(associated(np)) then
       call extractDataAttribute(thisnode,"nstfsp",getstructfermisurfaceplot%nstfsp)
       call removeAttribute(thisnode,"nstfsp")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"separate")
if(associated(np)) then
       call extractDataAttribute(thisnode,"separate",getstructfermisurfaceplot%separate)
       call removeAttribute(thisnode,"separate")  
endif

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructEFG(thisnode)

implicit none
type(Node),pointer::thisnode
type(EFG_type),pointer::getstructEFG

integer::len=1,i=0
allocate(getstructEFG)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at EFG"
#endif
      
      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructmossbauer(thisnode)

implicit none
type(Node),pointer::thisnode
type(mossbauer_type),pointer::getstructmossbauer

integer::len=1,i=0
allocate(getstructmossbauer)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at mossbauer"
#endif
      
      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructmomentummatrix(thisnode)

implicit none
type(Node),pointer::thisnode
type(momentummatrix_type),pointer::getstructmomentummatrix
type(Node),pointer::np


integer::len=1,i=0
allocate(getstructmomentummatrix)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at momentummatrix"
#endif
      
nullify(np)  
np=>getAttributeNode(thisnode,"fastpmat")
getstructmomentummatrix%fastpmat= .false.
if(associated(np)) then
       call extractDataAttribute(thisnode,"fastpmat",getstructmomentummatrix%fastpmat)
       call removeAttribute(thisnode,"fastpmat")  
endif

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructdielectric(thisnode)

implicit none
type(Node),pointer::thisnode
type(dielectric_type),pointer::getstructdielectric
type(Node),pointer::np


integer::len=1,i=0
allocate(getstructdielectric)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at dielectric"
#endif
      
nullify(np)  
np=>getAttributeNode(thisnode,"scissor")
getstructdielectric%scissor=0.0d0
if(associated(np)) then
       call extractDataAttribute(thisnode,"scissor",getstructdielectric%scissor)
       call removeAttribute(thisnode,"scissor")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"intraband")
getstructdielectric%intraband= .false.
if(associated(np)) then
       call extractDataAttribute(thisnode,"intraband",getstructdielectric%intraband)
       call removeAttribute(thisnode,"intraband")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"usegdft")
getstructdielectric%usegdft= .false.
if(associated(np)) then
       call extractDataAttribute(thisnode,"usegdft",getstructdielectric%usegdft)
       call removeAttribute(thisnode,"usegdft")  
endif

      len= countChildEmentsWithName (thisnode,"optcomp")           
allocate(getstructdielectric%optcomp(3,len))
Do i=1,len

getstructdielectric%optcomp(:,i)=getvalueofoptcomp(&
      removechild(thisnode,item(getElementsByTagname(thisnode,&
      "optcomp"),0)))
end do

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructmoke(thisnode)

implicit none
type(Node),pointer::thisnode
type(moke_type),pointer::getstructmoke

integer::len=1,i=0
allocate(getstructmoke)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at moke"
#endif
      
      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructexpiqr(thisnode)

implicit none
type(Node),pointer::thisnode
type(expiqr_type),pointer::getstructexpiqr

integer::len=1,i=0
allocate(getstructexpiqr)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at expiqr"
#endif
      
      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructelnes(thisnode)

implicit none
type(Node),pointer::thisnode
type(elnes_type),pointer::getstructelnes
type(Node),pointer::np


integer::len=1,i=0
allocate(getstructelnes)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at elnes"
#endif
      
nullify(np)  
np=>getAttributeNode(thisnode,"vecql")
getstructelnes%vecql=(/0.0d0,0.0d0,0.0d0/)
if(associated(np)) then
       call extractDataAttribute(thisnode,"vecql",getstructelnes%vecql)
       call removeAttribute(thisnode,"vecql")  
endif

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructeliashberg(thisnode)

implicit none
type(Node),pointer::thisnode
type(eliashberg_type),pointer::getstructeliashberg
type(Node),pointer::np


integer::len=1,i=0
allocate(getstructeliashberg)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at eliashberg"
#endif
      
nullify(np)  
np=>getAttributeNode(thisnode,"mustar")
getstructeliashberg%mustar=0.15d0
if(associated(np)) then
       call extractDataAttribute(thisnode,"mustar",getstructeliashberg%mustar)
       call removeAttribute(thisnode,"mustar")  
endif

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructphonons(thisnode)

implicit none
type(Node),pointer::thisnode
type(phonons_type),pointer::getstructphonons
type(Node),pointer::np


integer::len=1,i=0
allocate(getstructphonons)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at phonons"
#endif
      
nullify(np)  
np=>getAttributeNode(thisnode,"do")
getstructphonons%do= "fromscratch"
if(associated(np)) then
       call extractDataAttribute(thisnode,"do",getstructphonons%do)
       call removeAttribute(thisnode,"do")  
endif
getstructphonons%donumber=stringtonumberphononsdo(getstructphonons%do)

nullify(np)  
np=>getAttributeNode(thisnode,"ngridq")
if(associated(np)) then
       call extractDataAttribute(thisnode,"ngridq",getstructphonons%ngridq)
       call removeAttribute(thisnode,"ngridq")  
        else
        write(*,*)"Parser ERROR: The element 'phonons' requires the attribute 'ngridq' to be defined."
        write(*,*)"stopped"
        stop
        
endif

nullify(np)  
np=>getAttributeNode(thisnode,"reduceq")
getstructphonons%reduceq= .true.
if(associated(np)) then
       call extractDataAttribute(thisnode,"reduceq",getstructphonons%reduceq)
       call removeAttribute(thisnode,"reduceq")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"deltaph")
getstructphonons%deltaph=0.03d0
if(associated(np)) then
       call extractDataAttribute(thisnode,"deltaph",getstructphonons%deltaph)
       call removeAttribute(thisnode,"deltaph")  
endif

            len= countChildEmentsWithName(thisnode,"qpointset")
getstructphonons%qpointset=>null()
Do i=0,len-1
getstructphonons%qpointset=>getstructqpointset(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"qpointset"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"phonondos")
getstructphonons%phonondos=>null()
Do i=0,len-1
getstructphonons%phonondos=>getstructphonondos(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"phonondos"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"phonondispplot")
getstructphonons%phonondispplot=>null()
Do i=0,len-1
getstructphonons%phonondispplot=>getstructphonondispplot(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"phonondispplot"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"reformatdynmat")
getstructphonons%reformatdynmat=>null()
Do i=0,len-1
getstructphonons%reformatdynmat=>getstructreformatdynmat(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"reformatdynmat"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"interpolate")
getstructphonons%interpolate=>null()
Do i=0,len-1
getstructphonons%interpolate=>getstructinterpolate(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"interpolate"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"parts")
getstructphonons%parts=>null()
Do i=0,len-1
getstructphonons%parts=>getstructparts(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"parts"),0)) ) 
enddo

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructphonondos(thisnode)

implicit none
type(Node),pointer::thisnode
type(phonondos_type),pointer::getstructphonondos
type(Node),pointer::np


integer::len=1,i=0
allocate(getstructphonondos)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at phonondos"
#endif
      
nullify(np)  
np=>getAttributeNode(thisnode,"nwdos")
getstructphonondos%nwdos=500
if(associated(np)) then
       call extractDataAttribute(thisnode,"nwdos",getstructphonondos%nwdos)
       call removeAttribute(thisnode,"nwdos")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"ngrdos")
getstructphonondos%ngrdos=100
if(associated(np)) then
       call extractDataAttribute(thisnode,"ngrdos",getstructphonondos%ngrdos)
       call removeAttribute(thisnode,"ngrdos")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"nsmdos")
getstructphonondos%nsmdos=0
if(associated(np)) then
       call extractDataAttribute(thisnode,"nsmdos",getstructphonondos%nsmdos)
       call removeAttribute(thisnode,"nsmdos")  
endif

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructphonondispplot(thisnode)

implicit none
type(Node),pointer::thisnode
type(phonondispplot_type),pointer::getstructphonondispplot

integer::len=1,i=0
allocate(getstructphonondispplot)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at phonondispplot"
#endif
      
            len= countChildEmentsWithName(thisnode,"plot1d")
getstructphonondispplot%plot1d=>null()
Do i=0,len-1
getstructphonondispplot%plot1d=>getstructplot1d(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"plot1d"),0)) ) 
enddo

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructreformatdynmat(thisnode)

implicit none
type(Node),pointer::thisnode
type(reformatdynmat_type),pointer::getstructreformatdynmat

integer::len=1,i=0
allocate(getstructreformatdynmat)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at reformatdynmat"
#endif
      getstructreformatdynmat%exists=.false.
      if (associated(thisnode))  getstructreformatdynmat%exists=.true.
      
      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructinterpolate(thisnode)

implicit none
type(Node),pointer::thisnode
type(interpolate_type),pointer::getstructinterpolate
type(Node),pointer::np


integer::len=1,i=0
allocate(getstructinterpolate)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at interpolate"
#endif
      
nullify(np)  
np=>getAttributeNode(thisnode,"ngridq")
if(associated(np)) then
       call extractDataAttribute(thisnode,"ngridq",getstructinterpolate%ngridq)
       call removeAttribute(thisnode,"ngridq")  
        else
        write(*,*)"Parser ERROR: The element 'interpolate' requires the attribute 'ngridq' to be defined."
        write(*,*)"stopped"
        stop
        
endif

nullify(np)  
np=>getAttributeNode(thisnode,"vqloff")
getstructinterpolate%vqloff=(/0.0d0,0.0d0,0.0d0/)
if(associated(np)) then
       call extractDataAttribute(thisnode,"vqloff",getstructinterpolate%vqloff)
       call removeAttribute(thisnode,"vqloff")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"writeeigenvectors")
getstructinterpolate%writeeigenvectors= .false.
if(associated(np)) then
       call extractDataAttribute(thisnode,"writeeigenvectors",getstructinterpolate%writeeigenvectors)
       call removeAttribute(thisnode,"writeeigenvectors")  
endif

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructxs(thisnode)

implicit none
type(Node),pointer::thisnode
type(xs_type),pointer::getstructxs
type(Node),pointer::np


integer::len=1,i=0
allocate(getstructxs)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at xs"
#endif
      
nullify(np)  
np=>getAttributeNode(thisnode,"emattype")
getstructxs%emattype=1
if(associated(np)) then
       call extractDataAttribute(thisnode,"emattype",getstructxs%emattype)
       call removeAttribute(thisnode,"emattype")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"dfoffdiag")
getstructxs%dfoffdiag= .false.
if(associated(np)) then
       call extractDataAttribute(thisnode,"dfoffdiag",getstructxs%dfoffdiag)
       call removeAttribute(thisnode,"dfoffdiag")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"lmaxapwwf")
getstructxs%lmaxapwwf=-1
if(associated(np)) then
       call extractDataAttribute(thisnode,"lmaxapwwf",getstructxs%lmaxapwwf)
       call removeAttribute(thisnode,"lmaxapwwf")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"lmaxemat")
getstructxs%lmaxemat=3
if(associated(np)) then
       call extractDataAttribute(thisnode,"lmaxemat",getstructxs%lmaxemat)
       call removeAttribute(thisnode,"lmaxemat")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"emaxdf")
getstructxs%emaxdf=1.0d10
if(associated(np)) then
       call extractDataAttribute(thisnode,"emaxdf",getstructxs%emaxdf)
       call removeAttribute(thisnode,"emaxdf")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"broad")
getstructxs%broad=0.01d0
if(associated(np)) then
       call extractDataAttribute(thisnode,"broad",getstructxs%broad)
       call removeAttribute(thisnode,"broad")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"epsdfde")
getstructxs%epsdfde=1.0d-8
if(associated(np)) then
       call extractDataAttribute(thisnode,"epsdfde",getstructxs%epsdfde)
       call removeAttribute(thisnode,"epsdfde")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"tevout")
getstructxs%tevout= .false.
if(associated(np)) then
       call extractDataAttribute(thisnode,"tevout",getstructxs%tevout)
       call removeAttribute(thisnode,"tevout")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"xstype")
if(associated(np)) then
       call extractDataAttribute(thisnode,"xstype",getstructxs%xstype)
       call removeAttribute(thisnode,"xstype")  
        else
        write(*,*)"Parser ERROR: The element 'xs' requires the attribute 'xstype' to be defined."
        write(*,*)"stopped"
        stop
        
endif
getstructxs%xstypenumber=stringtonumberxsxstype(getstructxs%xstype)

nullify(np)  
np=>getAttributeNode(thisnode,"fastpmat")
getstructxs%fastpmat= .true.
if(associated(np)) then
       call extractDataAttribute(thisnode,"fastpmat",getstructxs%fastpmat)
       call removeAttribute(thisnode,"fastpmat")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"fastemat")
getstructxs%fastemat= .true.
if(associated(np)) then
       call extractDataAttribute(thisnode,"fastemat",getstructxs%fastemat)
       call removeAttribute(thisnode,"fastemat")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"tappinfo")
getstructxs%tappinfo= .false.
if(associated(np)) then
       call extractDataAttribute(thisnode,"tappinfo",getstructxs%tappinfo)
       call removeAttribute(thisnode,"tappinfo")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"dbglev")
getstructxs%dbglev=0
if(associated(np)) then
       call extractDataAttribute(thisnode,"dbglev",getstructxs%dbglev)
       call removeAttribute(thisnode,"dbglev")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"gqmaxtype")
getstructxs%gqmaxtype= "|G+q|"
if(associated(np)) then
       call extractDataAttribute(thisnode,"gqmaxtype",getstructxs%gqmaxtype)
       call removeAttribute(thisnode,"gqmaxtype")  
endif
getstructxs%gqmaxtypenumber=stringtonumberxsgqmaxtype(getstructxs%gqmaxtype)

nullify(np)  
np=>getAttributeNode(thisnode,"gqmax")
getstructxs%gqmax=0.0d0
if(associated(np)) then
       call extractDataAttribute(thisnode,"gqmax",getstructxs%gqmax)
       call removeAttribute(thisnode,"gqmax")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"nosym")
getstructxs%nosym= .false.
if(associated(np)) then
       call extractDataAttribute(thisnode,"nosym",getstructxs%nosym)
       call removeAttribute(thisnode,"nosym")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"ngridk")
getstructxs%ngridk=(/1,1,1/)
if(associated(np)) then
       call extractDataAttribute(thisnode,"ngridk",getstructxs%ngridk)
       call removeAttribute(thisnode,"ngridk")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"vkloff")
getstructxs%vkloff=(/0.0d0,0.0d0,0.0d0/)
if(associated(np)) then
       call extractDataAttribute(thisnode,"vkloff",getstructxs%vkloff)
       call removeAttribute(thisnode,"vkloff")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"reducek")
getstructxs%reducek= .false.
if(associated(np)) then
       call extractDataAttribute(thisnode,"reducek",getstructxs%reducek)
       call removeAttribute(thisnode,"reducek")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"ngridq")
getstructxs%ngridq=(/1,1,1/)
if(associated(np)) then
       call extractDataAttribute(thisnode,"ngridq",getstructxs%ngridq)
       call removeAttribute(thisnode,"ngridq")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"reduceq")
getstructxs%reduceq= .true.
if(associated(np)) then
       call extractDataAttribute(thisnode,"reduceq",getstructxs%reduceq)
       call removeAttribute(thisnode,"reduceq")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"rgkmax")
getstructxs%rgkmax=7.0d0
if(associated(np)) then
       call extractDataAttribute(thisnode,"rgkmax",getstructxs%rgkmax)
       call removeAttribute(thisnode,"rgkmax")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"swidth")
getstructxs%swidth=0.001d0
if(associated(np)) then
       call extractDataAttribute(thisnode,"swidth",getstructxs%swidth)
       call removeAttribute(thisnode,"swidth")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"lmaxapw")
getstructxs%lmaxapw=10
if(associated(np)) then
       call extractDataAttribute(thisnode,"lmaxapw",getstructxs%lmaxapw)
       call removeAttribute(thisnode,"lmaxapw")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"lmaxmat")
getstructxs%lmaxmat=5
if(associated(np)) then
       call extractDataAttribute(thisnode,"lmaxmat",getstructxs%lmaxmat)
       call removeAttribute(thisnode,"lmaxmat")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"nempty")
getstructxs%nempty=5
if(associated(np)) then
       call extractDataAttribute(thisnode,"nempty",getstructxs%nempty)
       call removeAttribute(thisnode,"nempty")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"scissor")
getstructxs%scissor=0.0d0
if(associated(np)) then
       call extractDataAttribute(thisnode,"scissor",getstructxs%scissor)
       call removeAttribute(thisnode,"scissor")  
endif

            len= countChildEmentsWithName(thisnode,"tddft")
getstructxs%tddft=>null()
Do i=0,len-1
getstructxs%tddft=>getstructtddft(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"tddft"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"screening")
getstructxs%screening=>null()
Do i=0,len-1
getstructxs%screening=>getstructscreening(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"screening"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"BSE")
getstructxs%BSE=>null()
Do i=0,len-1
getstructxs%BSE=>getstructBSE(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"BSE"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"transitions")
getstructxs%transitions=>null()
Do i=0,len-1
getstructxs%transitions=>getstructtransitions(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"transitions"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"qpointset")
getstructxs%qpointset=>null()
Do i=0,len-1
getstructxs%qpointset=>getstructqpointset(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"qpointset"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"tetra")
getstructxs%tetra=>null()
Do i=0,len-1
getstructxs%tetra=>getstructtetra(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"tetra"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"energywindow")
getstructxs%energywindow=>null()
Do i=0,len-1
getstructxs%energywindow=>getstructenergywindow(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"energywindow"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"plan")
getstructxs%plan=>null()
Do i=0,len-1
getstructxs%plan=>getstructplan(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"plan"),0)) ) 
enddo

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructtddft(thisnode)

implicit none
type(Node),pointer::thisnode
type(tddft_type),pointer::getstructtddft
type(Node),pointer::np


integer::len=1,i=0
allocate(getstructtddft)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at tddft"
#endif
      
nullify(np)  
np=>getAttributeNode(thisnode,"intraband")
getstructtddft%intraband= .false.
if(associated(np)) then
       call extractDataAttribute(thisnode,"intraband",getstructtddft%intraband)
       call removeAttribute(thisnode,"intraband")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"torddf")
getstructtddft%torddf= .false.
if(associated(np)) then
       call extractDataAttribute(thisnode,"torddf",getstructtddft%torddf)
       call removeAttribute(thisnode,"torddf")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"tordfxc")
getstructtddft%tordfxc= .false.
if(associated(np)) then
       call extractDataAttribute(thisnode,"tordfxc",getstructtddft%tordfxc)
       call removeAttribute(thisnode,"tordfxc")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"aresdf")
getstructtddft%aresdf= .true.
if(associated(np)) then
       call extractDataAttribute(thisnode,"aresdf",getstructtddft%aresdf)
       call removeAttribute(thisnode,"aresdf")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"aresfxc")
getstructtddft%aresfxc= .true.
if(associated(np)) then
       call extractDataAttribute(thisnode,"aresfxc",getstructtddft%aresfxc)
       call removeAttribute(thisnode,"aresfxc")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"fxcbsesplit")
getstructtddft%fxcbsesplit=1.0d-5
if(associated(np)) then
       call extractDataAttribute(thisnode,"fxcbsesplit",getstructtddft%fxcbsesplit)
       call removeAttribute(thisnode,"fxcbsesplit")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"acont")
getstructtddft%acont= .false.
if(associated(np)) then
       call extractDataAttribute(thisnode,"acont",getstructtddft%acont)
       call removeAttribute(thisnode,"acont")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"nwacont")
getstructtddft%nwacont=0
if(associated(np)) then
       call extractDataAttribute(thisnode,"nwacont",getstructtddft%nwacont)
       call removeAttribute(thisnode,"nwacont")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"lindhard")
getstructtddft%lindhard= .false.
if(associated(np)) then
       call extractDataAttribute(thisnode,"lindhard",getstructtddft%lindhard)
       call removeAttribute(thisnode,"lindhard")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"kerndiag")
getstructtddft%kerndiag= .false.
if(associated(np)) then
       call extractDataAttribute(thisnode,"kerndiag",getstructtddft%kerndiag)
       call removeAttribute(thisnode,"kerndiag")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"lmaxalda")
getstructtddft%lmaxalda=3
if(associated(np)) then
       call extractDataAttribute(thisnode,"lmaxalda",getstructtddft%lmaxalda)
       call removeAttribute(thisnode,"lmaxalda")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"alphalrc")
getstructtddft%alphalrc=0.0d0
if(associated(np)) then
       call extractDataAttribute(thisnode,"alphalrc",getstructtddft%alphalrc)
       call removeAttribute(thisnode,"alphalrc")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"alphalrcdyn")
getstructtddft%alphalrcdyn=0.0d0
if(associated(np)) then
       call extractDataAttribute(thisnode,"alphalrcdyn",getstructtddft%alphalrcdyn)
       call removeAttribute(thisnode,"alphalrcdyn")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"betalrcdyn")
if(associated(np)) then
       call extractDataAttribute(thisnode,"betalrcdyn",getstructtddft%betalrcdyn)
       call removeAttribute(thisnode,"betalrcdyn")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"mdfqtype")
getstructtddft%mdfqtype=0
if(associated(np)) then
       call extractDataAttribute(thisnode,"mdfqtype",getstructtddft%mdfqtype)
       call removeAttribute(thisnode,"mdfqtype")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"fxctype")
getstructtddft%fxctype= "RPA"
if(associated(np)) then
       call extractDataAttribute(thisnode,"fxctype",getstructtddft%fxctype)
       call removeAttribute(thisnode,"fxctype")  
endif
getstructtddft%fxctypenumber=stringtonumbertddftfxctype(getstructtddft%fxctype)

nullify(np)  
np=>getAttributeNode(thisnode,"do")
getstructtddft%do= "fromscratch"
if(associated(np)) then
       call extractDataAttribute(thisnode,"do",getstructtddft%do)
       call removeAttribute(thisnode,"do")  
endif
getstructtddft%donumber=stringtonumbertddftdo(getstructtddft%do)

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructscreening(thisnode)

implicit none
type(Node),pointer::thisnode
type(screening_type),pointer::getstructscreening
type(Node),pointer::np


integer::len=1,i=0
allocate(getstructscreening)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at screening"
#endif
      
nullify(np)  
np=>getAttributeNode(thisnode,"do")
getstructscreening%do= "fromscratch"
if(associated(np)) then
       call extractDataAttribute(thisnode,"do",getstructscreening%do)
       call removeAttribute(thisnode,"do")  
endif
getstructscreening%donumber=stringtonumberscreeningdo(getstructscreening%do)

nullify(np)  
np=>getAttributeNode(thisnode,"nosym")
getstructscreening%nosym= .false.
if(associated(np)) then
       call extractDataAttribute(thisnode,"nosym",getstructscreening%nosym)
       call removeAttribute(thisnode,"nosym")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"ngridk")
getstructscreening%ngridk=(/0,0,0/)
if(associated(np)) then
       call extractDataAttribute(thisnode,"ngridk",getstructscreening%ngridk)
       call removeAttribute(thisnode,"ngridk")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"reducek")
getstructscreening%reducek= .false.
if(associated(np)) then
       call extractDataAttribute(thisnode,"reducek",getstructscreening%reducek)
       call removeAttribute(thisnode,"reducek")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"vkloff")
getstructscreening%vkloff=(/-1.0d0,-1.0d0,-1.0d0/)
if(associated(np)) then
       call extractDataAttribute(thisnode,"vkloff",getstructscreening%vkloff)
       call removeAttribute(thisnode,"vkloff")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"rgkmax")
getstructscreening%rgkmax=0.0d0
if(associated(np)) then
       call extractDataAttribute(thisnode,"rgkmax",getstructscreening%rgkmax)
       call removeAttribute(thisnode,"rgkmax")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"nempty")
getstructscreening%nempty=0
if(associated(np)) then
       call extractDataAttribute(thisnode,"nempty",getstructscreening%nempty)
       call removeAttribute(thisnode,"nempty")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"screentype")
getstructscreening%screentype= "full"
if(associated(np)) then
       call extractDataAttribute(thisnode,"screentype",getstructscreening%screentype)
       call removeAttribute(thisnode,"screentype")  
endif
getstructscreening%screentypenumber=stringtonumberscreeningscreentype(getstructscreening%screentype)

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructBSE(thisnode)

implicit none
type(Node),pointer::thisnode
type(BSE_type),pointer::getstructBSE
type(Node),pointer::np


integer::len=1,i=0
allocate(getstructBSE)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at BSE"
#endif
      
nullify(np)  
np=>getAttributeNode(thisnode,"nosym")
getstructBSE%nosym= .false.
if(associated(np)) then
       call extractDataAttribute(thisnode,"nosym",getstructBSE%nosym)
       call removeAttribute(thisnode,"nosym")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"reducek")
getstructBSE%reducek= .false.
if(associated(np)) then
       call extractDataAttribute(thisnode,"reducek",getstructBSE%reducek)
       call removeAttribute(thisnode,"reducek")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"vkloff")
getstructBSE%vkloff=(/-1.0d0,-1.0d0,-1.0d0/)
if(associated(np)) then
       call extractDataAttribute(thisnode,"vkloff",getstructBSE%vkloff)
       call removeAttribute(thisnode,"vkloff")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"rgkmax")
getstructBSE%rgkmax=0.0d0
if(associated(np)) then
       call extractDataAttribute(thisnode,"rgkmax",getstructBSE%rgkmax)
       call removeAttribute(thisnode,"rgkmax")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"scrherm")
getstructBSE%scrherm=0
if(associated(np)) then
       call extractDataAttribute(thisnode,"scrherm",getstructBSE%scrherm)
       call removeAttribute(thisnode,"scrherm")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"fbzq")
getstructBSE%fbzq= .false.
if(associated(np)) then
       call extractDataAttribute(thisnode,"fbzq",getstructBSE%fbzq)
       call removeAttribute(thisnode,"fbzq")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"sciavtype")
getstructBSE%sciavtype= "spherical"
if(associated(np)) then
       call extractDataAttribute(thisnode,"sciavtype",getstructBSE%sciavtype)
       call removeAttribute(thisnode,"sciavtype")  
endif
getstructBSE%sciavtypenumber=stringtonumberBSEsciavtype(getstructBSE%sciavtype)

nullify(np)  
np=>getAttributeNode(thisnode,"sciavbd")
getstructBSE%sciavbd= .false.
if(associated(np)) then
       call extractDataAttribute(thisnode,"sciavbd",getstructBSE%sciavbd)
       call removeAttribute(thisnode,"sciavbd")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"sciavqhd")
getstructBSE%sciavqhd= .false.
if(associated(np)) then
       call extractDataAttribute(thisnode,"sciavqhd",getstructBSE%sciavqhd)
       call removeAttribute(thisnode,"sciavqhd")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"sciavqwg")
getstructBSE%sciavqwg= .false.
if(associated(np)) then
       call extractDataAttribute(thisnode,"sciavqwg",getstructBSE%sciavqwg)
       call removeAttribute(thisnode,"sciavqwg")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"sciavqbd")
getstructBSE%sciavqbd= .false.
if(associated(np)) then
       call extractDataAttribute(thisnode,"sciavqbd",getstructBSE%sciavqbd)
       call removeAttribute(thisnode,"sciavqbd")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"bsedirsing")
getstructBSE%bsedirsing= .false.
if(associated(np)) then
       call extractDataAttribute(thisnode,"bsedirsing",getstructBSE%bsedirsing)
       call removeAttribute(thisnode,"bsedirsing")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"lmaxdielt")
getstructBSE%lmaxdielt=14
if(associated(np)) then
       call extractDataAttribute(thisnode,"lmaxdielt",getstructBSE%lmaxdielt)
       call removeAttribute(thisnode,"lmaxdielt")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"nleblaik")
getstructBSE%nleblaik=5810
if(associated(np)) then
       call extractDataAttribute(thisnode,"nleblaik",getstructBSE%nleblaik)
       call removeAttribute(thisnode,"nleblaik")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"nexcitmax")
getstructBSE%nexcitmax=100
if(associated(np)) then
       call extractDataAttribute(thisnode,"nexcitmax",getstructBSE%nexcitmax)
       call removeAttribute(thisnode,"nexcitmax")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"nstlbsemat")
getstructBSE%nstlbsemat=(/0,0,0,0/)
if(associated(np)) then
       call extractDataAttribute(thisnode,"nstlbsemat",getstructBSE%nstlbsemat)
       call removeAttribute(thisnode,"nstlbsemat")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"nstlbse")
getstructBSE%nstlbse=(/0,0,0,0/)
if(associated(np)) then
       call extractDataAttribute(thisnode,"nstlbse",getstructBSE%nstlbse)
       call removeAttribute(thisnode,"nstlbse")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"aresbse")
getstructBSE%aresbse= .true.
if(associated(np)) then
       call extractDataAttribute(thisnode,"aresbse",getstructBSE%aresbse)
       call removeAttribute(thisnode,"aresbse")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"bsetype")
getstructBSE%bsetype= "singlet"
if(associated(np)) then
       call extractDataAttribute(thisnode,"bsetype",getstructBSE%bsetype)
       call removeAttribute(thisnode,"bsetype")  
endif
getstructBSE%bsetypenumber=stringtonumberBSEbsetype(getstructBSE%bsetype)

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructtransitions(thisnode)

implicit none
type(Node),pointer::thisnode
type(transitions_type),pointer::getstructtransitions

integer::len=1,i=0
allocate(getstructtransitions)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at transitions"
#endif
      
            len= countChildEmentsWithName(thisnode,"individual")
getstructtransitions%individual=>null()
Do i=0,len-1
getstructtransitions%individual=>getstructindividual(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"individual"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"ranges")
getstructtransitions%ranges=>null()
Do i=0,len-1
getstructtransitions%ranges=>getstructranges(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"ranges"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"lists")
getstructtransitions%lists=>null()
Do i=0,len-1
getstructtransitions%lists=>getstructlists(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"lists"),0)) ) 
enddo

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructindividual(thisnode)

implicit none
type(Node),pointer::thisnode
type(individual_type),pointer::getstructindividual

integer::len=1,i=0
allocate(getstructindividual)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at individual"
#endif
      
            len= countChildEmentsWithName(thisnode,"trans")
     
allocate(getstructindividual%transarray(len))
Do i=0,len-1
getstructindividual%transarray(i+1)%trans=>getstructtrans(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"trans"),0)) ) 
enddo

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructtrans(thisnode)

implicit none
type(Node),pointer::thisnode
type(trans_type),pointer::getstructtrans
type(Node),pointer::np


integer::len=1,i=0
allocate(getstructtrans)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at trans"
#endif
      
nullify(np)  
np=>getAttributeNode(thisnode,"action")
getstructtrans%action= "include"
if(associated(np)) then
       call extractDataAttribute(thisnode,"action",getstructtrans%action)
       call removeAttribute(thisnode,"action")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"kpointnumber")
getstructtrans%kpointnumber=0
if(associated(np)) then
       call extractDataAttribute(thisnode,"kpointnumber",getstructtrans%kpointnumber)
       call removeAttribute(thisnode,"kpointnumber")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"initial")
getstructtrans%initial=0
if(associated(np)) then
       call extractDataAttribute(thisnode,"initial",getstructtrans%initial)
       call removeAttribute(thisnode,"initial")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"final")
getstructtrans%final=0
if(associated(np)) then
       call extractDataAttribute(thisnode,"final",getstructtrans%final)
       call removeAttribute(thisnode,"final")  
endif

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructranges(thisnode)

implicit none
type(Node),pointer::thisnode
type(ranges_type),pointer::getstructranges

integer::len=1,i=0
allocate(getstructranges)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at ranges"
#endif
      
            len= countChildEmentsWithName(thisnode,"range")
     
allocate(getstructranges%rangearray(len))
Do i=0,len-1
getstructranges%rangearray(i+1)%range=>getstructrange(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"range"),0)) ) 
enddo

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructrange(thisnode)

implicit none
type(Node),pointer::thisnode
type(range_type),pointer::getstructrange
type(Node),pointer::np


integer::len=1,i=0
allocate(getstructrange)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at range"
#endif
      
nullify(np)  
np=>getAttributeNode(thisnode,"action")
getstructrange%action= "include"
if(associated(np)) then
       call extractDataAttribute(thisnode,"action",getstructrange%action)
       call removeAttribute(thisnode,"action")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"statestype")
if(associated(np)) then
       call extractDataAttribute(thisnode,"statestype",getstructrange%statestype)
       call removeAttribute(thisnode,"statestype")  
        else
        write(*,*)"Parser ERROR: The element 'range' requires the attribute 'statestype' to be defined."
        write(*,*)"stopped"
        stop
        
endif

nullify(np)  
np=>getAttributeNode(thisnode,"kpointnumber")
getstructrange%kpointnumber=0
if(associated(np)) then
       call extractDataAttribute(thisnode,"kpointnumber",getstructrange%kpointnumber)
       call removeAttribute(thisnode,"kpointnumber")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"start")
getstructrange%start=0
if(associated(np)) then
       call extractDataAttribute(thisnode,"start",getstructrange%start)
       call removeAttribute(thisnode,"start")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"stop")
getstructrange%stop=0
if(associated(np)) then
       call extractDataAttribute(thisnode,"stop",getstructrange%stop)
       call removeAttribute(thisnode,"stop")  
endif

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructlists(thisnode)

implicit none
type(Node),pointer::thisnode
type(lists_type),pointer::getstructlists

integer::len=1,i=0
allocate(getstructlists)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at lists"
#endif
      
            len= countChildEmentsWithName(thisnode,"istate")
     
allocate(getstructlists%istatearray(len))
Do i=0,len-1
getstructlists%istatearray(i+1)%istate=>getstructistate(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"istate"),0)) ) 
enddo

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructistate(thisnode)

implicit none
type(Node),pointer::thisnode
type(istate_type),pointer::getstructistate
type(Node),pointer::np


integer::len=1,i=0
allocate(getstructistate)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at istate"
#endif
      
nullify(np)  
np=>getAttributeNode(thisnode,"action")
getstructistate%action= "include"
if(associated(np)) then
       call extractDataAttribute(thisnode,"action",getstructistate%action)
       call removeAttribute(thisnode,"action")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"statestype")
if(associated(np)) then
       call extractDataAttribute(thisnode,"statestype",getstructistate%statestype)
       call removeAttribute(thisnode,"statestype")  
        else
        write(*,*)"Parser ERROR: The element 'istate' requires the attribute 'statestype' to be defined."
        write(*,*)"stopped"
        stop
        
endif

nullify(np)  
np=>getAttributeNode(thisnode,"kpointnumber")
getstructistate%kpointnumber=0
if(associated(np)) then
       call extractDataAttribute(thisnode,"kpointnumber",getstructistate%kpointnumber)
       call removeAttribute(thisnode,"kpointnumber")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"state")
getstructistate%state=0
if(associated(np)) then
       call extractDataAttribute(thisnode,"state",getstructistate%state)
       call removeAttribute(thisnode,"state")  
endif

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructtetra(thisnode)

implicit none
type(Node),pointer::thisnode
type(tetra_type),pointer::getstructtetra
type(Node),pointer::np


integer::len=1,i=0
allocate(getstructtetra)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at tetra"
#endif
      
nullify(np)  
np=>getAttributeNode(thisnode,"tetraocc")
getstructtetra%tetraocc= .false.
if(associated(np)) then
       call extractDataAttribute(thisnode,"tetraocc",getstructtetra%tetraocc)
       call removeAttribute(thisnode,"tetraocc")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"tetradf")
getstructtetra%tetradf= .false.
if(associated(np)) then
       call extractDataAttribute(thisnode,"tetradf",getstructtetra%tetradf)
       call removeAttribute(thisnode,"tetradf")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"kordexc")
getstructtetra%kordexc= .false.
if(associated(np)) then
       call extractDataAttribute(thisnode,"kordexc",getstructtetra%kordexc)
       call removeAttribute(thisnode,"kordexc")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"cw1k")
getstructtetra%cw1k= .false.
if(associated(np)) then
       call extractDataAttribute(thisnode,"cw1k",getstructtetra%cw1k)
       call removeAttribute(thisnode,"cw1k")  
endif

nullify(np)  
np=>getAttributeNode(thisnode,"qweights")
getstructtetra%qweights=1
if(associated(np)) then
       call extractDataAttribute(thisnode,"qweights",getstructtetra%qweights)
       call removeAttribute(thisnode,"qweights")  
endif

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructplan(thisnode)

implicit none
type(Node),pointer::thisnode
type(plan_type),pointer::getstructplan

integer::len=1,i=0
allocate(getstructplan)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at plan"
#endif
      
            len= countChildEmentsWithName(thisnode,"doonly")
     
allocate(getstructplan%doonlyarray(len))
Do i=0,len-1
getstructplan%doonlyarray(i+1)%doonly=>getstructdoonly(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"doonly"),0)) ) 
enddo

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructdoonly(thisnode)

implicit none
type(Node),pointer::thisnode
type(doonly_type),pointer::getstructdoonly
type(Node),pointer::np


integer::len=1,i=0
allocate(getstructdoonly)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at doonly"
#endif
      
nullify(np)  
np=>getAttributeNode(thisnode,"task")
if(associated(np)) then
       call extractDataAttribute(thisnode,"task",getstructdoonly%task)
       call removeAttribute(thisnode,"task")  
        else
        write(*,*)"Parser ERROR: The element 'doonly' requires the attribute 'task' to be defined."
        write(*,*)"stopped"
        stop
        
endif
getstructdoonly%tasknumber=stringtonumberdoonlytask(getstructdoonly%task)

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructqpointset(thisnode)

implicit none
type(Node),pointer::thisnode
type(qpointset_type),pointer::getstructqpointset

integer::len=1,i=0
allocate(getstructqpointset)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at qpointset"
#endif
      
      len= countChildEmentsWithName (thisnode,"qpoint")           
allocate(getstructqpointset%qpoint(3,len))
Do i=1,len

getstructqpointset%qpoint(:,i)=getvalueofqpoint(&
      removechild(thisnode,item(getElementsByTagname(thisnode,&
      "qpoint"),0)))
end do

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructparts(thisnode)

implicit none
type(Node),pointer::thisnode
type(parts_type),pointer::getstructparts

integer::len=1,i=0
allocate(getstructparts)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at parts"
#endif
      
            len= countChildEmentsWithName(thisnode,"dopart")
     
allocate(getstructparts%dopartarray(len))
Do i=0,len-1
getstructparts%dopartarray(i+1)%dopart=>getstructdopart(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"dopart"),0)) ) 
enddo

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructdopart(thisnode)

implicit none
type(Node),pointer::thisnode
type(dopart_type),pointer::getstructdopart
type(Node),pointer::np


integer::len=1,i=0
allocate(getstructdopart)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at dopart"
#endif
      
nullify(np)  
np=>getAttributeNode(thisnode,"id")
if(associated(np)) then
       call extractDataAttribute(thisnode,"id",getstructdopart%id)
       call removeAttribute(thisnode,"id")  
        else
        write(*,*)"Parser ERROR: The element 'dopart' requires the attribute 'id' to be defined."
        write(*,*)"stopped"
        stop
        
endif

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function
 
function getvalueofpointstatepair(thisnode)
implicit none
type(Node),pointer::thisnode
 integer::getvalueofpointstatepair(2)

#ifdef INPUTDEBUG
  write(*,*)"we are at pointstatepair"
#endif  
   call extractDataContent(thisnode,  getvalueofpointstatepair)
end function 
function getvalueoftitle(thisnode)
implicit none
type(Node),pointer::thisnode
 character(512)::getvalueoftitle

#ifdef INPUTDEBUG
  write(*,*)"we are at title"
#endif  
   call extractDataContent(thisnode,  getvalueoftitle)
end function 
function getvalueofbasevect(thisnode)
implicit none
type(Node),pointer::thisnode
 real(8)::getvalueofbasevect(3)

#ifdef INPUTDEBUG
  write(*,*)"we are at basevect"
#endif  
   call extractDataContent(thisnode,  getvalueofbasevect)
end function 
function getvalueofoptcomp(thisnode)
implicit none
type(Node),pointer::thisnode
 integer::getvalueofoptcomp(3)

#ifdef INPUTDEBUG
  write(*,*)"we are at optcomp"
#endif  
   call extractDataContent(thisnode,  getvalueofoptcomp)
end function 
function getvalueofkeywords(thisnode)
implicit none
type(Node),pointer::thisnode
 character(512)::getvalueofkeywords

#ifdef INPUTDEBUG
  write(*,*)"we are at keywords"
#endif  
   call extractDataContent(thisnode,  getvalueofkeywords)
end function 
function getvalueofqpoint(thisnode)
implicit none
type(Node),pointer::thisnode
 real(8)::getvalueofqpoint(3)

#ifdef INPUTDEBUG
  write(*,*)"we are at qpoint"
#endif  
   call extractDataContent(thisnode,  getvalueofqpoint)
end function
 integer function  stringtonumberdo(string) 
 character(80),intent(in)::string
 select case(trim(adjustl(string)))
case('fromscratch')
 stringtonumberdo=-1
case('fromfile')
 stringtonumberdo=-1
case('skip')
 stringtonumberdo=-1
case('')
 stringtonumberdo=0
case default
write(*,*) "Parser ERROR: '", string,"' is not valid selection fordo "
stop 
end select
end function

 
 integer function  stringtonumberaction(string) 
 character(80),intent(in)::string
 select case(trim(adjustl(string)))
case('include')
 stringtonumberaction=-1
case('exclude')
 stringtonumberaction=-1
case('')
 stringtonumberaction=0
case default
write(*,*) "Parser ERROR: '", string,"' is not valid selection foraction "
stop 
end select
end function

 
 integer function  stringtonumberstatestype(string) 
 character(80),intent(in)::string
 select case(trim(adjustl(string)))
case('initialstates')
 stringtonumberstatestype=-1
case('finalstates')
 stringtonumberstatestype=-1
case('')
 stringtonumberstatestype=0
case default
write(*,*) "Parser ERROR: '", string,"' is not valid selection forstatestype "
stop 
end select
end function

 
 integer function  stringtonumberspinfixspin(string) 
 character(80),intent(in)::string
 select case(trim(adjustl(string)))
case('none')
 stringtonumberspinfixspin=0
case('total FSM')
 stringtonumberspinfixspin=1
case('localmt FSM')
 stringtonumberspinfixspin=2
case('both')
 stringtonumberspinfixspin=3
case('')
 stringtonumberspinfixspin=0
case default
write(*,*) "Parser ERROR: '", string,"' is not valid selection forfixspin "
stop 
end select
end function

 
 integer function  stringtonumbersolvertype(string) 
 character(80),intent(in)::string
 select case(trim(adjustl(string)))
case('Lapack')
 stringtonumbersolvertype=1
case('Arpack')
 stringtonumbersolvertype=2
case('DIIS')
 stringtonumbersolvertype=3
case('')
 stringtonumbersolvertype=0
case default
write(*,*) "Parser ERROR: '", string,"' is not valid selection fortype "
stop 
end select
end function

 
 integer function  stringtonumberoutputstate(string) 
 character(80),intent(in)::string
 select case(trim(adjustl(string)))
case('binary')
 stringtonumberoutputstate=-1
case('XML')
 stringtonumberoutputstate=-1
case('')
 stringtonumberoutputstate=0
case default
write(*,*) "Parser ERROR: '", string,"' is not valid selection forstate "
stop 
end select
end function

 
 integer function  stringtonumberlibxcexchange(string) 
 character(80),intent(in)::string
 select case(trim(adjustl(string)))
case('XC_LDA_X')
 stringtonumberlibxcexchange=1
case('XC_LDA_X_2D')
 stringtonumberlibxcexchange=19
case('XC_GGA_X_PBE')
 stringtonumberlibxcexchange=101
case('XC_GGA_X_PBE_R')
 stringtonumberlibxcexchange=102
case('XC_GGA_X_B86')
 stringtonumberlibxcexchange=103
case('XC_GGA_X_B86_R')
 stringtonumberlibxcexchange=104
case('XC_GGA_X_B86_MGC')
 stringtonumberlibxcexchange=105
case('XC_GGA_X_B88')
 stringtonumberlibxcexchange=106
case('XC_GGA_X_G96')
 stringtonumberlibxcexchange=107
case('XC_GGA_X_PW86')
 stringtonumberlibxcexchange=108
case('XC_GGA_X_PW91')
 stringtonumberlibxcexchange=109
case('XC_GGA_X_OPTX')
 stringtonumberlibxcexchange=110
case('XC_GGA_X_DK87_R1')
 stringtonumberlibxcexchange=111
case('XC_GGA_X_DK87_R2')
 stringtonumberlibxcexchange=112
case('XC_GGA_X_LG93')
 stringtonumberlibxcexchange=113
case('XC_GGA_X_FT97_A')
 stringtonumberlibxcexchange=114
case('XC_GGA_X_FT97_B')
 stringtonumberlibxcexchange=115
case('XC_GGA_X_PBE_SOL')
 stringtonumberlibxcexchange=116
case('XC_GGA_X_RPBE')
 stringtonumberlibxcexchange=117
case('XC_GGA_X_WC')
 stringtonumberlibxcexchange=118
case('XC_GGA_X_mPW91')
 stringtonumberlibxcexchange=119
case('XC_GGA_X_AM05')
 stringtonumberlibxcexchange=120
case('XC_GGA_X_PBEA')
 stringtonumberlibxcexchange=121
case('XC_GGA_X_MPBE')
 stringtonumberlibxcexchange=122
case('XC_GGA_X_XPBE')
 stringtonumberlibxcexchange=123
case('XC_GGA_X_2D_B86_MGC')
 stringtonumberlibxcexchange=124
case('XC_GGA_X_BAYESIAN')
 stringtonumberlibxcexchange=125
case('XC_GGA_X_PBE_JSJR')
 stringtonumberlibxcexchange=126
case('')
 stringtonumberlibxcexchange=0
case default
write(*,*) "Parser ERROR: '", string,"' is not valid selection forexchange "
stop 
end select
end function

 
 integer function  stringtonumberlibxccorrelation(string) 
 character(80),intent(in)::string
 select case(trim(adjustl(string)))
case('XC_LDA_C_WIGNER')
 stringtonumberlibxccorrelation=2
case('XC_LDA_C_RPA')
 stringtonumberlibxccorrelation=3
case('XC_LDA_C_HL')
 stringtonumberlibxccorrelation=4
case('XC_LDA_C_GL')
 stringtonumberlibxccorrelation=5
case('XC_LDA_C_XALPHA')
 stringtonumberlibxccorrelation=6
case('XC_LDA_C_VWN')
 stringtonumberlibxccorrelation=7
case('XC_LDA_C_VWN_RPA')
 stringtonumberlibxccorrelation=8
case('XC_LDA_C_PZ')
 stringtonumberlibxccorrelation=9
case('XC_LDA_C_PZ_MOD')
 stringtonumberlibxccorrelation=10
case('XC_LDA_C_OB_PZ')
 stringtonumberlibxccorrelation=11
case('XC_LDA_C_PW')
 stringtonumberlibxccorrelation=12
case('XC_LDA_C_PW_MOD')
 stringtonumberlibxccorrelation=13
case('XC_LDA_C_OB_PW')
 stringtonumberlibxccorrelation=14
case('XC_LDA_C_2D_AMGB')
 stringtonumberlibxccorrelation=15
case('XC_LDA_C_2D_PRM')
 stringtonumberlibxccorrelation=16
case('XC_LDA_C_vBH')
 stringtonumberlibxccorrelation=17
case('XC_LDA_C_1D_CSC')
 stringtonumberlibxccorrelation=18
case('XC_GGA_C_PBE')
 stringtonumberlibxccorrelation=130
case('XC_GGA_C_LYP')
 stringtonumberlibxccorrelation=131
case('XC_GGA_C_P86')
 stringtonumberlibxccorrelation=132
case('XC_GGA_C_PBE_SOL')
 stringtonumberlibxccorrelation=133
case('XC_GGA_C_PW91')
 stringtonumberlibxccorrelation=134
case('XC_GGA_C_AM05')
 stringtonumberlibxccorrelation=135
case('XC_GGA_C_XPBE')
 stringtonumberlibxccorrelation=136
case('XC_GGA_C_LM')
 stringtonumberlibxccorrelation=137
case('XC_GGA_C_PBE_JRGX')
 stringtonumberlibxccorrelation=138
case('')
 stringtonumberlibxccorrelation=0
case default
write(*,*) "Parser ERROR: '", string,"' is not valid selection forcorrelation "
stop 
end select
end function

 
 integer function  stringtonumberlibxcxc(string) 
 character(80),intent(in)::string
 select case(trim(adjustl(string)))
case('none')
 stringtonumberlibxcxc=0
case('XC_GGA_XC_LB')
 stringtonumberlibxcxc=160
case('XC_GGA_XC_HCTH_93')
 stringtonumberlibxcxc=161
case('XC_GGA_XC_HCTH_120')
 stringtonumberlibxcxc=162
case('XC_GGA_XC_HCTH_147')
 stringtonumberlibxcxc=163
case('XC_GGA_XC_HCTH_407')
 stringtonumberlibxcxc=164
case('XC_GGA_XC_EDF1')
 stringtonumberlibxcxc=165
case('XC_GGA_XC_XLYP')
 stringtonumberlibxcxc=166
case('XC_GGA_XC_B97')
 stringtonumberlibxcxc=167
case('XC_GGA_XC_B97_1')
 stringtonumberlibxcxc=168
case('XC_GGA_XC_B97_2')
 stringtonumberlibxcxc=169
case('XC_GGA_XC_B97_D')
 stringtonumberlibxcxc=170
case('XC_GGA_XC_B97_K')
 stringtonumberlibxcxc=171
case('XC_GGA_XC_B97_3')
 stringtonumberlibxcxc=172
case('XC_GGA_XC_PBE1W')
 stringtonumberlibxcxc=173
case('XC_GGA_XC_MPWLYP1W')
 stringtonumberlibxcxc=174
case('XC_GGA_XC_PBELYP1W')
 stringtonumberlibxcxc=175
case('XC_GGA_XC_SB98_1a')
 stringtonumberlibxcxc=176
case('XC_GGA_XC_SB98_1b')
 stringtonumberlibxcxc=177
case('XC_GGA_XC_SB98_1c')
 stringtonumberlibxcxc=178
case('XC_GGA_XC_SB98_2a')
 stringtonumberlibxcxc=179
case('XC_GGA_XC_SB98_2b')
 stringtonumberlibxcxc=180
case('XC_GGA_XC_SB98_2c')
 stringtonumberlibxcxc=181
case('XC_HYB_GGA_XC_B3PW91')
 stringtonumberlibxcxc=401
case('XC_HYB_GGA_XC_B3LYP')
 stringtonumberlibxcxc=402
case('XC_HYB_GGA_XC_B3P86')
 stringtonumberlibxcxc=403
case('XC_HYB_GGA_XC_O3LYP')
 stringtonumberlibxcxc=404
case('XC_HYB_GGA_XC_mPW1K')
 stringtonumberlibxcxc=405
case('XC_HYB_GGA_XC_PBEH')
 stringtonumberlibxcxc=406
case('XC_HYB_GGA_XC_B97')
 stringtonumberlibxcxc=407
case('XC_HYB_GGA_XC_B97_1')
 stringtonumberlibxcxc=408
case('XC_HYB_GGA_XC_B97_2')
 stringtonumberlibxcxc=410
case('XC_HYB_GGA_XC_X3LYP')
 stringtonumberlibxcxc=411
case('XC_HYB_GGA_XC_B1WC')
 stringtonumberlibxcxc=412
case('XC_HYB_GGA_XC_B97_K')
 stringtonumberlibxcxc=413
case('XC_HYB_GGA_XC_B97_3')
 stringtonumberlibxcxc=414
case('XC_HYB_GGA_XC_mPW3PW')
 stringtonumberlibxcxc=415
case('XC_HYB_GGA_XC_B1LYP')
 stringtonumberlibxcxc=416
case('XC_HYB_GGA_XC_B1PW91')
 stringtonumberlibxcxc=417
case('XC_HYB_GGA_XC_mPW1PW')
 stringtonumberlibxcxc=418
case('XC_HYB_GGA_XC_mPW3LYP')
 stringtonumberlibxcxc=419
case('XC_HYB_GGA_XC_SB98_1a')
 stringtonumberlibxcxc=420
case('XC_HYB_GGA_XC_SB98_1b')
 stringtonumberlibxcxc=421
case('XC_HYB_GGA_XC_SB98_1c')
 stringtonumberlibxcxc=422
case('XC_HYB_GGA_XC_SB98_2a')
 stringtonumberlibxcxc=423
case('XC_HYB_GGA_XC_SB98_2b')
 stringtonumberlibxcxc=424
case('XC_HYB_GGA_XC_SB98_2c')
 stringtonumberlibxcxc=425
case('')
 stringtonumberlibxcxc=0
case default
write(*,*) "Parser ERROR: '", string,"' is not valid selection forxc "
stop 
end select
end function

 
 integer function  stringtonumbergroundstatestype(string) 
 character(80),intent(in)::string
 select case(trim(adjustl(string)))
case('Gaussian')
 stringtonumbergroundstatestype=0
case('Methfessel-Paxton 1')
 stringtonumbergroundstatestype=1
case('Methfessel-Paxton 2')
 stringtonumbergroundstatestype=2
case('Fermi Dirac')
 stringtonumbergroundstatestype=3
case('Square-wave impulse')
 stringtonumbergroundstatestype=4
case('')
 stringtonumbergroundstatestype=0
case default
write(*,*) "Parser ERROR: '", string,"' is not valid selection forstype "
stop 
end select
end function

 
 integer function  stringtonumbergroundstatefindlinentype(string) 
 character(80),intent(in)::string
 select case(trim(adjustl(string)))
case('simple')
 stringtonumbergroundstatefindlinentype=-1
case('advanced')
 stringtonumbergroundstatefindlinentype=-1
case('')
 stringtonumbergroundstatefindlinentype=0
case default
write(*,*) "Parser ERROR: '", string,"' is not valid selection forfindlinentype "
stop 
end select
end function

 
 integer function  stringtonumbergroundstatemixer(string) 
 character(80),intent(in)::string
 select case(trim(adjustl(string)))
case('lin')
 stringtonumbergroundstatemixer=1
case('msec')
 stringtonumbergroundstatemixer=2
case('pulay')
 stringtonumbergroundstatemixer=3
case('')
 stringtonumbergroundstatemixer=0
case default
write(*,*) "Parser ERROR: '", string,"' is not valid selection formixer "
stop 
end select
end function

 
 integer function  stringtonumbergroundstatexctype(string) 
 character(80),intent(in)::string
 select case(trim(adjustl(string)))
case('LDAPerdew-Zunger')
 stringtonumbergroundstatexctype=2
case('LSDAPerdew-Wang')
 stringtonumbergroundstatexctype=3
case('LDA-X-alpha')
 stringtonumbergroundstatexctype=4
case('LSDA-Barth-Hedin')
 stringtonumbergroundstatexctype=5
case('GGAPerdew-Burke-Ernzerhof')
 stringtonumbergroundstatexctype=20
case('GGArevPBE')
 stringtonumbergroundstatexctype=21
case('GGAPBEsol')
 stringtonumbergroundstatexctype=22
case('GGA-Wu-Cohen')
 stringtonumbergroundstatexctype=26
case('GGAArmiento-Mattsson')
 stringtonumbergroundstatexctype=30
case('EXX')
 stringtonumbergroundstatexctype=-2
case('none')
 stringtonumbergroundstatexctype=1
case('')
 stringtonumbergroundstatexctype=0
case default
write(*,*) "Parser ERROR: '", string,"' is not valid selection forxctype "
stop 
end select
end function

 
 integer function  stringtonumbergroundstateldapu(string) 
 character(80),intent(in)::string
 select case(trim(adjustl(string)))
case('none')
 stringtonumbergroundstateldapu=0
case('FullyLocalisedLimit')
 stringtonumbergroundstateldapu=1
case('AroundMeanField')
 stringtonumbergroundstateldapu=2
case('FFL-AMF-interpolation')
 stringtonumbergroundstateldapu=3
case('')
 stringtonumbergroundstateldapu=0
case default
write(*,*) "Parser ERROR: '", string,"' is not valid selection forldapu "
stop 
end select
end function

 
 integer function  stringtonumberphononsdo(string) 
 character(80),intent(in)::string
 select case(trim(adjustl(string)))
case('fromscratch')
 stringtonumberphononsdo=-1
case('skip')
 stringtonumberphononsdo=-1
case('')
 stringtonumberphononsdo=0
case default
write(*,*) "Parser ERROR: '", string,"' is not valid selection fordo "
stop 
end select
end function

 
 integer function  stringtonumbertddftfxctype(string) 
 character(80),intent(in)::string
 select case(trim(adjustl(string)))
case('RPA')
 stringtonumbertddftfxctype=0
case('LRCstatic_NLF')
 stringtonumbertddftfxctype=1
case('LRCstatic')
 stringtonumbertddftfxctype=2
case('LRCdyn_NLF')
 stringtonumbertddftfxctype=3
case('LRCdyn')
 stringtonumbertddftfxctype=4
case('ALDA')
 stringtonumbertddftfxctype=5
case('MB1_NLF')
 stringtonumbertddftfxctype=7
case('MB1')
 stringtonumbertddftfxctype=8
case('')
 stringtonumbertddftfxctype=0
case default
write(*,*) "Parser ERROR: '", string,"' is not valid selection forfxctype "
stop 
end select
end function

 
 integer function  stringtonumbertddftdo(string) 
 character(80),intent(in)::string
 select case(trim(adjustl(string)))
case('fromscratch')
 stringtonumbertddftdo=-1
case('fromkernel')
 stringtonumbertddftdo=-1
case('')
 stringtonumbertddftdo=0
case default
write(*,*) "Parser ERROR: '", string,"' is not valid selection fordo "
stop 
end select
end function

 
 integer function  stringtonumberscreeningdo(string) 
 character(80),intent(in)::string
 select case(trim(adjustl(string)))
case('fromscratch')
 stringtonumberscreeningdo=-1
case('skip')
 stringtonumberscreeningdo=-1
case('')
 stringtonumberscreeningdo=0
case default
write(*,*) "Parser ERROR: '", string,"' is not valid selection fordo "
stop 
end select
end function

 
 integer function  stringtonumberscreeningscreentype(string) 
 character(80),intent(in)::string
 select case(trim(adjustl(string)))
case('full')
 stringtonumberscreeningscreentype=-1
case('diag')
 stringtonumberscreeningscreentype=-1
case('noinvdiag')
 stringtonumberscreeningscreentype=-1
case('longrange')
 stringtonumberscreeningscreentype=-1
case('')
 stringtonumberscreeningscreentype=0
case default
write(*,*) "Parser ERROR: '", string,"' is not valid selection forscreentype "
stop 
end select
end function

 
 integer function  stringtonumberBSEsciavtype(string) 
 character(80),intent(in)::string
 select case(trim(adjustl(string)))
case('spherical')
 stringtonumberBSEsciavtype=-1
case('screendiag')
 stringtonumberBSEsciavtype=-1
case('invscreendiag')
 stringtonumberBSEsciavtype=-1
case('')
 stringtonumberBSEsciavtype=0
case default
write(*,*) "Parser ERROR: '", string,"' is not valid selection forsciavtype "
stop 
end select
end function

 
 integer function  stringtonumberBSEbsetype(string) 
 character(80),intent(in)::string
 select case(trim(adjustl(string)))
case('IP')
 stringtonumberBSEbsetype=-1
case('RPA')
 stringtonumberBSEbsetype=-1
case('singlet')
 stringtonumberBSEbsetype=-1
case('triplet')
 stringtonumberBSEbsetype=-1
case('')
 stringtonumberBSEbsetype=0
case default
write(*,*) "Parser ERROR: '", string,"' is not valid selection forbsetype "
stop 
end select
end function

 
 integer function  stringtonumberdoonlytask(string) 
 character(80),intent(in)::string
 select case(trim(adjustl(string)))
case('xsgeneigvec')
 stringtonumberdoonlytask=301
case('tetcalccw')
 stringtonumberdoonlytask=310
case('writepmatxs')
 stringtonumberdoonlytask=320
case('writeemat')
 stringtonumberdoonlytask=330
case('df')
 stringtonumberdoonlytask=340
case('df2')
 stringtonumberdoonlytask=345
case('idf')
 stringtonumberdoonlytask=350
case('scrgeneigvec')
 stringtonumberdoonlytask=401
case('scrtetcalccw')
 stringtonumberdoonlytask=410
case('scrwritepmat')
 stringtonumberdoonlytask=420
case('screen')
 stringtonumberdoonlytask=430
case('scrcoulint')
 stringtonumberdoonlytask=440
case('exccoulint')
 stringtonumberdoonlytask=441
case('bse')
 stringtonumberdoonlytask=445
case('kernxc_bse')
 stringtonumberdoonlytask=450
case('writebandgapgrid')
 stringtonumberdoonlytask=23
case('writepmat')
 stringtonumberdoonlytask=120
case('dielectric')
 stringtonumberdoonlytask=121
case('writepmatasc')
 stringtonumberdoonlytask=321
case('pmatxs2orig')
 stringtonumberdoonlytask=322
case('writeematasc')
 stringtonumberdoonlytask=331
case('writepwmat')
 stringtonumberdoonlytask=335
case('emattest')
 stringtonumberdoonlytask=339
case('x0toasc')
 stringtonumberdoonlytask=341
case('x0tobin')
 stringtonumberdoonlytask=342
case('fxc_alda_check')
 stringtonumberdoonlytask=398
case('kernxc_bse3')
 stringtonumberdoonlytask=451
case('testxs')
 stringtonumberdoonlytask=499
case('xsestimate')
 stringtonumberdoonlytask=700
case('xstiming')
 stringtonumberdoonlytask=701
case('testmain')
 stringtonumberdoonlytask=999
case('portstate(1)')
 stringtonumberdoonlytask=900
case('portstate(2)')
 stringtonumberdoonlytask=901
case('portstate(-1)')
 stringtonumberdoonlytask=910
case('portstate(-2)')
 stringtonumberdoonlytask=911
case('')
 stringtonumberdoonlytask=0
case default
write(*,*) "Parser ERROR: '", string,"' is not valid selection fortask "
stop 
end select
end function

 
 integer function  stringtonumberxsxstype(string) 
 character(80),intent(in)::string
 select case(trim(adjustl(string)))
case('TDDFT')
 stringtonumberxsxstype=-1
case('BSE')
 stringtonumberxsxstype=-1
case('')
 stringtonumberxsxstype=0
case default
write(*,*) "Parser ERROR: '", string,"' is not valid selection forxstype "
stop 
end select
end function

 
 integer function  stringtonumberxsgqmaxtype(string) 
 character(80),intent(in)::string
 select case(trim(adjustl(string)))
case('|G+q|')
 stringtonumberxsgqmaxtype=-1
case('|G|')
 stringtonumberxsgqmaxtype=-1
case('')
 stringtonumberxsgqmaxtype=0
case default
write(*,*) "Parser ERROR: '", string,"' is not valid selection forgqmaxtype "
stop 
end select
end function

 

function countChildEmentsWithName(nodep,name)
  implicit none
  integer::countChildEmentsWithName
  type(Node),pointer ::nodep
  character(len=*),intent(in)::name
  type(NodeList),pointer::children
  type(Node),pointer::child
  
  integer::i
  children=>getChildNodes(nodep)
  countChildEmentsWithName=0
  do i=0,getlength(children)-1
    child=>item(children,i)
    if(name.eq.getNodeName(child)) countChildEmentsWithName=countChildEmentsWithName+1
  end do

end function  
    
! these are some transient helper functions to simplify the port (should not be used)
function isspinorb()
logical::isspinorb
isspinorb=.false.
if(associated(input%groundstate%spin))then
if (input%groundstate%spin%spinorb) then
isspinorb=.true.
endif
endif
end function
function isspinspiral() 
logical::isspinspiral
isspinspiral=.false.
if(associated(input%groundstate%spin))then
if (input%groundstate%spin%spinsprl) then
isspinspiral=.true.
endif
endif
end function

function getfixspinnumber()
implicit none
integer::getfixspinnumber
getfixspinnumber=0
if(associated(input%groundstate%spin))then
getfixspinnumber=input%groundstate%spin%fixspinnumber
endif
end function

function istetraocc()
  implicit none
  logical ::istetraocc
  istetraocc =.false.
  if(associated(input%xs)) then
    if(associated(input%xs%tetra)) then
      istetraocc=input%xs%tetra%tetraocc
    endif
  endif
end function

end module

