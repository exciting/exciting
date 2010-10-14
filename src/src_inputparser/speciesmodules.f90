
module modsp
use inputdom
implicit none
type wf_type
 real(8)::trialEnergy
 integer::matchingOrder
 logical::searchE
end type

type wf_type_array
type(wf_type),pointer::wf
 end type
    type sp_type
 character(1024)::chemicalSymbol
 real(8)::z
 real(8)::mass
 character(512)::name
  type(muffinTin_type),pointer::muffinTin
  type(atomicState_type_array),pointer::atomicStatearray(:)
  type(basis_type),pointer::basis
  type(lorb_type_array),pointer::lorbarray(:)
end type

type sp_type_array
type(sp_type),pointer::sp
 end type
    type muffinTin_type
 real(8)::rmin
 real(8)::rinf
 real(8)::radius
 integer::radialmeshPoints
end type
type atomicState_type
 integer::n
 integer::l
 integer::kappa
 real(8)::occ
 logical::core
end type

type atomicState_type_array
type(atomicState_type),pointer::atomicState
 end type
    type basis_type
 integer::order
  type(wf_type_array),pointer::wfarray(:)
  type(exception_type_array),pointer::exceptionarray(:)
end type
type exception_type
 integer::l
  type(wf_type_array),pointer::wfarray(:)
end type

type exception_type_array
type(exception_type),pointer::exception
 end type
    type lorb_type
 integer::l
  type(wf_type_array),pointer::wfarray(:)
end type

type lorb_type_array
type(lorb_type),pointer::lorb
 end type
    type spdb_type
  type(sp_type_array),pointer::sparray(:)
end type

   type(sp_type)::sp
contains

function getstructwf(thisnode)

implicit none
type(Node),pointer::thisnode
type(wf_type),pointer::getstructwf
type(Node),pointer::np


integer::len=1,i=0
allocate(getstructwf)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at wf"
#endif
      
nullify(np)  
np=>getAttributeNode(thisnode,"trialEnergy")
if(associated(np)) then
       call extractDataAttribute(thisnode,"trialEnergy",getstructwf%trialEnergy)
       call removeAttribute(thisnode,"trialEnergy")  
        else
        write(*,*)"Parser ERROR: The element 'wf' requires the attribute 'trialEnergy' to be defined."
        write(*,*)"stopped"
        stop
        
endif

nullify(np)  
np=>getAttributeNode(thisnode,"matchingOrder")
if(associated(np)) then
       call extractDataAttribute(thisnode,"matchingOrder",getstructwf%matchingOrder)
       call removeAttribute(thisnode,"matchingOrder")  
        else
        write(*,*)"Parser ERROR: The element 'wf' requires the attribute 'matchingOrder' to be defined."
        write(*,*)"stopped"
        stop
        
endif

nullify(np)  
np=>getAttributeNode(thisnode,"searchE")
if(associated(np)) then
       call extractDataAttribute(thisnode,"searchE",getstructwf%searchE)
       call removeAttribute(thisnode,"searchE")  
        else
        write(*,*)"Parser ERROR: The element 'wf' requires the attribute 'searchE' to be defined."
        write(*,*)"stopped"
        stop
        
endif

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructsp(thisnode)

implicit none
type(Node),pointer::thisnode
type(sp_type),pointer::getstructsp
type(Node),pointer::np


integer::len=1,i=0
allocate(getstructsp)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at sp"
#endif
      
nullify(np)  
np=>getAttributeNode(thisnode,"chemicalSymbol")
if(associated(np)) then
       call extractDataAttribute(thisnode,"chemicalSymbol",getstructsp%chemicalSymbol)
       call removeAttribute(thisnode,"chemicalSymbol")  
        else
        write(*,*)"Parser ERROR: The element 'sp' requires the attribute 'chemicalSymbol' to be defined."
        write(*,*)"stopped"
        stop
        
endif

nullify(np)  
np=>getAttributeNode(thisnode,"z")
if(associated(np)) then
       call extractDataAttribute(thisnode,"z",getstructsp%z)
       call removeAttribute(thisnode,"z")  
        else
        write(*,*)"Parser ERROR: The element 'sp' requires the attribute 'z' to be defined."
        write(*,*)"stopped"
        stop
        
endif

nullify(np)  
np=>getAttributeNode(thisnode,"mass")
if(associated(np)) then
       call extractDataAttribute(thisnode,"mass",getstructsp%mass)
       call removeAttribute(thisnode,"mass")  
        else
        write(*,*)"Parser ERROR: The element 'sp' requires the attribute 'mass' to be defined."
        write(*,*)"stopped"
        stop
        
endif

nullify(np)  
np=>getAttributeNode(thisnode,"name")
if(associated(np)) then
       call extractDataAttribute(thisnode,"name",getstructsp%name)
       call removeAttribute(thisnode,"name")  
endif

            len= countChildEmentsWithName(thisnode,"muffinTin")
getstructsp%muffinTin=>null()
Do i=0,len-1
getstructsp%muffinTin=>getstructmuffinTin(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"muffinTin"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"atomicState")
     
allocate(getstructsp%atomicStatearray(len))
Do i=0,len-1
getstructsp%atomicStatearray(i+1)%atomicState=>getstructatomicState(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"atomicState"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"basis")
getstructsp%basis=>null()
Do i=0,len-1
getstructsp%basis=>getstructbasis(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"basis"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"lorb")
     
allocate(getstructsp%lorbarray(len))
Do i=0,len-1
getstructsp%lorbarray(i+1)%lorb=>getstructlorb(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"lorb"),0)) ) 
enddo

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructmuffinTin(thisnode)

implicit none
type(Node),pointer::thisnode
type(muffinTin_type),pointer::getstructmuffinTin
type(Node),pointer::np


integer::len=1,i=0
allocate(getstructmuffinTin)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at muffinTin"
#endif
      
nullify(np)  
np=>getAttributeNode(thisnode,"rmin")
if(associated(np)) then
       call extractDataAttribute(thisnode,"rmin",getstructmuffinTin%rmin)
       call removeAttribute(thisnode,"rmin")  
        else
        write(*,*)"Parser ERROR: The element 'muffinTin' requires the attribute 'rmin' to be defined."
        write(*,*)"stopped"
        stop
        
endif

nullify(np)  
np=>getAttributeNode(thisnode,"rinf")
if(associated(np)) then
       call extractDataAttribute(thisnode,"rinf",getstructmuffinTin%rinf)
       call removeAttribute(thisnode,"rinf")  
        else
        write(*,*)"Parser ERROR: The element 'muffinTin' requires the attribute 'rinf' to be defined."
        write(*,*)"stopped"
        stop
        
endif

nullify(np)  
np=>getAttributeNode(thisnode,"radius")
if(associated(np)) then
       call extractDataAttribute(thisnode,"radius",getstructmuffinTin%radius)
       call removeAttribute(thisnode,"radius")  
        else
        write(*,*)"Parser ERROR: The element 'muffinTin' requires the attribute 'radius' to be defined."
        write(*,*)"stopped"
        stop
        
endif

nullify(np)  
np=>getAttributeNode(thisnode,"radialmeshPoints")
if(associated(np)) then
       call extractDataAttribute(thisnode,"radialmeshPoints",getstructmuffinTin%radialmeshPoints)
       call removeAttribute(thisnode,"radialmeshPoints")  
        else
        write(*,*)"Parser ERROR: The element 'muffinTin' requires the attribute 'radialmeshPoints' to be defined."
        write(*,*)"stopped"
        stop
        
endif

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructatomicState(thisnode)

implicit none
type(Node),pointer::thisnode
type(atomicState_type),pointer::getstructatomicState
type(Node),pointer::np


integer::len=1,i=0
allocate(getstructatomicState)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at atomicState"
#endif
      
nullify(np)  
np=>getAttributeNode(thisnode,"n")
if(associated(np)) then
       call extractDataAttribute(thisnode,"n",getstructatomicState%n)
       call removeAttribute(thisnode,"n")  
        else
        write(*,*)"Parser ERROR: The element 'atomicState' requires the attribute 'n' to be defined."
        write(*,*)"stopped"
        stop
        
endif

nullify(np)  
np=>getAttributeNode(thisnode,"l")
if(associated(np)) then
       call extractDataAttribute(thisnode,"l",getstructatomicState%l)
       call removeAttribute(thisnode,"l")  
        else
        write(*,*)"Parser ERROR: The element 'atomicState' requires the attribute 'l' to be defined."
        write(*,*)"stopped"
        stop
        
endif

nullify(np)  
np=>getAttributeNode(thisnode,"kappa")
if(associated(np)) then
       call extractDataAttribute(thisnode,"kappa",getstructatomicState%kappa)
       call removeAttribute(thisnode,"kappa")  
        else
        write(*,*)"Parser ERROR: The element 'atomicState' requires the attribute 'kappa' to be defined."
        write(*,*)"stopped"
        stop
        
endif

nullify(np)  
np=>getAttributeNode(thisnode,"occ")
if(associated(np)) then
       call extractDataAttribute(thisnode,"occ",getstructatomicState%occ)
       call removeAttribute(thisnode,"occ")  
        else
        write(*,*)"Parser ERROR: The element 'atomicState' requires the attribute 'occ' to be defined."
        write(*,*)"stopped"
        stop
        
endif

nullify(np)  
np=>getAttributeNode(thisnode,"core")
if(associated(np)) then
       call extractDataAttribute(thisnode,"core",getstructatomicState%core)
       call removeAttribute(thisnode,"core")  
        else
        write(*,*)"Parser ERROR: The element 'atomicState' requires the attribute 'core' to be defined."
        write(*,*)"stopped"
        stop
        
endif

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructbasis(thisnode)

implicit none
type(Node),pointer::thisnode
type(basis_type),pointer::getstructbasis
type(Node),pointer::np


integer::len=1,i=0
allocate(getstructbasis)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at basis"
#endif
      
nullify(np)  
np=>getAttributeNode(thisnode,"order")
if(associated(np)) then
       call extractDataAttribute(thisnode,"order",getstructbasis%order)
       call removeAttribute(thisnode,"order")  
        else
        write(*,*)"Parser ERROR: The element 'basis' requires the attribute 'order' to be defined."
        write(*,*)"stopped"
        stop
        
endif

            len= countChildEmentsWithName(thisnode,"wf")
     
allocate(getstructbasis%wfarray(len))
Do i=0,len-1
getstructbasis%wfarray(i+1)%wf=>getstructwf(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"wf"),0)) ) 
enddo

            len= countChildEmentsWithName(thisnode,"exception")
     
allocate(getstructbasis%exceptionarray(len))
Do i=0,len-1
getstructbasis%exceptionarray(i+1)%exception=>getstructexception(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"exception"),0)) ) 
enddo

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructexception(thisnode)

implicit none
type(Node),pointer::thisnode
type(exception_type),pointer::getstructexception
type(Node),pointer::np


integer::len=1,i=0
allocate(getstructexception)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at exception"
#endif
      
nullify(np)  
np=>getAttributeNode(thisnode,"l")
if(associated(np)) then
       call extractDataAttribute(thisnode,"l",getstructexception%l)
       call removeAttribute(thisnode,"l")  
endif

            len= countChildEmentsWithName(thisnode,"wf")
     
allocate(getstructexception%wfarray(len))
Do i=0,len-1
getstructexception%wfarray(i+1)%wf=>getstructwf(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"wf"),0)) ) 
enddo

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructlorb(thisnode)

implicit none
type(Node),pointer::thisnode
type(lorb_type),pointer::getstructlorb
type(Node),pointer::np


integer::len=1,i=0
allocate(getstructlorb)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at lorb"
#endif
      
nullify(np)  
np=>getAttributeNode(thisnode,"l")
if(associated(np)) then
       call extractDataAttribute(thisnode,"l",getstructlorb%l)
       call removeAttribute(thisnode,"l")  
        else
        write(*,*)"Parser ERROR: The element 'lorb' requires the attribute 'l' to be defined."
        write(*,*)"stopped"
        stop
        
endif

            len= countChildEmentsWithName(thisnode,"wf")
     
allocate(getstructlorb%wfarray(len))
Do i=0,len-1
getstructlorb%wfarray(i+1)%wf=>getstructwf(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"wf"),0)) ) 
enddo

      i=0
      len=0
      call  handleunknownnodes(thisnode)
end function

function getstructspdb(thisnode)

implicit none
type(Node),pointer::thisnode
type(spdb_type),pointer::getstructspdb

integer::len=1,i=0
allocate(getstructspdb)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at spdb"
#endif
      
            len= countChildEmentsWithName(thisnode,"sp")
     
allocate(getstructspdb%sparray(len))
Do i=0,len-1
getstructspdb%sparray(i+1)%sp=>getstructsp(&
removeChild(thisnode,item(getElementsByTagname(thisnode,&
"sp"),0)) ) 
enddo

      i=0
      len=0
      call  handleunknownnodes(thisnode)
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
    
end module

