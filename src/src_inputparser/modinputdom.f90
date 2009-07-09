module inputdom
	use FoX_dom

  	implicit none
 	type(Node), pointer :: doc,inputnp,nullnode,emptynode,dummy
    type(DOMConfiguration),pointer :: config
    logical::parseerror
! Request full canonicalization
! ie convert CDATA sections to text sections, remove all entity references etc.



contains

subroutine loadinputDOM()
  config => newDOMConfig()
  parseerror=.false.
  call setParameter(config, "validate-if-schema", .true.)
  doc => parseFile("input.xml",config)
  inputnp=>getDocumentElement(doc)
  nullnode =>getattributenode(inputnp,"schouldneverexist")
  parseerror=.false.
   dummy => createDocument(getImplementation(), "", "info", null())
   emptynode=>createElementNS(dummy, "", "empty")
end subroutine

subroutine handleunknownnodes(np)
	type(Node),pointer::np,unknownnode
	type(NodeList),pointer::nl
	type( NamedNodeMap),pointer::nnm
	integer:: len,i
	len=0
	nnm=>getAttributes(np)
	len=getLength(nnm)
	if(.not. getNodeName(np) .eq. "input")then
	if(len.gt.0) then
		Do i=0, len-1
				unknownnode=>item(nnm,i)
				write(*,*) ">>>>> unrecognized attribute:"
				write(*,*) 	"      ",getName(unknownnode),'="',getValue(unknownnode),'" in element ',getNodeName(np)
				parseerror=.true.
		end do
	endif
	endif
	nl=>getChildNodes(np)
	len=getLength(nl)
	if(len.gt.0) then

		Do i=0, len-1
			unknownnode=>item(nl,i)
			select case(getnodetype(unknownnode))
			case(ELEMENT_NODE)
			write(*,*) ">>>>> unrecognized element:"
			write(*,*) 	"     <",getNodeName(unknownnode),"> in element <",getNodeName(np),">"
			parseerror=.true.
			case(TEXT_NODE)
			end select
		end do
	endif
end subroutine

subroutine ifparseerrorstop()
if(parseerror) then
write(*,*)"Stopping because of parse error"
stop
endif
end subroutine


subroutine destroyDOM()
 call destroy(doc)
 call destroy(config)
end subroutine

end module
