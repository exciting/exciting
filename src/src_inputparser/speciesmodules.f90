!
Module modsp
      Use inputdom
      Implicit None
      Type wf_type
         Real (8) :: trialEnergy
         Integer :: matchingOrder
         Logical :: searchE
      End Type
!
      Type wf_type_array
         Type (wf_type), Pointer :: wf
      End Type
      Type sp_type
         Character (1024) :: chemicalSymbol
         Real (8) :: z
         Real (8) :: mass
         Character (512) :: name
         Type (muffinTin_type), Pointer :: muffinTin
         Type (atomicState_type_array), Pointer :: atomicStatearray (:)
         Type (basis_type), Pointer :: basis
         Type (lorb_type_array), Pointer :: lorbarray (:)
      End Type
!
      Type sp_type_array
         Type (sp_type), Pointer :: sp
      End Type
      Type muffinTin_type
         Real (8) :: rmin
         Real (8) :: rinf
         Real (8) :: radius
         Integer :: radialmeshPoints
      End Type
      Type atomicState_type
         Integer :: n
         Integer :: l
         Integer :: kappa
         Real (8) :: occ
         Logical :: core
      End Type
!
      Type atomicState_type_array
         Type (atomicState_type), Pointer :: atomicState
      End Type
      Type basis_type
         Integer :: order
         Type (wf_type_array), Pointer :: wfarray (:)
         Type (exception_type_array), Pointer :: exceptionarray (:)
      End Type
      Type exception_type
         Integer :: l
         Type (wf_type_array), Pointer :: wfarray (:)
      End Type
!
      Type exception_type_array
         Type (exception_type), Pointer :: exception
      End Type
      Type lorb_type
         Integer :: l
         Type (wf_type_array), Pointer :: wfarray (:)
      End Type
!
      Type lorb_type_array
         Type (lorb_type), Pointer :: lorb
      End Type
      Type spdb_type
         Type (sp_type_array), Pointer :: sparray (:)
      End Type
!
      Type (sp_type) :: sp
Contains
!
      Function getstructwf (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (wf_type), Pointer :: getstructwf
         Type (Node), Pointer :: np
!
!
         Integer :: Len = 1, i = 0
         Allocate (getstructwf)
#ifdef INPUTDEBUG
         Write (*,*) "we are at wf"
#endif
!
         Nullify (np)
         np => getAttributeNode (thisnode, "trialEnergy")
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "trialEnergy", &
           & getstructwf%trialEnergy)
            Call removeAttribute (thisnode, "trialEnergy")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "matchingOrder")
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "matchingOrder", &
           & getstructwf%matchingOrder)
            Call removeAttribute (thisnode, "matchingOrder")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "searchE")
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "searchE", &
           & getstructwf%searchE)
            Call removeAttribute (thisnode, "searchE")
         End If
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructsp (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (sp_type), Pointer :: getstructsp
         Type (Node), Pointer :: np
!
!
         Integer :: Len = 1, i = 0
         Allocate (getstructsp)
#ifdef INPUTDEBUG
         Write (*,*) "we are at sp"
#endif
!
         Nullify (np)
         np => getAttributeNode (thisnode, "chemicalSymbol")
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "chemicalSymbol", &
           & getstructsp%chemicalSymbol)
            Call removeAttribute (thisnode, "chemicalSymbol")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "z")
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "z", getstructsp%z)
            Call removeAttribute (thisnode, "z")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "mass")
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "mass", &
           & getstructsp%mass)
            Call removeAttribute (thisnode, "mass")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "name")
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "name", &
           & getstructsp%name)
            Call removeAttribute (thisnode, "name")
         End If
!
         Len = countChildEmentsWithName (thisnode, "muffinTin")
         getstructsp%muffinTin => null ()
         Do i = 0, len - 1
            getstructsp%muffinTin => getstructmuffinTin &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "muffinTin"), 0)))
         End Do
!
         Len = countChildEmentsWithName (thisnode, "atomicState")
!
         Allocate (getstructsp%atomicStatearray(len))
         Do i = 0, len - 1
            getstructsp%atomicStatearray(i+1)%atomicState => &
           & getstructatomicState (removeChild(thisnode, &
           & item(getElementsByTagname(thisnode, "atomicState"), 0)))
         End Do
!
         Len = countChildEmentsWithName (thisnode, "basis")
         getstructsp%basis => null ()
         Do i = 0, len - 1
            getstructsp%basis => getstructbasis (removeChild(thisnode, &
           & item(getElementsByTagname(thisnode, "basis"), 0)))
         End Do
!
         Len = countChildEmentsWithName (thisnode, "lorb")
!
         Allocate (getstructsp%lorbarray(len))
         Do i = 0, len - 1
            getstructsp%lorbarray(i+1)%lorb => getstructlorb &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "lorb"), 0)))
         End Do
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructmuffinTin (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (muffinTin_type), Pointer :: getstructmuffinTin
         Type (Node), Pointer :: np
!
!
         Integer :: Len = 1, i = 0
         Allocate (getstructmuffinTin)
#ifdef INPUTDEBUG
         Write (*,*) "we are at muffinTin"
#endif
!
         Nullify (np)
         np => getAttributeNode (thisnode, "rmin")
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "rmin", &
           & getstructmuffinTin%rmin)
            Call removeAttribute (thisnode, "rmin")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "rinf")
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "rinf", &
           & getstructmuffinTin%rinf)
            Call removeAttribute (thisnode, "rinf")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "radius")
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "radius", &
           & getstructmuffinTin%radius)
            Call removeAttribute (thisnode, "radius")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "radialmeshPoints")
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "radialmeshPoints", &
           & getstructmuffinTin%radialmeshPoints)
            Call removeAttribute (thisnode, "radialmeshPoints")
         End If
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructatomicState (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (atomicState_type), Pointer :: getstructatomicState
         Type (Node), Pointer :: np
!
!
         Integer :: Len = 1, i = 0
         Allocate (getstructatomicState)
#ifdef INPUTDEBUG
         Write (*,*) "we are at atomicState"
#endif
!
         Nullify (np)
         np => getAttributeNode (thisnode, "n")
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "n", &
           & getstructatomicState%n)
            Call removeAttribute (thisnode, "n")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "l")
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "l", &
           & getstructatomicState%l)
            Call removeAttribute (thisnode, "l")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "kappa")
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "kappa", &
           & getstructatomicState%kappa)
            Call removeAttribute (thisnode, "kappa")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "occ")
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "occ", &
           & getstructatomicState%occ)
            Call removeAttribute (thisnode, "occ")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "core")
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "core", &
           & getstructatomicState%core)
            Call removeAttribute (thisnode, "core")
         End If
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructbasis (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (basis_type), Pointer :: getstructbasis
         Type (Node), Pointer :: np
!
!
         Integer :: Len = 1, i = 0
         Allocate (getstructbasis)
#ifdef INPUTDEBUG
         Write (*,*) "we are at basis"
#endif
!
         Nullify (np)
         np => getAttributeNode (thisnode, "order")
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "order", &
           & getstructbasis%order)
            Call removeAttribute (thisnode, "order")
         End If
!
         Len = countChildEmentsWithName (thisnode, "wf")
!
         Allocate (getstructbasis%wfarray(len))
         Do i = 0, len - 1
            getstructbasis%wfarray(i+1)%wf => getstructwf &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "wf"), 0)))
         End Do
!
         Len = countChildEmentsWithName (thisnode, "exception")
!
         Allocate (getstructbasis%exceptionarray(len))
         Do i = 0, len - 1
            getstructbasis%exceptionarray(i+1)%exception => &
           & getstructexception (removeChild(thisnode, &
           & item(getElementsByTagname(thisnode, "exception"), 0)))
         End Do
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructexception (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (exception_type), Pointer :: getstructexception
         Type (Node), Pointer :: np
!
!
         Integer :: Len = 1, i = 0
         Allocate (getstructexception)
#ifdef INPUTDEBUG
         Write (*,*) "we are at exception"
#endif
!
         Nullify (np)
         np => getAttributeNode (thisnode, "l")
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "l", &
           & getstructexception%l)
            Call removeAttribute (thisnode, "l")
         End If
!
         Len = countChildEmentsWithName (thisnode, "wf")
!
         Allocate (getstructexception%wfarray(len))
         Do i = 0, len - 1
            getstructexception%wfarray(i+1)%wf => getstructwf &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "wf"), 0)))
         End Do
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructlorb (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (lorb_type), Pointer :: getstructlorb
         Type (Node), Pointer :: np
!
!
         Integer :: Len = 1, i = 0
         Allocate (getstructlorb)
#ifdef INPUTDEBUG
         Write (*,*) "we are at lorb"
#endif
!
         Nullify (np)
         np => getAttributeNode (thisnode, "l")
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "l", getstructlorb%l)
            Call removeAttribute (thisnode, "l")
         End If
!
         Len = countChildEmentsWithName (thisnode, "wf")
!
         Allocate (getstructlorb%wfarray(len))
         Do i = 0, len - 1
            getstructlorb%wfarray(i+1)%wf => getstructwf &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "wf"), 0)))
         End Do
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructspdb (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (spdb_type), Pointer :: getstructspdb
		
         Integer :: Len = 1, i = 0
         Allocate (getstructspdb)
#ifdef INPUTDEBUG
         Write (*,*) "we are at spdb"
#endif
!
         Len = countChildEmentsWithName (thisnode, "sp")
!
         Allocate (getstructspdb%sparray(len))
         Do i = 0, len - 1
            getstructspdb%sparray(i+1)%sp => getstructsp &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "sp"), 0)))
         End Do
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
!
      Function countChildEmentsWithName (nodep, name)
         Implicit None
         Integer :: countChildEmentsWithName
         Type (Node), Pointer, Intent (In) :: nodep
         Character (Len=*), Intent (In) :: name
         Type (NodeList), Pointer :: children
         Type (Node), Pointer :: child
!
         Integer :: i
         children => getChildNodes (nodep)
         countChildEmentsWithName = 0
         Do i = 0, getlength (children) - 1
            child => item (children, i)
            If (name .Eq. getNodeName(child)) countChildEmentsWithName &
           & = countChildEmentsWithName + 1
         End Do
!
      End Function
!
End Module
!
