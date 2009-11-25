!
Module modinput
      Use inputdom
      Implicit None
      Type origin_type
         Real (8) :: coord (3)
      End Type
      Type point_type
         Real (8) :: coord (3)
         Character (512) :: label
      End Type
!
      Type point_type_array
         Type (point_type), Pointer :: point
      End Type
      Type plot1d_type
         Type (path_type), Pointer :: path
      End Type
      Type path_type
         Integer :: steps
         Character (512) :: outfileprefix
         Type (point_type_array), Pointer :: pointarray (:)
      End Type
      Type plot2d_type
         Type (parallelogram_type), Pointer :: parallelogram
      End Type
      Type parallelogram_type
         Integer :: grid (2)
         Character (512) :: outfileprefix
         Type (origin_type), Pointer :: origin
         Type (point_type_array), Pointer :: pointarray (:)
      End Type
      Type plot3d_type
         Type (box_type), Pointer :: box
      End Type
      Type box_type
         Integer :: grid (3)
         Character (512) :: outfileprefix
         Type (origin_type), Pointer :: origin
         Type (point_type_array), Pointer :: pointarray (:)
      End Type
      Type kstlist_type
         Integer, Pointer :: pointstatepair (:, :)
      End Type
      Type inputset_type
         Type (input_type_array), Pointer :: inputarray (:)
      End Type
      Type input_type
         Character (1024) :: xsltpath
         Character (1024) :: scratchpath
         Character (1024) :: id
         Character (1024) :: depends
         Character (512) :: title
         Type (structure_type), Pointer :: structure
         Type (groundstate_type), Pointer :: groundstate
         Type (structureoptimization_type), Pointer :: &
        & structureoptimization
         Type (properties_type), Pointer :: properties
         Type (phonons_type), Pointer :: phonons
         Type (xs_type), Pointer :: xs
      End Type
!
      Type input_type_array
         Type (input_type), Pointer :: input
      End Type
      Type structure_type
         Character (1024) :: speciespath
         Logical :: molecule
         Real (8) :: vacuum
         Real (8) :: epslat
         Logical :: autormt
         Logical :: primcell
         Logical :: tshift
         Type (symmetries_type), Pointer :: symmetries
         Type (crystal_type), Pointer :: crystal
         Type (species_type_array), Pointer :: speciesarray (:)
      End Type
      Type symmetries_type
         Character (512) :: HermannMauguinSymbol
         Character (512) :: HallSymbol
         Character (512) :: SchoenfliesSymbol
         Character (512) :: spaceGroupNumber
         Type (lattice_type), Pointer :: lattice
         Type (WyckoffPositions_type), Pointer :: WyckoffPositions
      End Type
      Type lattice_type
         Real (8) :: a
         Real (8) :: b
         Real (8) :: c
         Real (8) :: ab
         Real (8) :: ac
         Real (8) :: bc
         Integer :: ncell (3)
      End Type
      Type WyckoffPositions_type
         Type (wspecies_type_array), Pointer :: wspeciesarray (:)
      End Type
      Type wspecies_type
         Character (512) :: speciesfile
         Type (wpos_type_array), Pointer :: wposarray (:)
      End Type
!
      Type wspecies_type_array
         Type (wspecies_type), Pointer :: wspecies
      End Type
      Type wpos_type
         Real (8) :: coord (3)
      End Type
!
      Type wpos_type_array
         Type (wpos_type), Pointer :: wpos
      End Type
      Type crystal_type
         Real (8) :: scale
         Real (8) :: stretch (3)
         Real (8), Pointer :: basevect (:, :)
      End Type
      Type species_type
         Character (1024) :: speciesfile
         Character (512) :: chemicalSymbol
         Integer :: atomicNumber
         Real (8) :: rmt
         Type (atom_type_array), Pointer :: atomarray (:)
         Type (LDAplusu_type), Pointer :: LDAplusu
      End Type
!
      Type species_type_array
         Type (species_type), Pointer :: species
      End Type
      Type atom_type
         Real (8) :: coord (3)
         Real (8) :: bfcmt (3)
         Real (8) :: mommtfix (3)
      End Type
!
      Type atom_type_array
         Type (atom_type), Pointer :: atom
      End Type
      Type LDAplusu_type
         Real (8) :: L
         Real (8) :: U
         Real (8) :: J
      End Type
      Type groundstate_type
         Character (512) :: do
         Integer :: donumber
         Integer :: ngridk (3)
         Real (8) :: rgkmax
         Real (8) :: epspot
         Real (8) :: rmtapm (2)
         Real (8) :: swidth
         Character (512) :: stype
         Integer :: stypenumber
         Character (512) :: findlinentype
         Integer :: findlinentypenumber
         Integer :: isgkmax
         Real (8) :: gmaxvr
         Integer :: nempty
         Logical :: nosym
         Logical :: frozencore
         Logical :: autokpt
         Real (8) :: radkpt
         Logical :: reducek
         Logical :: tfibs
         Logical :: tforce
         Integer :: lmaxapw
         Integer :: maxscl
         Real (8) :: chgexs
         Real (8) :: deband
         Real (8) :: epschg
         Real (8) :: epsocc
         Character (512) :: mixer
         Integer :: mixernumber
         Real (8) :: beta0
         Real (8) :: betainc
         Real (8) :: betadec
         Integer :: lradstep
         Integer :: nprad
         Character (512) :: xctype
         Integer :: xctypenumber
         Real (8) :: evalmin
         Integer :: lmaxvr
         Real (8) :: fracinr
         Integer :: lmaxinr
         Integer :: lmaxmat
         Real (8) :: vkloff (3)
         Integer :: npsden
         Real (8) :: cfdamp
         Logical :: nosource
         Logical :: tevecsv
         Integer :: nwrite
         Logical :: ptnucl
         Type (spin_type), Pointer :: spin
         Type (HartreeFock_type), Pointer :: HartreeFock
         Type (solver_type), Pointer :: solver
         Type (OEP_type), Pointer :: OEP
         Type (RDMFT_type), Pointer :: RDMFT
      End Type
      Type spin_type
         Real (8) :: momfix (3)
         Real (8) :: bfieldc (3)
         Logical :: spinorb
         Logical :: spinsprl
         Real (8) :: vqlss (3)
         Real (8) :: taufsm
         Real (8) :: reducebf
         Character (512) :: fixspin
         Integer :: fixspinnumber
      End Type
      Type HartreeFock_type
         Real (8) :: epsengy
      End Type
      Type solver_type
         Character (512) :: Type
         Integer :: typenumber
         Logical :: packedmatrixstorage
         Real (8) :: epsarpack
         Real (8) :: evaltol
      End Type
      Type OEP_type
         Integer :: maxitoep
         Real (8) :: tauoep (3)
      End Type
      Type RDMFT_type
         Integer :: rdmxctype
         Integer :: rdmmaxscl
         Integer :: maxitn
         Integer :: maxitc
         Real (8) :: taurdmn
         Real (8) :: taurdmc
         Real (8) :: rdmalpha
         Real (8) :: rdmtemp
      End Type
      Type structureoptimization_type
         Real (8) :: epsforce
         Real (8) :: tau0atm
         Logical :: resume
      End Type
      Type properties_type
         Type (bandstructure_type), Pointer :: bandstructure
         Type (STM_type), Pointer :: STM
         Type (wfplot_type), Pointer :: wfplot
         Type (dos_type), Pointer :: dos
         Type (LSJ_type), Pointer :: LSJ
         Type (masstensor_type), Pointer :: masstensor
         Type (chargedensityplot_type), Pointer :: chargedensityplot
         Type (exccplot_type), Pointer :: exccplot
         Type (elfplot_type), Pointer :: elfplot
         Type (mvecfield_type), Pointer :: mvecfield
         Type (xcmvecfield_type), Pointer :: xcmvecfield
         Type (electricfield_type), Pointer :: electricfield
         Type (gradmvecfield_type), Pointer :: gradmvecfield
         Type (fermisurfaceplot_type), Pointer :: fermisurfaceplot
         Type (EFG_type), Pointer :: EFG
         Type (momentummatrix_type), Pointer :: momentummatrix
         Type (linresponsetensor_type), Pointer :: linresponsetensor
         Type (mossbauer_type), Pointer :: mossbauer
         Type (dielectric_type), Pointer :: dielectric
         Type (expiqr_type), Pointer :: expiqr
         Type (elnes_type), Pointer :: elnes
         Type (eliashberg_type), Pointer :: eliashberg
      End Type
      Type bandstructure_type
         Real (8) :: scissor
         Logical :: Character
         Type (plot1d_type), Pointer :: plot1d
      End Type
      Type STM_type
         Type (plot2d_type), Pointer :: plot2d
      End Type
      Type wfplot_type
         Type (kstlist_type), Pointer :: kstlist
         Type (plot1d_type), Pointer :: plot1d
         Type (plot2d_type), Pointer :: plot2d
         Type (plot3d_type), Pointer :: plot3d
      End Type
      Type dos_type
         Real (8) :: sqados (3)
         Logical :: lmirep
         Integer :: nwdos
         Integer :: ngrdos
         Real (8) :: scissor
         Integer :: nsmdos
         Real (8) :: winddos (2)
      End Type
      Type LSJ_type
         Type (kstlist_type), Pointer :: kstlist
      End Type
      Type masstensor_type
         Real (8) :: deltaem
         Integer :: ndspem
         Real (8) :: vklem (3)
      End Type
      Type chargedensityplot_type
         Type (plot1d_type), Pointer :: plot1d
         Type (plot2d_type), Pointer :: plot2d
         Type (plot3d_type), Pointer :: plot3d
      End Type
      Type exccplot_type
         Type (plot1d_type), Pointer :: plot1d
         Type (plot2d_type), Pointer :: plot2d
         Type (plot3d_type), Pointer :: plot3d
      End Type
      Type elfplot_type
         Type (plot1d_type), Pointer :: plot1d
         Type (plot2d_type), Pointer :: plot2d
         Type (plot3d_type), Pointer :: plot3d
      End Type
      Type mvecfield_type
         Type (plot2d_type), Pointer :: plot2d
         Type (plot3d_type), Pointer :: plot3d
      End Type
      Type xcmvecfield_type
         Type (plot2d_type), Pointer :: plot2d
         Type (plot3d_type), Pointer :: plot3d
      End Type
      Type electricfield_type
         Type (plot2d_type), Pointer :: plot2d
         Type (plot3d_type), Pointer :: plot3d
      End Type
      Type gradmvecfield_type
         Type (plot1d_type), Pointer :: plot1d
         Type (plot2d_type), Pointer :: plot2d
         Type (plot3d_type), Pointer :: plot3d
      End Type
      Type fermisurfaceplot_type
         Integer :: nstfsp
         Logical :: separate
      End Type
!
      Type EFG_type
         Logical :: exists
      End Type
!
      Type momentummatrix_type
         Logical :: exists
      End Type
      Type linresponsetensor_type
         Real (8) :: scissor
         Integer, Pointer :: optcomp (:, :)
      End Type
!
      Type mossbauer_type
         Logical :: exists
      End Type
!
      Type dielectric_type
         Logical :: exists
      End Type
!
      Type expiqr_type
         Logical :: exists
      End Type
      Type elnes_type
         Real (8) :: vecql (3)
      End Type
      Type eliashberg_type
         Real (8) :: mustar
      End Type
      Type phonons_type
         Logical :: reduceq
         Real (8) :: deltaph
         Type (qpointset_type), Pointer :: qpointset
         Type (phonondos_type), Pointer :: phonondos
         Type (phonondispplot_type), Pointer :: phonondispplot
      End Type
!
      Type phonondos_type
         Logical :: exists
      End Type
      Type phonondispplot_type
         Type (plot1d_type), Pointer :: plot1d
      End Type
      Type xs_type
         Integer :: emattype
         Logical :: dfoffdiag
         Integer :: lmaxapwwf
         Integer :: lmaxemat
         Real (8) :: emaxdf
         Real (8) :: broad
         Logical :: tevout
         Character (512) :: xstype
         Integer :: xstypenumber
         Logical :: symmorph
         Logical :: fastpmat
         Logical :: fastemat
         Logical :: gather
         Logical :: tappinfo
         Integer :: dbglev
         Logical :: usegdft
         Real (8) :: gqmax
         Logical :: nosym
         Integer :: ngridk (3)
         Real (8) :: vkloff (3)
         Logical :: reducek
         Integer :: ngridq (3)
         Logical :: reduceq
         Real (8) :: rgkmax
         Real (8) :: swidth
         Integer :: lmaxapw
         Integer :: lmaxmat
         Integer :: nempty
         Real (8) :: scissor
         Type (tddft_type), Pointer :: tddft
         Type (screening_type), Pointer :: screening
         Type (BSE_type), Pointer :: BSE
         Type (qpointset_type), Pointer :: qpointset
         Type (tetra_type), Pointer :: tetra
         Type (dosWindow_type), Pointer :: dosWindow
         Type (plan_type), Pointer :: plan
      End Type
      Type tddft_type
         Logical :: intraband
         Logical :: torddf
         Logical :: tordfxc
         Logical :: aresdf
         Logical :: aresfxc
         Real (8) :: fxcbsesplit
         Logical :: acont
         Integer :: nwacont
         Logical :: lindhard
         Real (8) :: epsdfde
         Logical :: kerndiag
         Integer :: lmaxalda
         Real (8) :: alphalrc
         Real (8) :: alphalrcdyn
         Real (8) :: betalrcdyn
         Integer :: mdfqtype
         Character (512) :: fxctype
         Integer :: fxctypenumber
         Logical :: resumefromkernel
         Type (dftrans_type), Pointer :: dftrans
      End Type
      Type dftrans_type
         Integer, Pointer :: trans (:, :)
      End Type
      Type screening_type
         Character (512) :: run
         Integer :: runnumber
         Logical :: nosym
         Integer :: ngridk (3)
         Logical :: reducek
         Real (8) :: vkloff (3)
         Real (8) :: rgkmax
         Integer :: nempty
         Character (512) :: screentype
         Integer :: screentypenumber
      End Type
      Type BSE_type
         Logical :: nosym
         Logical :: reducek
         Real (8) :: vkloff (3)
         Real (8) :: rgkmax
         Integer :: scrherm
         Logical :: fbzq
         Character (512) :: sciavtype
         Integer :: sciavtypenumber
         Logical :: sciavbd
         Logical :: sciavqhd
         Logical :: sciavqwg
         Logical :: sciavqbd
         Logical :: bsedirsing
         Integer :: lmaxdielt
         Integer :: nleblaik
         Integer :: nexcitmax
         Integer :: nstlbse (2)
         Integer :: nstlce (2)
         Character (512) :: bsetype
         Integer :: bsetypenumber
      End Type
      Type tetra_type
         Logical :: tetraocc
         Logical :: tetradf
         Logical :: kordexc
         Logical :: cw1k
         Integer :: qweights
      End Type
      Type dosWindow_type
         Integer :: points
         Real (8) :: intv (2)
         Integer :: nsmdos
      End Type
      Type plan_type
         Type (doonly_type_array), Pointer :: doonlyarray (:)
      End Type
      Type doonly_type
         Character (512) :: task
         Integer :: tasknumber
      End Type
!
      Type doonly_type_array
         Type (doonly_type), Pointer :: doonly
      End Type
      Type qpointset_type
         Real (8), Pointer :: qpoint (:, :)
      End Type
!
      Type (input_type) :: input
Contains
!
      Function getstructorigin (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (origin_type), Pointer :: getstructorigin
         Type (Node), Pointer :: np
!
!
         Integer :: Len = 1, i = 0
         Allocate (getstructorigin)
#ifdef INPUTDEBUG
         Write (*,*) "we are at origin"
#endif
!
         Nullify (np)
         np => getAttributeNode (thisnode, "coord")
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "coord", &
           & getstructorigin%coord)
            Call removeAttribute (thisnode, "coord")
         End If
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructpoint (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (point_type), Pointer :: getstructpoint
         Type (Node), Pointer :: np
!
!
         Integer :: Len = 1, i = 0
         Allocate (getstructpoint)
#ifdef INPUTDEBUG
         Write (*,*) "we are at point"
#endif
!
         Nullify (np)
         np => getAttributeNode (thisnode, "coord")
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "coord", &
           & getstructpoint%coord)
            Call removeAttribute (thisnode, "coord")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "label")
         getstructpoint%label = ""
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "label", &
           & getstructpoint%label)
            Call removeAttribute (thisnode, "label")
         End If
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructplot1d (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (plot1d_type), Pointer :: getstructplot1d
		
         Integer :: Len = 1, i = 0
         Allocate (getstructplot1d)
#ifdef INPUTDEBUG
         Write (*,*) "we are at plot1d"
#endif
!
         Len = countChildEmentsWithName (thisnode, "path")
         getstructplot1d%path => null ()
         Do i = 0, len - 1
            getstructplot1d%path => getstructpath &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "path"), 0)))
         End Do
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructpath (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (path_type), Pointer :: getstructpath
         Type (Node), Pointer :: np
!
!
         Integer :: Len = 1, i = 0
         Allocate (getstructpath)
#ifdef INPUTDEBUG
         Write (*,*) "we are at path"
#endif
!
         Nullify (np)
         np => getAttributeNode (thisnode, "steps")
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "steps", &
           & getstructpath%steps)
            Call removeAttribute (thisnode, "steps")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "outfileprefix")
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "outfileprefix", &
           & getstructpath%outfileprefix)
            Call removeAttribute (thisnode, "outfileprefix")
         End If
!
         Len = countChildEmentsWithName (thisnode, "point")
!
         Allocate (getstructpath%pointarray(len))
         Do i = 0, len - 1
            getstructpath%pointarray(i+1)%point => getstructpoint &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "point"), 0)))
         End Do
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructplot2d (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (plot2d_type), Pointer :: getstructplot2d
		
         Integer :: Len = 1, i = 0
         Allocate (getstructplot2d)
#ifdef INPUTDEBUG
         Write (*,*) "we are at plot2d"
#endif
!
         Len = countChildEmentsWithName (thisnode, "parallelogram")
         getstructplot2d%parallelogram => null ()
         Do i = 0, len - 1
            getstructplot2d%parallelogram => getstructparallelogram &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "parallelogram"), 0)))
         End Do
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructparallelogram (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (parallelogram_type), Pointer :: getstructparallelogram
         Type (Node), Pointer :: np
!
!
         Integer :: Len = 1, i = 0
         Allocate (getstructparallelogram)
#ifdef INPUTDEBUG
         Write (*,*) "we are at parallelogram"
#endif
!
         Nullify (np)
         np => getAttributeNode (thisnode, "grid")
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "grid", &
           & getstructparallelogram%grid)
            Call removeAttribute (thisnode, "grid")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "outfileprefix")
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "outfileprefix", &
           & getstructparallelogram%outfileprefix)
            Call removeAttribute (thisnode, "outfileprefix")
         End If
!
         Len = countChildEmentsWithName (thisnode, "origin")
         getstructparallelogram%origin => null ()
         Do i = 0, len - 1
            getstructparallelogram%origin => getstructorigin &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "origin"), 0)))
         End Do
!
         Len = countChildEmentsWithName (thisnode, "point")
!
         Allocate (getstructparallelogram%pointarray(len))
         Do i = 0, len - 1
            getstructparallelogram%pointarray(i+1)%point => &
           & getstructpoint (removeChild(thisnode, &
           & item(getElementsByTagname(thisnode, "point"), 0)))
         End Do
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructplot3d (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (plot3d_type), Pointer :: getstructplot3d
		
         Integer :: Len = 1, i = 0
         Allocate (getstructplot3d)
#ifdef INPUTDEBUG
         Write (*,*) "we are at plot3d"
#endif
!
         Len = countChildEmentsWithName (thisnode, "box")
         getstructplot3d%box => null ()
         Do i = 0, len - 1
            getstructplot3d%box => getstructbox (removeChild(thisnode, &
           & item(getElementsByTagname(thisnode, "box"), 0)))
         End Do
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructbox (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (box_type), Pointer :: getstructbox
         Type (Node), Pointer :: np
!
!
         Integer :: Len = 1, i = 0
         Allocate (getstructbox)
#ifdef INPUTDEBUG
         Write (*,*) "we are at box"
#endif
!
         Nullify (np)
         np => getAttributeNode (thisnode, "grid")
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "grid", &
           & getstructbox%grid)
            Call removeAttribute (thisnode, "grid")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "outfileprefix")
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "outfileprefix", &
           & getstructbox%outfileprefix)
            Call removeAttribute (thisnode, "outfileprefix")
         End If
!
         Len = countChildEmentsWithName (thisnode, "origin")
         getstructbox%origin => null ()
         Do i = 0, len - 1
            getstructbox%origin => getstructorigin &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "origin"), 0)))
         End Do
!
         Len = countChildEmentsWithName (thisnode, "point")
!
         Allocate (getstructbox%pointarray(len))
         Do i = 0, len - 1
            getstructbox%pointarray(i+1)%point => getstructpoint &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "point"), 0)))
         End Do
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructkstlist (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (kstlist_type), Pointer :: getstructkstlist
		
         Integer :: Len = 1, i = 0
         Allocate (getstructkstlist)
#ifdef INPUTDEBUG
         Write (*,*) "we are at kstlist"
#endif
!
         Len = countChildEmentsWithName (thisnode, "pointstatepair")
         Allocate (getstructkstlist%pointstatepair(2, len))
         Do i = 1, len
!
            getstructkstlist%pointstatepair (:, i) = &
           & getvalueofpointstatepair (removeChild(thisnode, &
           & item(getElementsByTagname(thisnode, "pointstatepair"), &
           & 0)))
         End Do
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructinputset (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (inputset_type), Pointer :: getstructinputset
		
         Integer :: Len = 1, i = 0
         Allocate (getstructinputset)
#ifdef INPUTDEBUG
         Write (*,*) "we are at inputset"
#endif
!
         Len = countChildEmentsWithName (thisnode, "input")
!
         Allocate (getstructinputset%inputarray(len))
         Do i = 0, len - 1
            getstructinputset%inputarray(i+1)%input => getstructinput &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "input"), 0)))
         End Do
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructinput (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (input_type), Pointer :: getstructinput
         Type (Node), Pointer :: np
!
!
         Integer :: Len = 1, i = 0
         Allocate (getstructinput)
#ifdef INPUTDEBUG
         Write (*,*) "we are at input"
#endif
!
         Nullify (np)
         np => getAttributeNode (thisnode, "xsltpath")
         getstructinput%xsltpath = "../../xml"
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "xsltpath", &
           & getstructinput%xsltpath)
            Call removeAttribute (thisnode, "xsltpath")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "scratchpath")
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "scratchpath", &
           & getstructinput%scratchpath)
            Call removeAttribute (thisnode, "scratchpath")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "id")
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "id", &
           & getstructinput%id)
            Call removeAttribute (thisnode, "id")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "depends")
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "depends", &
           & getstructinput%depends)
            Call removeAttribute (thisnode, "depends")
         End If
!
         Len = countChildEmentsWithName (thisnode, "structure")
         getstructinput%structure => null ()
         Do i = 0, len - 1
            getstructinput%structure => getstructstructure &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "structure"), 0)))
         End Do
!
         Len = countChildEmentsWithName (thisnode, "groundstate")
         getstructinput%groundstate => null ()
         Do i = 0, len - 1
            getstructinput%groundstate => getstructgroundstate &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "groundstate"), 0)))
         End Do
!
         Len = countChildEmentsWithName (thisnode, "structureoptimizati&
        &on")
         getstructinput%structureoptimization => null ()
         Do i = 0, len - 1
            getstructinput%structureoptimization => &
           & getstructstructureoptimization (removeChild(thisnode, &
           & item(getElementsByTagname(thisnode, "structureoptimization&
           &"), 0)))
         End Do
!
         Len = countChildEmentsWithName (thisnode, "properties")
         getstructinput%properties => null ()
         Do i = 0, len - 1
            getstructinput%properties => getstructproperties &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "properties"), 0)))
         End Do
!
         Len = countChildEmentsWithName (thisnode, "phonons")
         getstructinput%phonons => null ()
         Do i = 0, len - 1
            getstructinput%phonons => getstructphonons &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "phonons"), 0)))
         End Do
!
         Len = countChildEmentsWithName (thisnode, "xs")
         getstructinput%xs => null ()
         Do i = 0, len - 1
            getstructinput%xs => getstructxs (removeChild(thisnode, &
           & item(getElementsByTagname(thisnode, "xs"), 0)))
         End Do
!
         Len = countChildEmentsWithName (thisnode, "title")
         Do i = 1, len
!
            getstructinput%title = getvalueoftitle &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "title"), 0)))
         End Do
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructstructure (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (structure_type), Pointer :: getstructstructure
         Type (Node), Pointer :: np
!
!
         Integer :: Len = 1, i = 0
         Allocate (getstructstructure)
#ifdef INPUTDEBUG
         Write (*,*) "we are at structure"
#endif
!
         Nullify (np)
         np => getAttributeNode (thisnode, "speciespath")
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "speciespath", &
           & getstructstructure%speciespath)
            Call removeAttribute (thisnode, "speciespath")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "molecule")
         getstructstructure%molecule = .False.
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "molecule", &
           & getstructstructure%molecule)
            Call removeAttribute (thisnode, "molecule")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "vacuum")
         getstructstructure%vacuum = 10
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "vacuum", &
           & getstructstructure%vacuum)
            Call removeAttribute (thisnode, "vacuum")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "epslat")
         getstructstructure%epslat = 1e-6
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "epslat", &
           & getstructstructure%epslat)
            Call removeAttribute (thisnode, "epslat")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "autormt")
         getstructstructure%autormt = .False.
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "autormt", &
           & getstructstructure%autormt)
            Call removeAttribute (thisnode, "autormt")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "primcell")
         getstructstructure%primcell = .False.
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "primcell", &
           & getstructstructure%primcell)
            Call removeAttribute (thisnode, "primcell")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "tshift")
         getstructstructure%tshift = .True.
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "tshift", &
           & getstructstructure%tshift)
            Call removeAttribute (thisnode, "tshift")
         End If
!
         Len = countChildEmentsWithName (thisnode, "symmetries")
         getstructstructure%symmetries => null ()
         Do i = 0, len - 1
            getstructstructure%symmetries => getstructsymmetries &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "symmetries"), 0)))
         End Do
!
         Len = countChildEmentsWithName (thisnode, "crystal")
         getstructstructure%crystal => null ()
         Do i = 0, len - 1
            getstructstructure%crystal => getstructcrystal &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "crystal"), 0)))
         End Do
!
         Len = countChildEmentsWithName (thisnode, "species")
!
         Allocate (getstructstructure%speciesarray(len))
         Do i = 0, len - 1
            getstructstructure%speciesarray(i+1)%species => &
           & getstructspecies (removeChild(thisnode, &
           & item(getElementsByTagname(thisnode, "species"), 0)))
         End Do
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructsymmetries (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (symmetries_type), Pointer :: getstructsymmetries
         Type (Node), Pointer :: np
!
!
         Integer :: Len = 1, i = 0
         Allocate (getstructsymmetries)
#ifdef INPUTDEBUG
         Write (*,*) "we are at symmetries"
#endif
!
         Nullify (np)
         np => getAttributeNode (thisnode, "HermannMauguinSymbol")
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "HermannMauguinSymbol",&
           &  getstructsymmetries%HermannMauguinSymbol)
            Call removeAttribute (thisnode, "HermannMauguinSymbol")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "HallSymbol")
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "HallSymbol", &
           & getstructsymmetries%HallSymbol)
            Call removeAttribute (thisnode, "HallSymbol")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "SchoenfliesSymbol")
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "SchoenfliesSymbol", &
           & getstructsymmetries%SchoenfliesSymbol)
            Call removeAttribute (thisnode, "SchoenfliesSymbol")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "spaceGroupNumber")
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "spaceGroupNumber", &
           & getstructsymmetries%spaceGroupNumber)
            Call removeAttribute (thisnode, "spaceGroupNumber")
         End If
!
         Len = countChildEmentsWithName (thisnode, "lattice")
         getstructsymmetries%lattice => null ()
         Do i = 0, len - 1
            getstructsymmetries%lattice => getstructlattice &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "lattice"), 0)))
         End Do
!
         Len = countChildEmentsWithName (thisnode, "WyckoffPositions")
         getstructsymmetries%WyckoffPositions => null ()
         Do i = 0, len - 1
            getstructsymmetries%WyckoffPositions => &
           & getstructWyckoffPositions (removeChild(thisnode, &
           & item(getElementsByTagname(thisnode, "WyckoffPositions"), &
           & 0)))
         End Do
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructlattice (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (lattice_type), Pointer :: getstructlattice
         Type (Node), Pointer :: np
!
!
         Integer :: Len = 1, i = 0
         Allocate (getstructlattice)
#ifdef INPUTDEBUG
         Write (*,*) "we are at lattice"
#endif
!
         Nullify (np)
         np => getAttributeNode (thisnode, "a")
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "a", &
           & getstructlattice%a)
            Call removeAttribute (thisnode, "a")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "b")
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "b", &
           & getstructlattice%b)
            Call removeAttribute (thisnode, "b")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "c")
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "c", &
           & getstructlattice%c)
            Call removeAttribute (thisnode, "c")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "ab")
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "ab", &
           & getstructlattice%ab)
            Call removeAttribute (thisnode, "ab")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "ac")
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "ac", &
           & getstructlattice%ac)
            Call removeAttribute (thisnode, "ac")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "bc")
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "bc", &
           & getstructlattice%bc)
            Call removeAttribute (thisnode, "bc")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "ncell")
         getstructlattice%ncell = (/ 1, 1, 1 /)
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "ncell", &
           & getstructlattice%ncell)
            Call removeAttribute (thisnode, "ncell")
         End If
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructWyckoffPositions (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (WyckoffPositions_type), Pointer :: &
        & getstructWyckoffPositions
		
         Integer :: Len = 1, i = 0
         Allocate (getstructWyckoffPositions)
#ifdef INPUTDEBUG
         Write (*,*) "we are at WyckoffPositions"
#endif
!
         Len = countChildEmentsWithName (thisnode, "wspecies")
!
         Allocate (getstructWyckoffPositions%wspeciesarray(len))
         Do i = 0, len - 1
            getstructWyckoffPositions%wspeciesarray(i+1)%wspecies => &
           & getstructwspecies (removeChild(thisnode, &
           & item(getElementsByTagname(thisnode, "wspecies"), 0)))
         End Do
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructwspecies (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (wspecies_type), Pointer :: getstructwspecies
         Type (Node), Pointer :: np
!
!
         Integer :: Len = 1, i = 0
         Allocate (getstructwspecies)
#ifdef INPUTDEBUG
         Write (*,*) "we are at wspecies"
#endif
!
         Nullify (np)
         np => getAttributeNode (thisnode, "speciesfile")
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "speciesfile", &
           & getstructwspecies%speciesfile)
            Call removeAttribute (thisnode, "speciesfile")
         End If
!
         Len = countChildEmentsWithName (thisnode, "wpos")
!
         Allocate (getstructwspecies%wposarray(len))
         Do i = 0, len - 1
            getstructwspecies%wposarray(i+1)%wpos => getstructwpos &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "wpos"), 0)))
         End Do
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructwpos (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (wpos_type), Pointer :: getstructwpos
         Type (Node), Pointer :: np
!
!
         Integer :: Len = 1, i = 0
         Allocate (getstructwpos)
#ifdef INPUTDEBUG
         Write (*,*) "we are at wpos"
#endif
!
         Nullify (np)
         np => getAttributeNode (thisnode, "coord")
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "coord", &
           & getstructwpos%coord)
            Call removeAttribute (thisnode, "coord")
         End If
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructcrystal (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (crystal_type), Pointer :: getstructcrystal
         Type (Node), Pointer :: np
!
!
         Integer :: Len = 1, i = 0
         Allocate (getstructcrystal)
#ifdef INPUTDEBUG
         Write (*,*) "we are at crystal"
#endif
!
         Nullify (np)
         np => getAttributeNode (thisnode, "scale")
         getstructcrystal%scale = 1
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "scale", &
           & getstructcrystal%scale)
            Call removeAttribute (thisnode, "scale")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "stretch")
         getstructcrystal%stretch = (/ 1.0d0, 1.0d0, 1.0d0 /)
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "stretch", &
           & getstructcrystal%stretch)
            Call removeAttribute (thisnode, "stretch")
         End If
!
         Len = countChildEmentsWithName (thisnode, "basevect")
         Allocate (getstructcrystal%basevect(3, len))
         Do i = 1, len
!
            getstructcrystal%basevect (:, i) = getvalueofbasevect &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "basevect"), 0)))
         End Do
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructspecies (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (species_type), Pointer :: getstructspecies
         Type (Node), Pointer :: np
!
!
         Integer :: Len = 1, i = 0
         Allocate (getstructspecies)
#ifdef INPUTDEBUG
         Write (*,*) "we are at species"
#endif
!
         Nullify (np)
         np => getAttributeNode (thisnode, "speciesfile")
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "speciesfile", &
           & getstructspecies%speciesfile)
            Call removeAttribute (thisnode, "speciesfile")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "chemicalSymbol")
         getstructspecies%chemicalSymbol = ""
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "chemicalSymbol", &
           & getstructspecies%chemicalSymbol)
            Call removeAttribute (thisnode, "chemicalSymbol")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "atomicNumber")
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "atomicNumber", &
           & getstructspecies%atomicNumber)
            Call removeAttribute (thisnode, "atomicNumber")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "rmt")
         getstructspecies%rmt = - 1
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "rmt", &
           & getstructspecies%rmt)
            Call removeAttribute (thisnode, "rmt")
         End If
!
         Len = countChildEmentsWithName (thisnode, "atom")
!
         Allocate (getstructspecies%atomarray(len))
         Do i = 0, len - 1
            getstructspecies%atomarray(i+1)%atom => getstructatom &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "atom"), 0)))
         End Do
!
         Len = countChildEmentsWithName (thisnode, "LDAplusu")
         getstructspecies%LDAplusu => null ()
         Do i = 0, len - 1
            getstructspecies%LDAplusu => getstructLDAplusu &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "LDAplusu"), 0)))
         End Do
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructatom (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (atom_type), Pointer :: getstructatom
         Type (Node), Pointer :: np
!
!
         Integer :: Len = 1, i = 0
         Allocate (getstructatom)
#ifdef INPUTDEBUG
         Write (*,*) "we are at atom"
#endif
!
         Nullify (np)
         np => getAttributeNode (thisnode, "coord")
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "coord", &
           & getstructatom%coord)
            Call removeAttribute (thisnode, "coord")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "bfcmt")
         getstructatom%bfcmt = (/ 0, 0, 0 /)
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "bfcmt", &
           & getstructatom%bfcmt)
            Call removeAttribute (thisnode, "bfcmt")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "mommtfix")
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "mommtfix", &
           & getstructatom%mommtfix)
            Call removeAttribute (thisnode, "mommtfix")
         End If
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructLDAplusu (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (LDAplusu_type), Pointer :: getstructLDAplusu
         Type (Node), Pointer :: np
!
!
         Integer :: Len = 1, i = 0
         Allocate (getstructLDAplusu)
#ifdef INPUTDEBUG
         Write (*,*) "we are at LDAplusu"
#endif
!
         Nullify (np)
         np => getAttributeNode (thisnode, "L")
         getstructLDAplusu%L = 0
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "L", &
           & getstructLDAplusu%L)
            Call removeAttribute (thisnode, "L")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "U")
         getstructLDAplusu%U = 0
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "U", &
           & getstructLDAplusu%U)
            Call removeAttribute (thisnode, "U")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "J")
         getstructLDAplusu%J = 0
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "J", &
           & getstructLDAplusu%J)
            Call removeAttribute (thisnode, "J")
         End If
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructgroundstate (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (groundstate_type), Pointer :: getstructgroundstate
         Type (Node), Pointer :: np
!
!
         Integer :: Len = 1, i = 0
         Allocate (getstructgroundstate)
#ifdef INPUTDEBUG
         Write (*,*) "we are at groundstate"
#endif
!
         Nullify (np)
         np => getAttributeNode (thisnode, "do")
         getstructgroundstate%do = "fromscratch"
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "do", &
           & getstructgroundstate%do)
            Call removeAttribute (thisnode, "do")
         End If
         getstructgroundstate%donumber = stringtonumberdo &
        & (getstructgroundstate%do)
!
         Nullify (np)
         np => getAttributeNode (thisnode, "ngridk")
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "ngridk", &
           & getstructgroundstate%ngridk)
            Call removeAttribute (thisnode, "ngridk")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "rgkmax")
         getstructgroundstate%rgkmax = 7
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "rgkmax", &
           & getstructgroundstate%rgkmax)
            Call removeAttribute (thisnode, "rgkmax")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "epspot")
         getstructgroundstate%epspot = 1e-6
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "epspot", &
           & getstructgroundstate%epspot)
            Call removeAttribute (thisnode, "epspot")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "rmtapm")
         getstructgroundstate%rmtapm = (/ 0.25d0, 0.95d0 /)
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "rmtapm", &
           & getstructgroundstate%rmtapm)
            Call removeAttribute (thisnode, "rmtapm")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "swidth")
         getstructgroundstate%swidth = 0.001d0
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "swidth", &
           & getstructgroundstate%swidth)
            Call removeAttribute (thisnode, "swidth")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "stype")
         getstructgroundstate%stype = "Gaussian"
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "stype", &
           & getstructgroundstate%stype)
            Call removeAttribute (thisnode, "stype")
         End If
         getstructgroundstate%stypenumber = stringtonumberstype &
        & (getstructgroundstate%stype)
!
         Nullify (np)
         np => getAttributeNode (thisnode, "findlinentype")
         getstructgroundstate%findlinentype = "advanced"
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "findlinentype", &
           & getstructgroundstate%findlinentype)
            Call removeAttribute (thisnode, "findlinentype")
         End If
         getstructgroundstate%findlinentypenumber = &
        & stringtonumberfindlinentype &
        & (getstructgroundstate%findlinentype)
!
         Nullify (np)
         np => getAttributeNode (thisnode, "isgkmax")
         getstructgroundstate%isgkmax = - 1
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "isgkmax", &
           & getstructgroundstate%isgkmax)
            Call removeAttribute (thisnode, "isgkmax")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "gmaxvr")
         getstructgroundstate%gmaxvr = 12
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "gmaxvr", &
           & getstructgroundstate%gmaxvr)
            Call removeAttribute (thisnode, "gmaxvr")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "nempty")
         getstructgroundstate%nempty = 5
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "nempty", &
           & getstructgroundstate%nempty)
            Call removeAttribute (thisnode, "nempty")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "nosym")
         getstructgroundstate%nosym = .False.
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "nosym", &
           & getstructgroundstate%nosym)
            Call removeAttribute (thisnode, "nosym")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "frozencore")
         getstructgroundstate%frozencore = .False.
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "frozencore", &
           & getstructgroundstate%frozencore)
            Call removeAttribute (thisnode, "frozencore")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "autokpt")
         getstructgroundstate%autokpt = .False.
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "autokpt", &
           & getstructgroundstate%autokpt)
            Call removeAttribute (thisnode, "autokpt")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "radkpt")
         getstructgroundstate%radkpt = 40
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "radkpt", &
           & getstructgroundstate%radkpt)
            Call removeAttribute (thisnode, "radkpt")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "reducek")
         getstructgroundstate%reducek = .True.
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "reducek", &
           & getstructgroundstate%reducek)
            Call removeAttribute (thisnode, "reducek")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "tfibs")
         getstructgroundstate%tfibs = .True.
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "tfibs", &
           & getstructgroundstate%tfibs)
            Call removeAttribute (thisnode, "tfibs")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "tforce")
         getstructgroundstate%tforce = .False.
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "tforce", &
           & getstructgroundstate%tforce)
            Call removeAttribute (thisnode, "tforce")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "lmaxapw")
         getstructgroundstate%lmaxapw = 10
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "lmaxapw", &
           & getstructgroundstate%lmaxapw)
            Call removeAttribute (thisnode, "lmaxapw")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "maxscl")
         getstructgroundstate%maxscl = 200
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "maxscl", &
           & getstructgroundstate%maxscl)
            Call removeAttribute (thisnode, "maxscl")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "chgexs")
         getstructgroundstate%chgexs = 0
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "chgexs", &
           & getstructgroundstate%chgexs)
            Call removeAttribute (thisnode, "chgexs")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "deband")
         getstructgroundstate%deband = 0.0025d0
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "deband", &
           & getstructgroundstate%deband)
            Call removeAttribute (thisnode, "deband")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "epschg")
         getstructgroundstate%epschg = 1.0d-3
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "epschg", &
           & getstructgroundstate%epschg)
            Call removeAttribute (thisnode, "epschg")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "epsocc")
         getstructgroundstate%epsocc = 1e-8
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "epsocc", &
           & getstructgroundstate%epsocc)
            Call removeAttribute (thisnode, "epsocc")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "mixer")
         getstructgroundstate%mixer = "msec"
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "mixer", &
           & getstructgroundstate%mixer)
            Call removeAttribute (thisnode, "mixer")
         End If
         getstructgroundstate%mixernumber = stringtonumbermixer &
        & (getstructgroundstate%mixer)
!
         Nullify (np)
         np => getAttributeNode (thisnode, "beta0")
         getstructgroundstate%beta0 = 0.4
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "beta0", &
           & getstructgroundstate%beta0)
            Call removeAttribute (thisnode, "beta0")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "betainc")
         getstructgroundstate%betainc = 1.1
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "betainc", &
           & getstructgroundstate%betainc)
            Call removeAttribute (thisnode, "betainc")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "betadec")
         getstructgroundstate%betadec = 0.6
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "betadec", &
           & getstructgroundstate%betadec)
            Call removeAttribute (thisnode, "betadec")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "lradstep")
         getstructgroundstate%lradstep = 4
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "lradstep", &
           & getstructgroundstate%lradstep)
            Call removeAttribute (thisnode, "lradstep")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "nprad")
         getstructgroundstate%nprad = 4
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "nprad", &
           & getstructgroundstate%nprad)
            Call removeAttribute (thisnode, "nprad")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "xctype")
         getstructgroundstate%xctype = "LSDAPerdew-Wang"
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "xctype", &
           & getstructgroundstate%xctype)
            Call removeAttribute (thisnode, "xctype")
         End If
         getstructgroundstate%xctypenumber = stringtonumberxctype &
        & (getstructgroundstate%xctype)
!
         Nullify (np)
         np => getAttributeNode (thisnode, "evalmin")
         getstructgroundstate%evalmin = - 4.5d0
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "evalmin", &
           & getstructgroundstate%evalmin)
            Call removeAttribute (thisnode, "evalmin")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "lmaxvr")
         getstructgroundstate%lmaxvr = 6
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "lmaxvr", &
           & getstructgroundstate%lmaxvr)
            Call removeAttribute (thisnode, "lmaxvr")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "fracinr")
         getstructgroundstate%fracinr = 0.25d0
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "fracinr", &
           & getstructgroundstate%fracinr)
            Call removeAttribute (thisnode, "fracinr")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "lmaxinr")
         getstructgroundstate%lmaxinr = 2
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "lmaxinr", &
           & getstructgroundstate%lmaxinr)
            Call removeAttribute (thisnode, "lmaxinr")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "lmaxmat")
         getstructgroundstate%lmaxmat = 5
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "lmaxmat", &
           & getstructgroundstate%lmaxmat)
            Call removeAttribute (thisnode, "lmaxmat")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "vkloff")
         getstructgroundstate%vkloff = (/ 0, 0, 0 /)
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "vkloff", &
           & getstructgroundstate%vkloff)
            Call removeAttribute (thisnode, "vkloff")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "npsden")
         getstructgroundstate%npsden = 9
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "npsden", &
           & getstructgroundstate%npsden)
            Call removeAttribute (thisnode, "npsden")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "cfdamp")
         getstructgroundstate%cfdamp = 0
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "cfdamp", &
           & getstructgroundstate%cfdamp)
            Call removeAttribute (thisnode, "cfdamp")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "nosource")
         getstructgroundstate%nosource = .False.
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "nosource", &
           & getstructgroundstate%nosource)
            Call removeAttribute (thisnode, "nosource")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "tevecsv")
         getstructgroundstate%tevecsv = .False.
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "tevecsv", &
           & getstructgroundstate%tevecsv)
            Call removeAttribute (thisnode, "tevecsv")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "nwrite")
         getstructgroundstate%nwrite = 0
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "nwrite", &
           & getstructgroundstate%nwrite)
            Call removeAttribute (thisnode, "nwrite")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "ptnucl")
         getstructgroundstate%ptnucl = .True.
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "ptnucl", &
           & getstructgroundstate%ptnucl)
            Call removeAttribute (thisnode, "ptnucl")
         End If
!
         Len = countChildEmentsWithName (thisnode, "spin")
         getstructgroundstate%spin => null ()
         Do i = 0, len - 1
            getstructgroundstate%spin => getstructspin &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "spin"), 0)))
         End Do
!
         Len = countChildEmentsWithName (thisnode, "HartreeFock")
         getstructgroundstate%HartreeFock => null ()
         Do i = 0, len - 1
            getstructgroundstate%HartreeFock => getstructHartreeFock &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "HartreeFock"), 0)))
         End Do
!
         Len = countChildEmentsWithName (thisnode, "solver")
         getstructgroundstate%solver => null ()
         Do i = 0, len - 1
            getstructgroundstate%solver => getstructsolver &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "solver"), 0)))
         End Do
!
         Len = countChildEmentsWithName (thisnode, "OEP")
         getstructgroundstate%OEP => null ()
         Do i = 0, len - 1
            getstructgroundstate%OEP => getstructOEP &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "OEP"), 0)))
         End Do
!
         Len = countChildEmentsWithName (thisnode, "RDMFT")
         getstructgroundstate%RDMFT => null ()
         Do i = 0, len - 1
            getstructgroundstate%RDMFT => getstructRDMFT &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "RDMFT"), 0)))
         End Do
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructspin (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (spin_type), Pointer :: getstructspin
         Type (Node), Pointer :: np
!
!
         Integer :: Len = 1, i = 0
         Allocate (getstructspin)
#ifdef INPUTDEBUG
         Write (*,*) "we are at spin"
#endif
!
         Nullify (np)
         np => getAttributeNode (thisnode, "momfix")
         getstructspin%momfix = (/ 0, 0, 0 /)
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "momfix", &
           & getstructspin%momfix)
            Call removeAttribute (thisnode, "momfix")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "bfieldc")
         getstructspin%bfieldc = (/ 0, 0, 0 /)
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "bfieldc", &
           & getstructspin%bfieldc)
            Call removeAttribute (thisnode, "bfieldc")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "spinorb")
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "spinorb", &
           & getstructspin%spinorb)
            Call removeAttribute (thisnode, "spinorb")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "spinsprl")
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "spinsprl", &
           & getstructspin%spinsprl)
            Call removeAttribute (thisnode, "spinsprl")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "vqlss")
         getstructspin%vqlss = (/ 0, 0, 0 /)
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "vqlss", &
           & getstructspin%vqlss)
            Call removeAttribute (thisnode, "vqlss")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "taufsm")
         getstructspin%taufsm = 0.01d0
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "taufsm", &
           & getstructspin%taufsm)
            Call removeAttribute (thisnode, "taufsm")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "reducebf")
         getstructspin%reducebf = 1
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "reducebf", &
           & getstructspin%reducebf)
            Call removeAttribute (thisnode, "reducebf")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "fixspin")
         getstructspin%fixspin = "none"
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "fixspin", &
           & getstructspin%fixspin)
            Call removeAttribute (thisnode, "fixspin")
         End If
         getstructspin%fixspinnumber = stringtonumberfixspin &
        & (getstructspin%fixspin)
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructHartreeFock (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (HartreeFock_type), Pointer :: getstructHartreeFock
         Type (Node), Pointer :: np
!
!
         Integer :: Len = 1, i = 0
         Allocate (getstructHartreeFock)
#ifdef INPUTDEBUG
         Write (*,*) "we are at HartreeFock"
#endif
!
         Nullify (np)
         np => getAttributeNode (thisnode, "epsengy")
         getstructHartreeFock%epsengy = 1e-7
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "epsengy", &
           & getstructHartreeFock%epsengy)
            Call removeAttribute (thisnode, "epsengy")
         End If
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructsolver (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (solver_type), Pointer :: getstructsolver
         Type (Node), Pointer :: np
!
!
         Integer :: Len = 1, i = 0
         Allocate (getstructsolver)
#ifdef INPUTDEBUG
         Write (*,*) "we are at solver"
#endif
!
         Nullify (np)
         np => getAttributeNode (thisnode, "type")
         getstructsolver%type = "Lapack"
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "type", &
           & getstructsolver%type)
            Call removeAttribute (thisnode, "type")
         End If
         getstructsolver%typenumber = stringtonumbertype &
        & (getstructsolver%type)
!
         Nullify (np)
         np => getAttributeNode (thisnode, "packedmatrixstorage")
         getstructsolver%packedmatrixstorage = .True.
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "packedmatrixstorage", &
           & getstructsolver%packedmatrixstorage)
            Call removeAttribute (thisnode, "packedmatrixstorage")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "epsarpack")
         getstructsolver%epsarpack = 1.0e-8
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "epsarpack", &
           & getstructsolver%epsarpack)
            Call removeAttribute (thisnode, "epsarpack")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "evaltol")
         getstructsolver%evaltol = 1e-8
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "evaltol", &
           & getstructsolver%evaltol)
            Call removeAttribute (thisnode, "evaltol")
         End If
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructOEP (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (OEP_type), Pointer :: getstructOEP
         Type (Node), Pointer :: np
!
!
         Integer :: Len = 1, i = 0
         Allocate (getstructOEP)
#ifdef INPUTDEBUG
         Write (*,*) "we are at OEP"
#endif
!
         Nullify (np)
         np => getAttributeNode (thisnode, "maxitoep")
         getstructOEP%maxitoep = 120
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "maxitoep", &
           & getstructOEP%maxitoep)
            Call removeAttribute (thisnode, "maxitoep")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "tauoep")
         getstructOEP%tauoep = (/ 1.0, 0.2, 1.5 /)
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "tauoep", &
           & getstructOEP%tauoep)
            Call removeAttribute (thisnode, "tauoep")
         End If
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructRDMFT (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (RDMFT_type), Pointer :: getstructRDMFT
         Type (Node), Pointer :: np
!
!
         Integer :: Len = 1, i = 0
         Allocate (getstructRDMFT)
#ifdef INPUTDEBUG
         Write (*,*) "we are at RDMFT"
#endif
!
         Nullify (np)
         np => getAttributeNode (thisnode, "rdmxctype")
         getstructRDMFT%rdmxctype = 2
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "rdmxctype", &
           & getstructRDMFT%rdmxctype)
            Call removeAttribute (thisnode, "rdmxctype")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "rdmmaxscl")
         getstructRDMFT%rdmmaxscl = 1
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "rdmmaxscl", &
           & getstructRDMFT%rdmmaxscl)
            Call removeAttribute (thisnode, "rdmmaxscl")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "maxitn")
         getstructRDMFT%maxitn = 250
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "maxitn", &
           & getstructRDMFT%maxitn)
            Call removeAttribute (thisnode, "maxitn")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "maxitc")
         getstructRDMFT%maxitc = 10
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "maxitc", &
           & getstructRDMFT%maxitc)
            Call removeAttribute (thisnode, "maxitc")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "taurdmn")
         getstructRDMFT%taurdmn = 1.0
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "taurdmn", &
           & getstructRDMFT%taurdmn)
            Call removeAttribute (thisnode, "taurdmn")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "taurdmc")
         getstructRDMFT%taurdmc = 0.5
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "taurdmc", &
           & getstructRDMFT%taurdmc)
            Call removeAttribute (thisnode, "taurdmc")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "rdmalpha")
         getstructRDMFT%rdmalpha = 0.7
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "rdmalpha", &
           & getstructRDMFT%rdmalpha)
            Call removeAttribute (thisnode, "rdmalpha")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "rdmtemp")
         getstructRDMFT%rdmtemp = 0.0
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "rdmtemp", &
           & getstructRDMFT%rdmtemp)
            Call removeAttribute (thisnode, "rdmtemp")
         End If
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructstructureoptimization (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (structureoptimization_type), Pointer :: &
        & getstructstructureoptimization
         Type (Node), Pointer :: np
!
!
         Integer :: Len = 1, i = 0
         Allocate (getstructstructureoptimization)
#ifdef INPUTDEBUG
         Write (*,*) "we are at structureoptimization"
#endif
!
         Nullify (np)
         np => getAttributeNode (thisnode, "epsforce")
         getstructstructureoptimization%epsforce = 5e-5
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "epsforce", &
           & getstructstructureoptimization%epsforce)
            Call removeAttribute (thisnode, "epsforce")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "tau0atm")
         getstructstructureoptimization%tau0atm = 0.2d0
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "tau0atm", &
           & getstructstructureoptimization%tau0atm)
            Call removeAttribute (thisnode, "tau0atm")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "resume")
         getstructstructureoptimization%resume = .False.
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "resume", &
           & getstructstructureoptimization%resume)
            Call removeAttribute (thisnode, "resume")
         End If
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructproperties (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (properties_type), Pointer :: getstructproperties
		
         Integer :: Len = 1, i = 0
         Allocate (getstructproperties)
#ifdef INPUTDEBUG
         Write (*,*) "we are at properties"
#endif
!
         Len = countChildEmentsWithName (thisnode, "bandstructure")
         getstructproperties%bandstructure => null ()
         Do i = 0, len - 1
            getstructproperties%bandstructure => getstructbandstructure &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "bandstructure"), 0)))
         End Do
!
         Len = countChildEmentsWithName (thisnode, "STM")
         getstructproperties%STM => null ()
         Do i = 0, len - 1
            getstructproperties%STM => getstructSTM &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "STM"), 0)))
         End Do
!
         Len = countChildEmentsWithName (thisnode, "wfplot")
         getstructproperties%wfplot => null ()
         Do i = 0, len - 1
            getstructproperties%wfplot => getstructwfplot &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "wfplot"), 0)))
         End Do
!
         Len = countChildEmentsWithName (thisnode, "dos")
         getstructproperties%dos => null ()
         Do i = 0, len - 1
            getstructproperties%dos => getstructdos &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "dos"), 0)))
         End Do
!
         Len = countChildEmentsWithName (thisnode, "LSJ")
         getstructproperties%LSJ => null ()
         Do i = 0, len - 1
            getstructproperties%LSJ => getstructLSJ &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "LSJ"), 0)))
         End Do
!
         Len = countChildEmentsWithName (thisnode, "masstensor")
         getstructproperties%masstensor => null ()
         Do i = 0, len - 1
            getstructproperties%masstensor => getstructmasstensor &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "masstensor"), 0)))
         End Do
!
         Len = countChildEmentsWithName (thisnode, "chargedensityplot")
         getstructproperties%chargedensityplot => null ()
         Do i = 0, len - 1
            getstructproperties%chargedensityplot => &
           & getstructchargedensityplot (removeChild(thisnode, &
           & item(getElementsByTagname(thisnode, "chargedensityplot"), &
           & 0)))
         End Do
!
         Len = countChildEmentsWithName (thisnode, "exccplot")
         getstructproperties%exccplot => null ()
         Do i = 0, len - 1
            getstructproperties%exccplot => getstructexccplot &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "exccplot"), 0)))
         End Do
!
         Len = countChildEmentsWithName (thisnode, "elfplot")
         getstructproperties%elfplot => null ()
         Do i = 0, len - 1
            getstructproperties%elfplot => getstructelfplot &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "elfplot"), 0)))
         End Do
!
         Len = countChildEmentsWithName (thisnode, "mvecfield")
         getstructproperties%mvecfield => null ()
         Do i = 0, len - 1
            getstructproperties%mvecfield => getstructmvecfield &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "mvecfield"), 0)))
         End Do
!
         Len = countChildEmentsWithName (thisnode, "xcmvecfield")
         getstructproperties%xcmvecfield => null ()
         Do i = 0, len - 1
            getstructproperties%xcmvecfield => getstructxcmvecfield &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "xcmvecfield"), 0)))
         End Do
!
         Len = countChildEmentsWithName (thisnode, "electricfield")
         getstructproperties%electricfield => null ()
         Do i = 0, len - 1
            getstructproperties%electricfield => getstructelectricfield &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "electricfield"), 0)))
         End Do
!
         Len = countChildEmentsWithName (thisnode, "gradmvecfield")
         getstructproperties%gradmvecfield => null ()
         Do i = 0, len - 1
            getstructproperties%gradmvecfield => getstructgradmvecfield &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "gradmvecfield"), 0)))
         End Do
!
         Len = countChildEmentsWithName (thisnode, "fermisurfaceplot")
         getstructproperties%fermisurfaceplot => null ()
         Do i = 0, len - 1
            getstructproperties%fermisurfaceplot => &
           & getstructfermisurfaceplot (removeChild(thisnode, &
           & item(getElementsByTagname(thisnode, "fermisurfaceplot"), &
           & 0)))
         End Do
!
         Len = countChildEmentsWithName (thisnode, "EFG")
         getstructproperties%EFG => null ()
         Do i = 0, len - 1
            getstructproperties%EFG => getstructEFG &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "EFG"), 0)))
         End Do
!
         Len = countChildEmentsWithName (thisnode, "momentummatrix")
         getstructproperties%momentummatrix => null ()
         Do i = 0, len - 1
            getstructproperties%momentummatrix => &
           & getstructmomentummatrix (removeChild(thisnode, &
           & item(getElementsByTagname(thisnode, "momentummatrix"), &
           & 0)))
         End Do
!
         Len = countChildEmentsWithName (thisnode, "linresponsetensor")
         getstructproperties%linresponsetensor => null ()
         Do i = 0, len - 1
            getstructproperties%linresponsetensor => &
           & getstructlinresponsetensor (removeChild(thisnode, &
           & item(getElementsByTagname(thisnode, "linresponsetensor"), &
           & 0)))
         End Do
!
         Len = countChildEmentsWithName (thisnode, "mossbauer")
         getstructproperties%mossbauer => null ()
         Do i = 0, len - 1
            getstructproperties%mossbauer => getstructmossbauer &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "mossbauer"), 0)))
         End Do
!
         Len = countChildEmentsWithName (thisnode, "dielectric")
         getstructproperties%dielectric => null ()
         Do i = 0, len - 1
            getstructproperties%dielectric => getstructdielectric &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "dielectric"), 0)))
         End Do
!
         Len = countChildEmentsWithName (thisnode, "expiqr")
         getstructproperties%expiqr => null ()
         Do i = 0, len - 1
            getstructproperties%expiqr => getstructexpiqr &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "expiqr"), 0)))
         End Do
!
         Len = countChildEmentsWithName (thisnode, "elnes")
         getstructproperties%elnes => null ()
         Do i = 0, len - 1
            getstructproperties%elnes => getstructelnes &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "elnes"), 0)))
         End Do
!
         Len = countChildEmentsWithName (thisnode, "eliashberg")
         getstructproperties%eliashberg => null ()
         Do i = 0, len - 1
            getstructproperties%eliashberg => getstructeliashberg &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "eliashberg"), 0)))
         End Do
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructbandstructure (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (bandstructure_type), Pointer :: getstructbandstructure
         Type (Node), Pointer :: np
!
!
         Integer :: Len = 1, i = 0
         Allocate (getstructbandstructure)
#ifdef INPUTDEBUG
         Write (*,*) "we are at bandstructure"
#endif
!
         Nullify (np)
         np => getAttributeNode (thisnode, "scissor")
         getstructbandstructure%scissor = 0d0
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "scissor", &
           & getstructbandstructure%scissor)
            Call removeAttribute (thisnode, "scissor")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "character")
         getstructbandstructure%character = .False.
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "character", &
           & getstructbandstructure%character)
            Call removeAttribute (thisnode, "character")
         End If
!
         Len = countChildEmentsWithName (thisnode, "plot1d")
         getstructbandstructure%plot1d => null ()
         Do i = 0, len - 1
            getstructbandstructure%plot1d => getstructplot1d &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "plot1d"), 0)))
         End Do
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructSTM (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (STM_type), Pointer :: getstructSTM
		
         Integer :: Len = 1, i = 0
         Allocate (getstructSTM)
#ifdef INPUTDEBUG
         Write (*,*) "we are at STM"
#endif
!
         Len = countChildEmentsWithName (thisnode, "plot2d")
         getstructSTM%plot2d => null ()
         Do i = 0, len - 1
            getstructSTM%plot2d => getstructplot2d &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "plot2d"), 0)))
         End Do
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructwfplot (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (wfplot_type), Pointer :: getstructwfplot
		
         Integer :: Len = 1, i = 0
         Allocate (getstructwfplot)
#ifdef INPUTDEBUG
         Write (*,*) "we are at wfplot"
#endif
!
         Len = countChildEmentsWithName (thisnode, "kstlist")
         getstructwfplot%kstlist => null ()
         Do i = 0, len - 1
            getstructwfplot%kstlist => getstructkstlist &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "kstlist"), 0)))
         End Do
!
         Len = countChildEmentsWithName (thisnode, "plot1d")
         getstructwfplot%plot1d => null ()
         Do i = 0, len - 1
            getstructwfplot%plot1d => getstructplot1d &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "plot1d"), 0)))
         End Do
!
         Len = countChildEmentsWithName (thisnode, "plot2d")
         getstructwfplot%plot2d => null ()
         Do i = 0, len - 1
            getstructwfplot%plot2d => getstructplot2d &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "plot2d"), 0)))
         End Do
!
         Len = countChildEmentsWithName (thisnode, "plot3d")
         getstructwfplot%plot3d => null ()
         Do i = 0, len - 1
            getstructwfplot%plot3d => getstructplot3d &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "plot3d"), 0)))
         End Do
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructdos (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (dos_type), Pointer :: getstructdos
         Type (Node), Pointer :: np
!
!
         Integer :: Len = 1, i = 0
         Allocate (getstructdos)
#ifdef INPUTDEBUG
         Write (*,*) "we are at dos"
#endif
!
         Nullify (np)
         np => getAttributeNode (thisnode, "sqados")
         getstructdos%sqados = (/ 0.0, 0.0, 1.0 /)
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "sqados", &
           & getstructdos%sqados)
            Call removeAttribute (thisnode, "sqados")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "lmirep")
         getstructdos%lmirep = .False.
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "lmirep", &
           & getstructdos%lmirep)
            Call removeAttribute (thisnode, "lmirep")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "nwdos")
         getstructdos%nwdos = 500
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "nwdos", &
           & getstructdos%nwdos)
            Call removeAttribute (thisnode, "nwdos")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "ngrdos")
         getstructdos%ngrdos = 100
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "ngrdos", &
           & getstructdos%ngrdos)
            Call removeAttribute (thisnode, "ngrdos")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "scissor")
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "scissor", &
           & getstructdos%scissor)
            Call removeAttribute (thisnode, "scissor")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "nsmdos")
         getstructdos%nsmdos = 0
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "nsmdos", &
           & getstructdos%nsmdos)
            Call removeAttribute (thisnode, "nsmdos")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "winddos")
         getstructdos%winddos = (/ .5, .5 /)
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "winddos", &
           & getstructdos%winddos)
            Call removeAttribute (thisnode, "winddos")
         End If
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructLSJ (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (LSJ_type), Pointer :: getstructLSJ
		
         Integer :: Len = 1, i = 0
         Allocate (getstructLSJ)
#ifdef INPUTDEBUG
         Write (*,*) "we are at LSJ"
#endif
!
         Len = countChildEmentsWithName (thisnode, "kstlist")
         getstructLSJ%kstlist => null ()
         Do i = 0, len - 1
            getstructLSJ%kstlist => getstructkstlist &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "kstlist"), 0)))
         End Do
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructmasstensor (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (masstensor_type), Pointer :: getstructmasstensor
         Type (Node), Pointer :: np
!
!
         Integer :: Len = 1, i = 0
         Allocate (getstructmasstensor)
#ifdef INPUTDEBUG
         Write (*,*) "we are at masstensor"
#endif
!
         Nullify (np)
         np => getAttributeNode (thisnode, "deltaem")
         getstructmasstensor%deltaem = 0.025d0
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "deltaem", &
           & getstructmasstensor%deltaem)
            Call removeAttribute (thisnode, "deltaem")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "ndspem")
         getstructmasstensor%ndspem = 1
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "ndspem", &
           & getstructmasstensor%ndspem)
            Call removeAttribute (thisnode, "ndspem")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "vklem")
         getstructmasstensor%vklem = (/ 0, 0, 0 /)
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "vklem", &
           & getstructmasstensor%vklem)
            Call removeAttribute (thisnode, "vklem")
         End If
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructchargedensityplot (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (chargedensityplot_type), Pointer :: &
        & getstructchargedensityplot
		
         Integer :: Len = 1, i = 0
         Allocate (getstructchargedensityplot)
#ifdef INPUTDEBUG
         Write (*,*) "we are at chargedensityplot"
#endif
!
         Len = countChildEmentsWithName (thisnode, "plot1d")
         getstructchargedensityplot%plot1d => null ()
         Do i = 0, len - 1
            getstructchargedensityplot%plot1d => getstructplot1d &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "plot1d"), 0)))
         End Do
!
         Len = countChildEmentsWithName (thisnode, "plot2d")
         getstructchargedensityplot%plot2d => null ()
         Do i = 0, len - 1
            getstructchargedensityplot%plot2d => getstructplot2d &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "plot2d"), 0)))
         End Do
!
         Len = countChildEmentsWithName (thisnode, "plot3d")
         getstructchargedensityplot%plot3d => null ()
         Do i = 0, len - 1
            getstructchargedensityplot%plot3d => getstructplot3d &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "plot3d"), 0)))
         End Do
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructexccplot (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (exccplot_type), Pointer :: getstructexccplot
		
         Integer :: Len = 1, i = 0
         Allocate (getstructexccplot)
#ifdef INPUTDEBUG
         Write (*,*) "we are at exccplot"
#endif
!
         Len = countChildEmentsWithName (thisnode, "plot1d")
         getstructexccplot%plot1d => null ()
         Do i = 0, len - 1
            getstructexccplot%plot1d => getstructplot1d &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "plot1d"), 0)))
         End Do
!
         Len = countChildEmentsWithName (thisnode, "plot2d")
         getstructexccplot%plot2d => null ()
         Do i = 0, len - 1
            getstructexccplot%plot2d => getstructplot2d &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "plot2d"), 0)))
         End Do
!
         Len = countChildEmentsWithName (thisnode, "plot3d")
         getstructexccplot%plot3d => null ()
         Do i = 0, len - 1
            getstructexccplot%plot3d => getstructplot3d &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "plot3d"), 0)))
         End Do
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructelfplot (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (elfplot_type), Pointer :: getstructelfplot
		
         Integer :: Len = 1, i = 0
         Allocate (getstructelfplot)
#ifdef INPUTDEBUG
         Write (*,*) "we are at elfplot"
#endif
!
         Len = countChildEmentsWithName (thisnode, "plot1d")
         getstructelfplot%plot1d => null ()
         Do i = 0, len - 1
            getstructelfplot%plot1d => getstructplot1d &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "plot1d"), 0)))
         End Do
!
         Len = countChildEmentsWithName (thisnode, "plot2d")
         getstructelfplot%plot2d => null ()
         Do i = 0, len - 1
            getstructelfplot%plot2d => getstructplot2d &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "plot2d"), 0)))
         End Do
!
         Len = countChildEmentsWithName (thisnode, "plot3d")
         getstructelfplot%plot3d => null ()
         Do i = 0, len - 1
            getstructelfplot%plot3d => getstructplot3d &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "plot3d"), 0)))
         End Do
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructmvecfield (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (mvecfield_type), Pointer :: getstructmvecfield
		
         Integer :: Len = 1, i = 0
         Allocate (getstructmvecfield)
#ifdef INPUTDEBUG
         Write (*,*) "we are at mvecfield"
#endif
!
         Len = countChildEmentsWithName (thisnode, "plot2d")
         getstructmvecfield%plot2d => null ()
         Do i = 0, len - 1
            getstructmvecfield%plot2d => getstructplot2d &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "plot2d"), 0)))
         End Do
!
         Len = countChildEmentsWithName (thisnode, "plot3d")
         getstructmvecfield%plot3d => null ()
         Do i = 0, len - 1
            getstructmvecfield%plot3d => getstructplot3d &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "plot3d"), 0)))
         End Do
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructxcmvecfield (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (xcmvecfield_type), Pointer :: getstructxcmvecfield
		
         Integer :: Len = 1, i = 0
         Allocate (getstructxcmvecfield)
#ifdef INPUTDEBUG
         Write (*,*) "we are at xcmvecfield"
#endif
!
         Len = countChildEmentsWithName (thisnode, "plot2d")
         getstructxcmvecfield%plot2d => null ()
         Do i = 0, len - 1
            getstructxcmvecfield%plot2d => getstructplot2d &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "plot2d"), 0)))
         End Do
!
         Len = countChildEmentsWithName (thisnode, "plot3d")
         getstructxcmvecfield%plot3d => null ()
         Do i = 0, len - 1
            getstructxcmvecfield%plot3d => getstructplot3d &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "plot3d"), 0)))
         End Do
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructelectricfield (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (electricfield_type), Pointer :: getstructelectricfield
		
         Integer :: Len = 1, i = 0
         Allocate (getstructelectricfield)
#ifdef INPUTDEBUG
         Write (*,*) "we are at electricfield"
#endif
!
         Len = countChildEmentsWithName (thisnode, "plot2d")
         getstructelectricfield%plot2d => null ()
         Do i = 0, len - 1
            getstructelectricfield%plot2d => getstructplot2d &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "plot2d"), 0)))
         End Do
!
         Len = countChildEmentsWithName (thisnode, "plot3d")
         getstructelectricfield%plot3d => null ()
         Do i = 0, len - 1
            getstructelectricfield%plot3d => getstructplot3d &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "plot3d"), 0)))
         End Do
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructgradmvecfield (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (gradmvecfield_type), Pointer :: getstructgradmvecfield
		
         Integer :: Len = 1, i = 0
         Allocate (getstructgradmvecfield)
#ifdef INPUTDEBUG
         Write (*,*) "we are at gradmvecfield"
#endif
!
         Len = countChildEmentsWithName (thisnode, "plot1d")
         getstructgradmvecfield%plot1d => null ()
         Do i = 0, len - 1
            getstructgradmvecfield%plot1d => getstructplot1d &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "plot1d"), 0)))
         End Do
!
         Len = countChildEmentsWithName (thisnode, "plot2d")
         getstructgradmvecfield%plot2d => null ()
         Do i = 0, len - 1
            getstructgradmvecfield%plot2d => getstructplot2d &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "plot2d"), 0)))
         End Do
!
         Len = countChildEmentsWithName (thisnode, "plot3d")
         getstructgradmvecfield%plot3d => null ()
         Do i = 0, len - 1
            getstructgradmvecfield%plot3d => getstructplot3d &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "plot3d"), 0)))
         End Do
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructfermisurfaceplot (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (fermisurfaceplot_type), Pointer :: &
        & getstructfermisurfaceplot
         Type (Node), Pointer :: np
!
!
         Integer :: Len = 1, i = 0
         Allocate (getstructfermisurfaceplot)
#ifdef INPUTDEBUG
         Write (*,*) "we are at fermisurfaceplot"
#endif
!
         Nullify (np)
         np => getAttributeNode (thisnode, "nstfsp")
         getstructfermisurfaceplot%nstfsp = 6
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "nstfsp", &
           & getstructfermisurfaceplot%nstfsp)
            Call removeAttribute (thisnode, "nstfsp")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "separate")
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "separate", &
           & getstructfermisurfaceplot%separate)
            Call removeAttribute (thisnode, "separate")
         End If
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructEFG (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (EFG_type), Pointer :: getstructEFG
		
         Integer :: Len = 1, i = 0
         Allocate (getstructEFG)
#ifdef INPUTDEBUG
         Write (*,*) "we are at EFG"
#endif
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructmomentummatrix (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (momentummatrix_type), Pointer :: getstructmomentummatrix
		
         Integer :: Len = 1, i = 0
         Allocate (getstructmomentummatrix)
#ifdef INPUTDEBUG
         Write (*,*) "we are at momentummatrix"
#endif
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructlinresponsetensor (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (linresponsetensor_type), Pointer :: &
        & getstructlinresponsetensor
         Type (Node), Pointer :: np
!
!
         Integer :: Len = 1, i = 0
         Allocate (getstructlinresponsetensor)
#ifdef INPUTDEBUG
         Write (*,*) "we are at linresponsetensor"
#endif
!
         Nullify (np)
         np => getAttributeNode (thisnode, "scissor")
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "scissor", &
           & getstructlinresponsetensor%scissor)
            Call removeAttribute (thisnode, "scissor")
         End If
!
         Len = countChildEmentsWithName (thisnode, "optcomp")
         Allocate (getstructlinresponsetensor%optcomp(3, len))
         Do i = 1, len
!
            getstructlinresponsetensor%optcomp (:, i) = &
           & getvalueofoptcomp (removeChild(thisnode, &
           & item(getElementsByTagname(thisnode, "optcomp"), 0)))
         End Do
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructmossbauer (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (mossbauer_type), Pointer :: getstructmossbauer
		
         Integer :: Len = 1, i = 0
         Allocate (getstructmossbauer)
#ifdef INPUTDEBUG
         Write (*,*) "we are at mossbauer"
#endif
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructdielectric (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (dielectric_type), Pointer :: getstructdielectric
		
         Integer :: Len = 1, i = 0
         Allocate (getstructdielectric)
#ifdef INPUTDEBUG
         Write (*,*) "we are at dielectric"
#endif
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructexpiqr (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (expiqr_type), Pointer :: getstructexpiqr
		
         Integer :: Len = 1, i = 0
         Allocate (getstructexpiqr)
#ifdef INPUTDEBUG
         Write (*,*) "we are at expiqr"
#endif
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructelnes (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (elnes_type), Pointer :: getstructelnes
         Type (Node), Pointer :: np
!
!
         Integer :: Len = 1, i = 0
         Allocate (getstructelnes)
#ifdef INPUTDEBUG
         Write (*,*) "we are at elnes"
#endif
!
         Nullify (np)
         np => getAttributeNode (thisnode, "vecql")
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "vecql", &
           & getstructelnes%vecql)
            Call removeAttribute (thisnode, "vecql")
         End If
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructeliashberg (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (eliashberg_type), Pointer :: getstructeliashberg
         Type (Node), Pointer :: np
!
!
         Integer :: Len = 1, i = 0
         Allocate (getstructeliashberg)
#ifdef INPUTDEBUG
         Write (*,*) "we are at eliashberg"
#endif
!
         Nullify (np)
         np => getAttributeNode (thisnode, "mustar")
         getstructeliashberg%mustar = 0.15
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "mustar", &
           & getstructeliashberg%mustar)
            Call removeAttribute (thisnode, "mustar")
         End If
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructphonons (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (phonons_type), Pointer :: getstructphonons
         Type (Node), Pointer :: np
!
!
         Integer :: Len = 1, i = 0
         Allocate (getstructphonons)
#ifdef INPUTDEBUG
         Write (*,*) "we are at phonons"
#endif
!
         Nullify (np)
         np => getAttributeNode (thisnode, "reduceq")
         getstructphonons%reduceq = .True.
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "reduceq", &
           & getstructphonons%reduceq)
            Call removeAttribute (thisnode, "reduceq")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "deltaph")
         getstructphonons%deltaph = 0.03
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "deltaph", &
           & getstructphonons%deltaph)
            Call removeAttribute (thisnode, "deltaph")
         End If
!
         Len = countChildEmentsWithName (thisnode, "qpointset")
         getstructphonons%qpointset => null ()
         Do i = 0, len - 1
            getstructphonons%qpointset => getstructqpointset &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "qpointset"), 0)))
         End Do
!
         Len = countChildEmentsWithName (thisnode, "phonondos")
         getstructphonons%phonondos => null ()
         Do i = 0, len - 1
            getstructphonons%phonondos => getstructphonondos &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "phonondos"), 0)))
         End Do
!
         Len = countChildEmentsWithName (thisnode, "phonondispplot")
         getstructphonons%phonondispplot => null ()
         Do i = 0, len - 1
            getstructphonons%phonondispplot => getstructphonondispplot &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "phonondispplot"), 0)))
         End Do
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructphonondos (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (phonondos_type), Pointer :: getstructphonondos
		
         Integer :: Len = 1, i = 0
         Allocate (getstructphonondos)
#ifdef INPUTDEBUG
         Write (*,*) "we are at phonondos"
#endif
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructphonondispplot (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (phonondispplot_type), Pointer :: getstructphonondispplot
		
         Integer :: Len = 1, i = 0
         Allocate (getstructphonondispplot)
#ifdef INPUTDEBUG
         Write (*,*) "we are at phonondispplot"
#endif
!
         Len = countChildEmentsWithName (thisnode, "plot1d")
         getstructphonondispplot%plot1d => null ()
         Do i = 0, len - 1
            getstructphonondispplot%plot1d => getstructplot1d &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "plot1d"), 0)))
         End Do
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructxs (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (xs_type), Pointer :: getstructxs
         Type (Node), Pointer :: np
!
!
         Integer :: Len = 1, i = 0
         Allocate (getstructxs)
#ifdef INPUTDEBUG
         Write (*,*) "we are at xs"
#endif
!
         Nullify (np)
         np => getAttributeNode (thisnode, "emattype")
         getstructxs%emattype = 1
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "emattype", &
           & getstructxs%emattype)
            Call removeAttribute (thisnode, "emattype")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "dfoffdiag")
         getstructxs%dfoffdiag = .False.
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "dfoffdiag", &
           & getstructxs%dfoffdiag)
            Call removeAttribute (thisnode, "dfoffdiag")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "lmaxapwwf")
         getstructxs%lmaxapwwf = - 1
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "lmaxapwwf", &
           & getstructxs%lmaxapwwf)
            Call removeAttribute (thisnode, "lmaxapwwf")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "lmaxemat")
         getstructxs%lmaxemat = 3
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "lmaxemat", &
           & getstructxs%lmaxemat)
            Call removeAttribute (thisnode, "lmaxemat")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "emaxdf")
         getstructxs%emaxdf = 1d10
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "emaxdf", &
           & getstructxs%emaxdf)
            Call removeAttribute (thisnode, "emaxdf")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "broad")
         getstructxs%broad = 0.01d0
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "broad", &
           & getstructxs%broad)
            Call removeAttribute (thisnode, "broad")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "tevout")
         getstructxs%tevout = .False.
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "tevout", &
           & getstructxs%tevout)
            Call removeAttribute (thisnode, "tevout")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "xstype")
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "xstype", &
           & getstructxs%xstype)
            Call removeAttribute (thisnode, "xstype")
         End If
         getstructxs%xstypenumber = stringtonumberxstype &
        & (getstructxs%xstype)
!
         Nullify (np)
         np => getAttributeNode (thisnode, "symmorph")
         getstructxs%symmorph = .False.
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "symmorph", &
           & getstructxs%symmorph)
            Call removeAttribute (thisnode, "symmorph")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "fastpmat")
         getstructxs%fastpmat = .True.
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "fastpmat", &
           & getstructxs%fastpmat)
            Call removeAttribute (thisnode, "fastpmat")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "fastemat")
         getstructxs%fastemat = .True.
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "fastemat", &
           & getstructxs%fastemat)
            Call removeAttribute (thisnode, "fastemat")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "gather")
         getstructxs%gather = .False.
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "gather", &
           & getstructxs%gather)
            Call removeAttribute (thisnode, "gather")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "tappinfo")
         getstructxs%tappinfo = .False.
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "tappinfo", &
           & getstructxs%tappinfo)
            Call removeAttribute (thisnode, "tappinfo")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "dbglev")
         getstructxs%dbglev = 0
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "dbglev", &
           & getstructxs%dbglev)
            Call removeAttribute (thisnode, "dbglev")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "usegdft")
         getstructxs%usegdft = .False.
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "usegdft", &
           & getstructxs%usegdft)
            Call removeAttribute (thisnode, "usegdft")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "gqmax")
         getstructxs%gqmax = 0
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "gqmax", &
           & getstructxs%gqmax)
            Call removeAttribute (thisnode, "gqmax")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "nosym")
         getstructxs%nosym = .False.
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "nosym", &
           & getstructxs%nosym)
            Call removeAttribute (thisnode, "nosym")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "ngridk")
         getstructxs%ngridk = (/ 1, 1, 1 /)
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "ngridk", &
           & getstructxs%ngridk)
            Call removeAttribute (thisnode, "ngridk")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "vkloff")
         getstructxs%vkloff = (/ 0, 0, 0 /)
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "vkloff", &
           & getstructxs%vkloff)
            Call removeAttribute (thisnode, "vkloff")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "reducek")
         getstructxs%reducek = .True.
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "reducek", &
           & getstructxs%reducek)
            Call removeAttribute (thisnode, "reducek")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "ngridq")
         getstructxs%ngridq = (/ 1, 1, 1 /)
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "ngridq", &
           & getstructxs%ngridq)
            Call removeAttribute (thisnode, "ngridq")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "reduceq")
         getstructxs%reduceq = .True.
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "reduceq", &
           & getstructxs%reduceq)
            Call removeAttribute (thisnode, "reduceq")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "rgkmax")
         getstructxs%rgkmax = 7
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "rgkmax", &
           & getstructxs%rgkmax)
            Call removeAttribute (thisnode, "rgkmax")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "swidth")
         getstructxs%swidth = 0.001d0
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "swidth", &
           & getstructxs%swidth)
            Call removeAttribute (thisnode, "swidth")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "lmaxapw")
         getstructxs%lmaxapw = 10
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "lmaxapw", &
           & getstructxs%lmaxapw)
            Call removeAttribute (thisnode, "lmaxapw")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "lmaxmat")
         getstructxs%lmaxmat = 5
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "lmaxmat", &
           & getstructxs%lmaxmat)
            Call removeAttribute (thisnode, "lmaxmat")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "nempty")
         getstructxs%nempty = 5
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "nempty", &
           & getstructxs%nempty)
            Call removeAttribute (thisnode, "nempty")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "scissor")
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "scissor", &
           & getstructxs%scissor)
            Call removeAttribute (thisnode, "scissor")
         End If
!
         Len = countChildEmentsWithName (thisnode, "tddft")
         getstructxs%tddft => null ()
         Do i = 0, len - 1
            getstructxs%tddft => getstructtddft (removeChild(thisnode, &
           & item(getElementsByTagname(thisnode, "tddft"), 0)))
         End Do
!
         Len = countChildEmentsWithName (thisnode, "screening")
         getstructxs%screening => null ()
         Do i = 0, len - 1
            getstructxs%screening => getstructscreening &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "screening"), 0)))
         End Do
!
         Len = countChildEmentsWithName (thisnode, "BSE")
         getstructxs%BSE => null ()
         Do i = 0, len - 1
            getstructxs%BSE => getstructBSE (removeChild(thisnode, &
           & item(getElementsByTagname(thisnode, "BSE"), 0)))
         End Do
!
         Len = countChildEmentsWithName (thisnode, "qpointset")
         getstructxs%qpointset => null ()
         Do i = 0, len - 1
            getstructxs%qpointset => getstructqpointset &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "qpointset"), 0)))
         End Do
!
         Len = countChildEmentsWithName (thisnode, "tetra")
         getstructxs%tetra => null ()
         Do i = 0, len - 1
            getstructxs%tetra => getstructtetra (removeChild(thisnode, &
           & item(getElementsByTagname(thisnode, "tetra"), 0)))
         End Do
!
         Len = countChildEmentsWithName (thisnode, "dosWindow")
         getstructxs%dosWindow => null ()
         Do i = 0, len - 1
            getstructxs%dosWindow => getstructdosWindow &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "dosWindow"), 0)))
         End Do
!
         Len = countChildEmentsWithName (thisnode, "plan")
         getstructxs%plan => null ()
         Do i = 0, len - 1
            getstructxs%plan => getstructplan (removeChild(thisnode, &
           & item(getElementsByTagname(thisnode, "plan"), 0)))
         End Do
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructtddft (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (tddft_type), Pointer :: getstructtddft
         Type (Node), Pointer :: np
!
!
         Integer :: Len = 1, i = 0
         Allocate (getstructtddft)
#ifdef INPUTDEBUG
         Write (*,*) "we are at tddft"
#endif
!
         Nullify (np)
         np => getAttributeNode (thisnode, "intraband")
         getstructtddft%intraband = .False.
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "intraband", &
           & getstructtddft%intraband)
            Call removeAttribute (thisnode, "intraband")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "torddf")
         getstructtddft%torddf = .False.
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "torddf", &
           & getstructtddft%torddf)
            Call removeAttribute (thisnode, "torddf")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "tordfxc")
         getstructtddft%tordfxc = .False.
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "tordfxc", &
           & getstructtddft%tordfxc)
            Call removeAttribute (thisnode, "tordfxc")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "aresdf")
         getstructtddft%aresdf = .True.
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "aresdf", &
           & getstructtddft%aresdf)
            Call removeAttribute (thisnode, "aresdf")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "aresfxc")
         getstructtddft%aresfxc = .True.
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "aresfxc", &
           & getstructtddft%aresfxc)
            Call removeAttribute (thisnode, "aresfxc")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "fxcbsesplit")
         getstructtddft%fxcbsesplit = 1d-5
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "fxcbsesplit", &
           & getstructtddft%fxcbsesplit)
            Call removeAttribute (thisnode, "fxcbsesplit")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "acont")
         getstructtddft%acont = .False.
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "acont", &
           & getstructtddft%acont)
            Call removeAttribute (thisnode, "acont")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "nwacont")
         getstructtddft%nwacont = 0
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "nwacont", &
           & getstructtddft%nwacont)
            Call removeAttribute (thisnode, "nwacont")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "lindhard")
         getstructtddft%lindhard = .False.
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "lindhard", &
           & getstructtddft%lindhard)
            Call removeAttribute (thisnode, "lindhard")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "epsdfde")
         getstructtddft%epsdfde = 1.0d-8
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "epsdfde", &
           & getstructtddft%epsdfde)
            Call removeAttribute (thisnode, "epsdfde")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "kerndiag")
         getstructtddft%kerndiag = .False.
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "kerndiag", &
           & getstructtddft%kerndiag)
            Call removeAttribute (thisnode, "kerndiag")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "lmaxalda")
         getstructtddft%lmaxalda = 3
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "lmaxalda", &
           & getstructtddft%lmaxalda)
            Call removeAttribute (thisnode, "lmaxalda")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "alphalrc")
         getstructtddft%alphalrc = 0
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "alphalrc", &
           & getstructtddft%alphalrc)
            Call removeAttribute (thisnode, "alphalrc")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "alphalrcdyn")
         getstructtddft%alphalrcdyn = 0
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "alphalrcdyn", &
           & getstructtddft%alphalrcdyn)
            Call removeAttribute (thisnode, "alphalrcdyn")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "betalrcdyn")
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "betalrcdyn", &
           & getstructtddft%betalrcdyn)
            Call removeAttribute (thisnode, "betalrcdyn")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "mdfqtype")
         getstructtddft%mdfqtype = 0
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "mdfqtype", &
           & getstructtddft%mdfqtype)
            Call removeAttribute (thisnode, "mdfqtype")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "fxctype")
         getstructtddft%fxctype = "RPA"
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "fxctype", &
           & getstructtddft%fxctype)
            Call removeAttribute (thisnode, "fxctype")
         End If
         getstructtddft%fxctypenumber = stringtonumberfxctype &
        & (getstructtddft%fxctype)
!
         Nullify (np)
         np => getAttributeNode (thisnode, "resumefromkernel")
         getstructtddft%resumefromkernel = .False.
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "resumefromkernel", &
           & getstructtddft%resumefromkernel)
            Call removeAttribute (thisnode, "resumefromkernel")
         End If
!
         Len = countChildEmentsWithName (thisnode, "dftrans")
         getstructtddft%dftrans => null ()
         Do i = 0, len - 1
            getstructtddft%dftrans => getstructdftrans &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "dftrans"), 0)))
         End Do
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructdftrans (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (dftrans_type), Pointer :: getstructdftrans
		
         Integer :: Len = 1, i = 0
         Allocate (getstructdftrans)
#ifdef INPUTDEBUG
         Write (*,*) "we are at dftrans"
#endif
!
         Len = countChildEmentsWithName (thisnode, "trans")
         Allocate (getstructdftrans%trans(3, len))
         Do i = 1, len
!
            getstructdftrans%trans (:, i) = getvalueoftrans &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "trans"), 0)))
         End Do
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructscreening (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (screening_type), Pointer :: getstructscreening
         Type (Node), Pointer :: np
!
!
         Integer :: Len = 1, i = 0
         Allocate (getstructscreening)
#ifdef INPUTDEBUG
         Write (*,*) "we are at screening"
#endif
!
         Nullify (np)
         np => getAttributeNode (thisnode, "run")
         getstructscreening%run = "fromscratch"
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "run", &
           & getstructscreening%run)
            Call removeAttribute (thisnode, "run")
         End If
         getstructscreening%runnumber = stringtonumberrun &
        & (getstructscreening%run)
!
         Nullify (np)
         np => getAttributeNode (thisnode, "nosym")
         getstructscreening%nosym = .False.
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "nosym", &
           & getstructscreening%nosym)
            Call removeAttribute (thisnode, "nosym")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "ngridk")
         getstructscreening%ngridk = (/ 0, 0, 0 /)
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "ngridk", &
           & getstructscreening%ngridk)
            Call removeAttribute (thisnode, "ngridk")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "reducek")
         getstructscreening%reducek = .False.
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "reducek", &
           & getstructscreening%reducek)
            Call removeAttribute (thisnode, "reducek")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "vkloff")
         getstructscreening%vkloff = (/ - 1, - 1, - 1 /)
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "vkloff", &
           & getstructscreening%vkloff)
            Call removeAttribute (thisnode, "vkloff")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "rgkmax")
         getstructscreening%rgkmax = 0
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "rgkmax", &
           & getstructscreening%rgkmax)
            Call removeAttribute (thisnode, "rgkmax")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "nempty")
         getstructscreening%nempty = 0
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "nempty", &
           & getstructscreening%nempty)
            Call removeAttribute (thisnode, "nempty")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "screentype")
         getstructscreening%screentype = "full"
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "screentype", &
           & getstructscreening%screentype)
            Call removeAttribute (thisnode, "screentype")
         End If
         getstructscreening%screentypenumber = stringtonumberscreentype &
        & (getstructscreening%screentype)
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructBSE (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (BSE_type), Pointer :: getstructBSE
         Type (Node), Pointer :: np
!
!
         Integer :: Len = 1, i = 0
         Allocate (getstructBSE)
#ifdef INPUTDEBUG
         Write (*,*) "we are at BSE"
#endif
!
         Nullify (np)
         np => getAttributeNode (thisnode, "nosym")
         getstructBSE%nosym = .False.
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "nosym", &
           & getstructBSE%nosym)
            Call removeAttribute (thisnode, "nosym")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "reducek")
         getstructBSE%reducek = .False.
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "reducek", &
           & getstructBSE%reducek)
            Call removeAttribute (thisnode, "reducek")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "vkloff")
         getstructBSE%vkloff = (/ - 1, - 1, - 1 /)
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "vkloff", &
           & getstructBSE%vkloff)
            Call removeAttribute (thisnode, "vkloff")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "rgkmax")
         getstructBSE%rgkmax = 0
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "rgkmax", &
           & getstructBSE%rgkmax)
            Call removeAttribute (thisnode, "rgkmax")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "scrherm")
         getstructBSE%scrherm = 0
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "scrherm", &
           & getstructBSE%scrherm)
            Call removeAttribute (thisnode, "scrherm")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "fbzq")
         getstructBSE%fbzq = .False.
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "fbzq", &
           & getstructBSE%fbzq)
            Call removeAttribute (thisnode, "fbzq")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "sciavtype")
         getstructBSE%sciavtype = "spherical"
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "sciavtype", &
           & getstructBSE%sciavtype)
            Call removeAttribute (thisnode, "sciavtype")
         End If
         getstructBSE%sciavtypenumber = stringtonumbersciavtype &
        & (getstructBSE%sciavtype)
!
         Nullify (np)
         np => getAttributeNode (thisnode, "sciavbd")
         getstructBSE%sciavbd = .False.
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "sciavbd", &
           & getstructBSE%sciavbd)
            Call removeAttribute (thisnode, "sciavbd")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "sciavqhd")
         getstructBSE%sciavqhd = .False.
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "sciavqhd", &
           & getstructBSE%sciavqhd)
            Call removeAttribute (thisnode, "sciavqhd")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "sciavqwg")
         getstructBSE%sciavqwg = .False.
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "sciavqwg", &
           & getstructBSE%sciavqwg)
            Call removeAttribute (thisnode, "sciavqwg")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "sciavqbd")
         getstructBSE%sciavqbd = .False.
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "sciavqbd", &
           & getstructBSE%sciavqbd)
            Call removeAttribute (thisnode, "sciavqbd")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "bsedirsing")
         getstructBSE%bsedirsing = .False.
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "bsedirsing", &
           & getstructBSE%bsedirsing)
            Call removeAttribute (thisnode, "bsedirsing")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "lmaxdielt")
         getstructBSE%lmaxdielt = 14
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "lmaxdielt", &
           & getstructBSE%lmaxdielt)
            Call removeAttribute (thisnode, "lmaxdielt")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "nleblaik")
         getstructBSE%nleblaik = 5810
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "nleblaik", &
           & getstructBSE%nleblaik)
            Call removeAttribute (thisnode, "nleblaik")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "nexcitmax")
         getstructBSE%nexcitmax = 100
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "nexcitmax", &
           & getstructBSE%nexcitmax)
            Call removeAttribute (thisnode, "nexcitmax")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "nstlbse")
         getstructBSE%nstlbse = (/ 0, 0 /)
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "nstlbse", &
           & getstructBSE%nstlbse)
            Call removeAttribute (thisnode, "nstlbse")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "nstlce")
         getstructBSE%nstlce = (/ 0, 0 /)
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "nstlce", &
           & getstructBSE%nstlce)
            Call removeAttribute (thisnode, "nstlce")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "bsetype")
         getstructBSE%bsetype = "singlet"
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "bsetype", &
           & getstructBSE%bsetype)
            Call removeAttribute (thisnode, "bsetype")
         End If
         getstructBSE%bsetypenumber = stringtonumberbsetype &
        & (getstructBSE%bsetype)
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructtetra (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (tetra_type), Pointer :: getstructtetra
         Type (Node), Pointer :: np
!
!
         Integer :: Len = 1, i = 0
         Allocate (getstructtetra)
#ifdef INPUTDEBUG
         Write (*,*) "we are at tetra"
#endif
!
         Nullify (np)
         np => getAttributeNode (thisnode, "tetraocc")
         getstructtetra%tetraocc = .False.
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "tetraocc", &
           & getstructtetra%tetraocc)
            Call removeAttribute (thisnode, "tetraocc")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "tetradf")
         getstructtetra%tetradf = .False.
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "tetradf", &
           & getstructtetra%tetradf)
            Call removeAttribute (thisnode, "tetradf")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "kordexc")
         getstructtetra%kordexc = .False.
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "kordexc", &
           & getstructtetra%kordexc)
            Call removeAttribute (thisnode, "kordexc")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "cw1k")
         getstructtetra%cw1k = .False.
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "cw1k", &
           & getstructtetra%cw1k)
            Call removeAttribute (thisnode, "cw1k")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "qweights")
         getstructtetra%qweights = 1
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "qweights", &
           & getstructtetra%qweights)
            Call removeAttribute (thisnode, "qweights")
         End If
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructdosWindow (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (dosWindow_type), Pointer :: getstructdosWindow
         Type (Node), Pointer :: np
!
!
         Integer :: Len = 1, i = 0
         Allocate (getstructdosWindow)
#ifdef INPUTDEBUG
         Write (*,*) "we are at dosWindow"
#endif
!
         Nullify (np)
         np => getAttributeNode (thisnode, "points")
         getstructdosWindow%points = 500
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "points", &
           & getstructdosWindow%points)
            Call removeAttribute (thisnode, "points")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "intv")
         getstructdosWindow%intv = (/ - 0.5, 0.5 /)
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "intv", &
           & getstructdosWindow%intv)
            Call removeAttribute (thisnode, "intv")
         End If
!
         Nullify (np)
         np => getAttributeNode (thisnode, "nsmdos")
         getstructdosWindow%nsmdos = 0
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "nsmdos", &
           & getstructdosWindow%nsmdos)
            Call removeAttribute (thisnode, "nsmdos")
         End If
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructplan (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (plan_type), Pointer :: getstructplan
		
         Integer :: Len = 1, i = 0
         Allocate (getstructplan)
#ifdef INPUTDEBUG
         Write (*,*) "we are at plan"
#endif
!
         Len = countChildEmentsWithName (thisnode, "doonly")
!
         Allocate (getstructplan%doonlyarray(len))
         Do i = 0, len - 1
            getstructplan%doonlyarray(i+1)%doonly => getstructdoonly &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "doonly"), 0)))
         End Do
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructdoonly (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (doonly_type), Pointer :: getstructdoonly
         Type (Node), Pointer :: np
!
!
         Integer :: Len = 1, i = 0
         Allocate (getstructdoonly)
#ifdef INPUTDEBUG
         Write (*,*) "we are at doonly"
#endif
!
         Nullify (np)
         np => getAttributeNode (thisnode, "task")
         If (associated(np)) Then
            Call extractDataAttribute (thisnode, "task", &
           & getstructdoonly%task)
            Call removeAttribute (thisnode, "task")
         End If
         getstructdoonly%tasknumber = stringtonumbertask &
        & (getstructdoonly%task)
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getstructqpointset (thisnode)
!
         Implicit None
         Type (Node), Pointer :: thisnode
         Type (qpointset_type), Pointer :: getstructqpointset
		
         Integer :: Len = 1, i = 0
         Allocate (getstructqpointset)
#ifdef INPUTDEBUG
         Write (*,*) "we are at qpointset"
#endif
!
         Len = countChildEmentsWithName (thisnode, "qpoint")
         Allocate (getstructqpointset%qpoint(3, len))
         Do i = 1, len
!
            getstructqpointset%qpoint (:, i) = getvalueofqpoint &
           & (removeChild(thisnode, item(getElementsByTagname(thisnode, &
           & "qpoint"), 0)))
         End Do
!
         i = 0
         Len = 0
         Call handleunknownnodes (thisnode)
      End Function
!
      Function getvalueofpointstatepair (thisnode)
         Implicit None
         Type (Node), Pointer :: thisnode
         Integer :: getvalueofpointstatepair (2)
!
#ifdef INPUTDEBUG
         Write (*,*) "we are at pointstatepair"
#endif
         Call extractDataContent (thisnode, getvalueofpointstatepair)
      End Function
      Function getvalueoftitle (thisnode)
         Implicit None
         Type (Node), Pointer :: thisnode
         Character (512) :: getvalueoftitle
!
#ifdef INPUTDEBUG
         Write (*,*) "we are at title"
#endif
         Call extractDataContent (thisnode, getvalueoftitle)
      End Function
      Function getvalueofbasevect (thisnode)
         Implicit None
         Type (Node), Pointer :: thisnode
         Real (8) :: getvalueofbasevect (3)
!
#ifdef INPUTDEBUG
         Write (*,*) "we are at basevect"
#endif
         Call extractDataContent (thisnode, getvalueofbasevect)
      End Function
      Function getvalueofoptcomp (thisnode)
         Implicit None
         Type (Node), Pointer :: thisnode
         Integer :: getvalueofoptcomp (3)
!
#ifdef INPUTDEBUG
         Write (*,*) "we are at optcomp"
#endif
         Call extractDataContent (thisnode, getvalueofoptcomp)
      End Function
      Function getvalueoftrans (thisnode)
         Implicit None
         Type (Node), Pointer :: thisnode
         Integer :: getvalueoftrans (3)
!
#ifdef INPUTDEBUG
         Write (*,*) "we are at trans"
#endif
         Call extractDataContent (thisnode, getvalueoftrans)
      End Function
      Function getvalueofqpoint (thisnode)
         Implicit None
         Type (Node), Pointer :: thisnode
         Real (8) :: getvalueofqpoint (3)
!
#ifdef INPUTDEBUG
         Write (*,*) "we are at qpoint"
#endif
         Call extractDataContent (thisnode, getvalueofqpoint)
      End Function
      Integer Function stringtonumberfixspin (string)
         Character (80), Intent (In) :: string
         Select Case (trim(adjustl(string)))
         Case ('none')
            stringtonumberfixspin = 0
         Case ('total FSM')
            stringtonumberfixspin = 1
         Case ('localmt FSM')
            stringtonumberfixspin = 2
         Case ('both')
            stringtonumberfixspin = 3
         Case ('')
            stringtonumberfixspin = 0
         Case Default
            Write (*,*) "'", string, "' is not valid selection forfixsp&
           &in "
            Stop
         End Select
      End Function
!
!
      Integer Function stringtonumbertype (string)
         Character (80), Intent (In) :: string
         Select Case (trim(adjustl(string)))
         Case ('Lapack')
            stringtonumbertype = 1
         Case ('Arpack')
            stringtonumbertype = 2
         Case ('DIIS')
            stringtonumbertype = 3
         Case ('')
            stringtonumbertype = 0
         Case Default
            Write (*,*) "'", string, "' is not valid selection fortype &
           &"
            Stop
         End Select
      End Function
!
!
      Integer Function stringtonumberdo (string)
         Character (80), Intent (In) :: string
         Select Case (trim(adjustl(string)))
         Case ('fromscratch')
            stringtonumberdo = - 1
         Case ('fromfile')
            stringtonumberdo = - 1
         Case ('skip')
            stringtonumberdo = - 1
         Case ('')
            stringtonumberdo = 0
         Case Default
            Write (*,*) "'", string, "' is not valid selection fordo "
            Stop
         End Select
      End Function
!
!
      Integer Function stringtonumberstype (string)
         Character (80), Intent (In) :: string
         Select Case (trim(adjustl(string)))
         Case ('Gaussian')
            stringtonumberstype = 0
         Case ('Methfessel-Paxton 1')
            stringtonumberstype = 1
         Case ('Methfessel-Paxton 2')
            stringtonumberstype = 2
         Case ('Fermi Dirac')
            stringtonumberstype = 3
         Case ('Square-wave impulse')
            stringtonumberstype = 4
         Case ('')
            stringtonumberstype = 0
         Case Default
            Write (*,*) "'", string, "' is not valid selection forstype&
           & "
            Stop
         End Select
      End Function
!
!
      Integer Function stringtonumberfindlinentype (string)
         Character (80), Intent (In) :: string
         Select Case (trim(adjustl(string)))
         Case ('simple')
            stringtonumberfindlinentype = - 1
         Case ('advanced')
            stringtonumberfindlinentype = - 1
         Case ('')
            stringtonumberfindlinentype = 0
         Case Default
            Write (*,*) "'", string, "' is not valid selection forfindl&
           &inentype "
            Stop
         End Select
      End Function
!
!
      Integer Function stringtonumbermixer (string)
         Character (80), Intent (In) :: string
         Select Case (trim(adjustl(string)))
         Case ('lin')
            stringtonumbermixer = 1
         Case ('msec')
            stringtonumbermixer = 2
         Case ('pulay')
            stringtonumbermixer = 3
         Case ('')
            stringtonumbermixer = 0
         Case Default
            Write (*,*) "'", string, "' is not valid selection formixer&
           & "
            Stop
         End Select
      End Function
!
!
      Integer Function stringtonumberxctype (string)
         Character (80), Intent (In) :: string
         Select Case (trim(adjustl(string)))
         Case ('LDAPerdew-Zunger')
            stringtonumberxctype = 2
         Case ('LSDAPerdew-Wang')
            stringtonumberxctype = 3
         Case ('LDA-X-alpha')
            stringtonumberxctype = 4
         Case ('LSDA-Barth-Hedin')
            stringtonumberxctype = 5
         Case ('GGAPerdew-Burke-Ernzerhof')
            stringtonumberxctype = 20
         Case ('GGArevPBE')
            stringtonumberxctype = 21
         Case ('GGAPBEsol')
            stringtonumberxctype = 22
         Case ('GGA-Wu-Cohen')
            stringtonumberxctype = 26
         Case ('GGAArmiento-Mattsson')
            stringtonumberxctype = 30
         Case ('EXX')
            stringtonumberxctype = - 2
         Case ('none')
            stringtonumberxctype = 0
         Case ('')
            stringtonumberxctype = 0
         Case Default
            Write (*,*) "'", string, "' is not valid selection forxctyp&
           &e "
            Stop
         End Select
      End Function
!
!
      Integer Function stringtonumberfxctype (string)
         Character (80), Intent (In) :: string
         Select Case (trim(adjustl(string)))
         Case ('RPA')
            stringtonumberfxctype = 0
         Case ('LRCstatic_NLF')
            stringtonumberfxctype = 1
         Case ('LRCstatic')
            stringtonumberfxctype = 2
         Case ('LRCdyn_NLF')
            stringtonumberfxctype = 3
         Case ('LRCdyn')
            stringtonumberfxctype = 4
         Case ('ALDA')
            stringtonumberfxctype = 5
         Case ('MB1_NLF')
            stringtonumberfxctype = 7
         Case ('MB1')
            stringtonumberfxctype = 8
         Case ('')
            stringtonumberfxctype = 0
         Case Default
            Write (*,*) "'", string, "' is not valid selection forfxcty&
           &pe "
            Stop
         End Select
      End Function
!
!
      Integer Function stringtonumberrun (string)
         Character (80), Intent (In) :: string
         Select Case (trim(adjustl(string)))
         Case ('fromscratch')
            stringtonumberrun = - 1
         Case ('skip')
            stringtonumberrun = - 1
         Case ('')
            stringtonumberrun = 0
         Case Default
            Write (*,*) "'", string, "' is not valid selection forrun "
            Stop
         End Select
      End Function
!
!
      Integer Function stringtonumberscreentype (string)
         Character (80), Intent (In) :: string
         Select Case (trim(adjustl(string)))
         Case ('full')
            stringtonumberscreentype = - 1
         Case ('diag')
            stringtonumberscreentype = - 1
         Case ('noinvdiag')
            stringtonumberscreentype = - 1
         Case ('longrange')
            stringtonumberscreentype = - 1
         Case ('')
            stringtonumberscreentype = 0
         Case Default
            Write (*,*) "'", string, "' is not valid selection forscree&
           &ntype "
            Stop
         End Select
      End Function
!
!
      Integer Function stringtonumbersciavtype (string)
         Character (80), Intent (In) :: string
         Select Case (trim(adjustl(string)))
         Case ('spherical')
            stringtonumbersciavtype = - 1
         Case ('screendiag')
            stringtonumbersciavtype = - 1
         Case ('invscreendiag')
            stringtonumbersciavtype = - 1
         Case ('')
            stringtonumbersciavtype = 0
         Case Default
            Write (*,*) "'", string, "' is not valid selection forsciav&
           &type "
            Stop
         End Select
      End Function
!
!
      Integer Function stringtonumberbsetype (string)
         Character (80), Intent (In) :: string
         Select Case (trim(adjustl(string)))
         Case ('ip')
            stringtonumberbsetype = - 1
         Case ('rpa')
            stringtonumberbsetype = - 1
         Case ('singlet')
            stringtonumberbsetype = - 1
         Case ('triplet')
            stringtonumberbsetype = - 1
         Case ('')
            stringtonumberbsetype = 0
         Case Default
            Write (*,*) "'", string, "' is not valid selection forbsety&
           &pe "
            Stop
         End Select
      End Function
!
!
      Integer Function stringtonumbertask (string)
         Character (80), Intent (In) :: string
         Select Case (trim(adjustl(string)))
         Case ('xsgeneigvec')
            stringtonumbertask = 301
         Case ('tetcalccw')
            stringtonumbertask = 310
         Case ('writepmatxs')
            stringtonumbertask = 320
         Case ('writeemat')
            stringtonumbertask = 330
         Case ('df')
            stringtonumbertask = 340
         Case ('df2')
            stringtonumbertask = 345
         Case ('idf')
            stringtonumbertask = 350
         Case ('scrgeneigvec')
            stringtonumbertask = 401
         Case ('scrtetcalccw')
            stringtonumbertask = 410
         Case ('scrwritepmat')
            stringtonumbertask = 420
         Case ('screen')
            stringtonumbertask = 430
         Case ('scrcoulint')
            stringtonumbertask = 440
         Case ('exccoulint')
            stringtonumbertask = 441
         Case ('BSE')
            stringtonumbertask = 445
         Case ('kernxc_bse')
            stringtonumbertask = 450
         Case ('writebandgapgrid')
            stringtonumbertask = 23
         Case ('writepmat')
            stringtonumbertask = 120
         Case ('dielectric')
            stringtonumbertask = 121
         Case ('writepmatasc')
            stringtonumbertask = 321
         Case ('pmatxs2orig')
            stringtonumbertask = 322
         Case ('writeematasc')
            stringtonumbertask = 331
         Case ('writepwmat')
            stringtonumbertask = 335
         Case ('emattest')
            stringtonumbertask = 339
         Case ('x0toasc')
            stringtonumbertask = 341
         Case ('x0tobin')
            stringtonumbertask = 342
         Case ('epsconv')
            stringtonumbertask = 396
         Case ('fxc_alda_check')
            stringtonumbertask = 398
         Case ('kernxc_bse3')
            stringtonumbertask = 451
         Case ('testxs')
            stringtonumbertask = 499
         Case ('xsestimate')
            stringtonumbertask = 700
         Case ('xstiming')
            stringtonumbertask = 701
         Case ('testmain')
            stringtonumbertask = 999
         Case ('portstate(1)')
            stringtonumbertask = 900
         Case ('portstate(2)')
            stringtonumbertask = 901
         Case ('portstate(-1)')
            stringtonumbertask = 910
         Case ('portstate(-2)')
            stringtonumbertask = 911
         Case ('')
            stringtonumbertask = 0
         Case Default
            Write (*,*) "'", string, "' is not valid selection fortask &
           &"
            Stop
         End Select
      End Function
!
!
      Integer Function stringtonumberxstype (string)
         Character (80), Intent (In) :: string
         Select Case (trim(adjustl(string)))
         Case ('TDDFT')
            stringtonumberxstype = - 1
         Case ('BSE')
            stringtonumberxstype = - 1
         Case ('')
            stringtonumberxstype = 0
         Case Default
            Write (*,*) "'", string, "' is not valid selection forxstyp&
           &e "
            Stop
         End Select
      End Function
!
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
! these are some transient helper functions to simplify the port (should not be used)
      Function isspinorb ()
         Logical :: isspinorb
         isspinorb = .False.
         If (associated(input%groundstate%spin)) Then
            If (input%groundstate%spin%spinorb) Then
               isspinorb = .True.
            End If
         End If
      End Function
      Function isspinspiral ()
         Logical :: isspinspiral
         isspinspiral = .False.
         If (associated(input%groundstate%spin)) Then
            If (input%groundstate%spin%spinsprl) Then
               isspinspiral = .True.
            End If
         End If
      End Function
!
      Function getfixspinnumber ()
         Implicit None
         Integer :: getfixspinnumber
         getfixspinnumber = 0
         If (associated(input%groundstate%spin)) Then
            getfixspinnumber = input%groundstate%spin%fixspinnumber
         End If
      End Function
!
      Function istetraocc ()
         Implicit None
         Logical :: istetraocc
         istetraocc = .False.
         If (associated(input%xs)) Then
            If (associated(input%xs%tetra)) Then
               istetraocc = input%xs%tetra%tetraocc
            End If
         End If
      End Function
!
End Module
!
