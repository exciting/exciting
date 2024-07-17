""" Automatically generated file with the valid attributes from the schema. 
Do not manually change. Instead, run "utils/schema_parsing.py" to regenerate. """ 

# input information 
input_attribute_types = {"sharedfs": (bool, 1), "xsltpath": (str, 1)} 
input_valid_subtrees = ["title", "structure", "groundstate", "relax", "properties", "phonons", "xs", "gw", "MD", "eph",
                        "keywords"] 
input_mandatory_attributes = ["groundstate", "structure", "title"] 


# common information 
origin_attribute_types = {"coord": (float, 3)} 

point_attribute_types = {"breakafter": (bool, 1), "coord": (float, 3), "label": (str, 1)} 
point_mandatory_attributes = ["coord"] 

plot1d_valid_subtrees = ["path"] 
plot1d_mandatory_attributes = ["path"] 

path_attribute_types = {"outfileprefix": (str, 1), "steps": (int, 1)} 
path_valid_subtrees = ["point"] 
path_mandatory_attributes = ["point", "steps"] 
path_multiple_children = ["point"] 

plot2d_valid_subtrees = ["parallelogram"] 
plot2d_mandatory_attributes = ["parallelogram"] 

parallelogram_attribute_types = {"grid": (int, 2), "outfileprefix": (str, 1)} 
parallelogram_valid_subtrees = ["origin", "point"] 
parallelogram_mandatory_attributes = ["grid", "origin", "point"] 
parallelogram_multiple_children = ["point"] 

plot3d_attribute_types = {"usesym": (bool, 1)} 
plot3d_valid_subtrees = ["box"] 
plot3d_mandatory_attributes = ["box"] 

box_attribute_types = {"grid": (int, 3), "outfileprefix": (str, 1)} 
box_valid_subtrees = ["origin", "point"] 
box_mandatory_attributes = ["grid", "origin", "point"] 
box_multiple_children = ["point"] 

kstlist_valid_subtrees = ["pointstatepair"] 
kstlist_mandatory_attributes = ["pointstatepair"] 
kstlist_multiple_children = ["pointstatepair"] 

energywindow_attribute_types = {"intv": (float, 2), "points": (int, 1)} 

qpointset_valid_subtrees = ["qpoint"] 
qpointset_mandatory_attributes = ["qpoint"] 
qpointset_multiple_children = ["qpoint"] 

parts_valid_subtrees = ["dopart"] 
parts_multiple_children = ["dopart"] 

dopart_attribute_types = {"id": (str, 1)} 
dopart_mandatory_attributes = ["id"] 

qpoints_attribute_types = {"qf": (int, 1), "qi": (int, 1)} 

kpoints_attribute_types = {"kf": (int, 1), "ki": (int, 1)} 


# structure information 
structure_attribute_types = {"autormt": (bool, 1),
                             "autormtscaling": (float, 1),
                             "cartesian": (bool, 1),
                             "epslat": (float, 1),
                             "primcell": (bool, 1),
                             "speciespath": (str, 1),
                             "tshift": (bool, 1)} 
structure_valid_subtrees = ["crystal", "species", "symmetries"] 
structure_mandatory_attributes = ["speciespath"] 
structure_multiple_children = ["species"] 

crystal_attribute_types = {"scale": (float, 1), "stretch": (float, 3)} 
crystal_valid_subtrees = ["basevect"] 
crystal_multiple_children = ["basevect"] 

species_attribute_types = {"atomicNumber": (int, 1),
                           "chemicalSymbol": (str, 1),
                           "fixrmt": (bool, 1),
                           "rmt": (float, 1),
                           "speciesfile": (str, 1)} 
species_valid_subtrees = ["atom", "LDAplusU", "dfthalfparam"] 
species_mandatory_attributes = ["speciesfile"] 
species_multiple_children = ["atom"] 

atom_attribute_types = {"bfcmt": (float, 3),
                        "coord": (float, 3),
                        "lockxyz": (bool, 3),
                        "mommtfix": (float, 3),
                        "velocity": (float, 3)} 
atom_mandatory_attributes = ["coord"] 

LDAplusU_attribute_types = {"J": (float, 1), "U": (float, 1), "l": (int, 1)} 

dfthalfparam_attribute_types = {"ampl": (float, 1), "cut": (float, 1), "exponent": (int, 1)} 
dfthalfparam_valid_subtrees = ["shell"] 
dfthalfparam_mandatory_attributes = ["shell"] 
dfthalfparam_multiple_children = ["shell"] 

shell_attribute_types = {"ionization": (float, 1), "number": (int, 1)} 


# groundstate information 
groundstate_attribute_types = {"APWprecision": (float, 1),
                               "CoreRelativity": (str, ["dirac", "none"]),
                               "ExplicitKineticEnergy": (bool, 1),
                               "PrelimLinSteps": (int, 1),
                               "ValenceRelativity": (str, ["iora", "iora*", "kh", "kh*", "none", "zora"]),
                               "autokpt": (bool, 1),
                               "beta0": (float, 1),
                               "betadec": (float, 1),
                               "betainc": (float, 1),
                               "cfdamp": (float, 1),
                               "chgexs": (float, 1),
                               "deband": (float, 1),
                               "dipolecorrection": (bool, 1),
                               "dipoleposition": (float, 1),
                               "dlinengyfermi": (float, 1),
                               "do": (str, ["fromfile", "fromscratch", "skip"]),
                               "energyref": (float, 1),
                               "epsband": (float, 1),
                               "epschg": (float, 1),
                               "epsengy": (float, 1),
                               "epsforcescf": (float, 1),
                               "epsocc": (float, 1),
                               "epspot": (float, 1),
                               "fermilinengy": (bool, 1),
                               "findlinentype": (str, ["Wigner_Seitz", "lcharge", "logderiv", "no_search"]),
                               "fracinr": (float, 1),
                               "frozencore": (bool, 1),
                               "gmaxvr": (float, 1),
                               "isgkmax": (int, 1),
                               "ldapu": (str,
                                         ["AroundMeanField", "FFL-AMF-interpolation", "FullyLocalisedLimit", "none"]),
                               "lmaxapw": (int, 1),
                               "lmaxinr": (int, 1),
                               "lmaxmat": (int, 1),
                               "lmaxvr": (int, 1),
                               "lradstep": (int, 1),
                               "maxscl": (int, 1),
                               "mixer": (str, ["lin", "msec", "pulay"]),
                               "mixerswitch": (int, 1),
                               "modifiedsv": (bool, 1),
                               "msecStoredSteps": (int, 1),
                               "nempty": (int, 1),
                               "ngridk": (int, 3),
                               "niterconvcheck": (int, 1),
                               "nktot": (int, 1),
                               "nosource": (bool, 1),
                               "nosym": (bool, 1),
                               "nprad": (int, 1),
                               "npsden": (int, 1),
                               "nwrite": (int, 1),
                               "outputlevel": (str, ["high", "low", "none", "normal"]),
                               "ptnucl": (bool, 1),
                               "radialgridtype": (str, ["cubic", "cubic-2", "expocubic", "exponential"]),
                               "radkpt": (float, 1),
                               "reducek": (bool, 1),
                               "rgkmax": (float, 1),
                               "scfconv": (str, 1),
                               "stype": (str,
                                         ["Fermi Dirac", "Gaussian", "Methfessel-Paxton 1", "Methfessel-Paxton 2",
                                          "Square-wave impulse", "libbzint"]),
                               "swidth": (float, 1),
                               "symmorph": (bool, 1),
                               "tevecsv": (bool, 1),
                               "tfibs": (bool, 1),
                               "tforce": (bool, 1),
                               "tpartcharges": (bool, 1),
                               "useAPWprecision": (bool, 1),
                               "useDensityMatrix": (bool, 1),
                               "vdWcorrection": (str, ["DFTD2", "TSvdW", "none"]),
                               "vkloff": (float, 3),
                               "xctype": (str,
                                          ["EXX", "GGA_AC_PBE", "GGA_AM05", "GGA_PBE", "GGA_PBE_R", "GGA_PBE_SOL",
                                           "GGA_PBE_SR", "GGA_WC", "HYB_HSE", "HYB_LDA0", "HYB_PBE0", "LDA_PW",
                                           "LDA_PZ", "LDA_XALPHA", "LDA_vBH", "none"])} 
groundstate_valid_subtrees = ["DFTD2parameters", "TSvdWparameters", "spin", "HartreeFock", "dfthalf", "Hybrid",
                              "sirius", "solver", "OEP", "RDMFT", "output", "libxc", "xsLO", "lorecommendation"] 

DFTD2parameters_attribute_types = {"cutoff": (float, 1), "d": (float, 1), "s6": (float, 1), "sr6": (float, 1)} 

TSvdWparameters_attribute_types = {"cutoff": (float, 1),
                                   "d": (float, 1),
                                   "nr": (int, 1),
                                   "nsph": (int, 1),
                                   "s6": (float, 1),
                                   "sr6": (float, 1)} 

spin_attribute_types = {"bfieldc": (float, 3),
                        "fixspin": (str, ["both", "localmt FSM", "none", "total FSM"]),
                        "momfix": (float, 3),
                        "nosv": (bool, 1),
                        "realspace": (bool, 1),
                        "reducebf": (float, 1),
                        "spinorb": (bool, 1),
                        "spinsprl": (bool, 1),
                        "svlo": (bool, 1),
                        "taufsm": (float, 1),
                        "vqlss": (float, 3)} 

dfthalf_attribute_types = {"printVSfile": (bool, 1)} 

Hybrid_attribute_types = {"BasisBareCoulomb": (str, 1),
                          "HSEsingularity": (str, ["Exact", "Taylor"]),
                          "eccoeff": (float, 1),
                          "epsmb": (float, 1),
                          "exchangetype": (str, ["HF", "OEP"]),
                          "excoeff": (float, 1),
                          "gmb": (float, 1),
                          "lmaxmb": (int, 1),
                          "maxscl": (int, 1),
                          "mblksiz": (int, 1),
                          "omega": (float, 1),
                          "updateRadial": (bool, 1)} 

sirius_attribute_types = {"cfun": (bool, 1),
                          "density": (bool, 1),
                          "densityinit": (bool, 1),
                          "eigenstates": (bool, 1),
                          "sfacg": (bool, 1),
                          "vha": (bool, 1),
                          "xc": (bool, 1)} 

solver_attribute_types = {"constructHS": (bool, 1),
                          "evaltol": (float, 1),
                          "minenergy": (float, 1),
                          "packedmatrixstorage": (bool, 1),
                          "type": (str, ["Davidson", "Lapack"])} 

OEP_attribute_types = {"convoep": (float, 1), "maxitoep": (int, 1), "tauoep": (float, 3)} 

RDMFT_attribute_types = {"maxitc": (int, 1),
                         "maxitn": (int, 1),
                         "rdmalpha": (float, 1),
                         "rdmmaxscl": (int, 1),
                         "rdmtemp": (float, 1),
                         "rdmxctype": (int, 1),
                         "taurdmc": (float, 1),
                         "taurdmn": (float, 1)} 

output_attribute_types = {"state": (str, ["XML", "binary"])} 

libxc_attribute_types = {"correlation": (str,
                                         ["XC_GGA_C_AM05", "XC_GGA_C_APBE", "XC_GGA_C_FT97", "XC_GGA_C_LM",
                                          "XC_GGA_C_LYP", "XC_GGA_C_OPTC", "XC_GGA_C_OP_B88", "XC_GGA_C_OP_G96",
                                          "XC_GGA_C_OP_PBE", "XC_GGA_C_OP_XALPHA", "XC_GGA_C_P86", "XC_GGA_C_PBE",
                                          "XC_GGA_C_PBE_JRGX", "XC_GGA_C_PBE_SOL", "XC_GGA_C_PW91", "XC_GGA_C_REVTCA",
                                          "XC_GGA_C_RGE2", "XC_GGA_C_SOGGA11", "XC_GGA_C_SOGGA11_X", "XC_GGA_C_SPBE",
                                          "XC_GGA_C_TCA", "XC_GGA_C_WI", "XC_GGA_C_WI0", "XC_GGA_C_WL", "XC_GGA_C_XPBE",
                                          "XC_LDA_C_1D_CSC", "XC_LDA_C_1D_LOOS", "XC_LDA_C_2D_AMGB", "XC_LDA_C_2D_PRM",
                                          "XC_LDA_C_GL", "XC_LDA_C_GOMBAS", "XC_LDA_C_HL", "XC_LDA_C_ML1",
                                          "XC_LDA_C_ML2", "XC_LDA_C_OB_PW", "XC_LDA_C_OB_PZ", "XC_LDA_C_PW",
                                          "XC_LDA_C_PW_MOD", "XC_LDA_C_PW_RPA", "XC_LDA_C_PZ", "XC_LDA_C_PZ_MOD",
                                          "XC_LDA_C_RC04", "XC_LDA_C_RPA", "XC_LDA_C_VWN", "XC_LDA_C_VWN_1",
                                          "XC_LDA_C_VWN_2", "XC_LDA_C_VWN_3", "XC_LDA_C_VWN_4", "XC_LDA_C_VWN_RPA",
                                          "XC_LDA_C_WIGNER", "XC_LDA_C_XALPHA", "XC_LDA_C_vBH", "none"]),
                         "exchange": (str,
                                      ["XC_GGA_X_2D_B86", "XC_GGA_X_2D_B86_MGC", "XC_GGA_X_2D_B88", "XC_GGA_X_2D_PBE",
                                       "XC_GGA_X_AIRY", "XC_GGA_X_AM05", "XC_GGA_X_APBE", "XC_GGA_X_B86",
                                       "XC_GGA_X_B86_MGC", "XC_GGA_X_B88", "XC_GGA_X_BAYESIAN", "XC_GGA_X_BPCCAC",
                                       "XC_GGA_X_C09X", "XC_GGA_X_DK87_R1", "XC_GGA_X_DK87_R2", "XC_GGA_X_FT97_A",
                                       "XC_GGA_X_FT97_B", "XC_GGA_X_G96", "XC_GGA_X_HERMAN", "XC_GGA_X_HTBS",
                                       "XC_GGA_X_KT1", "XC_GGA_X_LAG", "XC_GGA_X_LB", "XC_GGA_X_LBM", "XC_GGA_X_LG93",
                                       "XC_GGA_X_MB88", "XC_GGA_X_MPBE", "XC_GGA_X_MPW91", "XC_GGA_X_OL2",
                                       "XC_GGA_X_OPTB88_VDW", "XC_GGA_X_OPTPBE_VDW", "XC_GGA_X_OPTX", "XC_GGA_X_PBE",
                                       "XC_GGA_X_PBEA", "XC_GGA_X_PBEK1_VDW", "XC_GGA_X_PBE_JSJR", "XC_GGA_X_PBE_R",
                                       "XC_GGA_X_PBE_SOL", "XC_GGA_X_PW86", "XC_GGA_X_PW91", "XC_GGA_X_RGE2",
                                       "XC_GGA_X_RPBE", "XC_GGA_X_RPW86", "XC_GGA_X_SOGGA", "XC_GGA_X_SOGGA11",
                                       "XC_GGA_X_SSB", "XC_GGA_X_SSB_D", "XC_GGA_X_SSB_SW", "XC_GGA_X_WC",
                                       "XC_GGA_X_XPBE", "XC_LDA_X", "XC_LDA_X_1D", "XC_LDA_X_2D", "none"]),
                         "xc": (str,
                                ["XC_GGA_XC_B97", "XC_GGA_XC_B97_1", "XC_GGA_XC_B97_2", "XC_GGA_XC_B97_3",
                                 "XC_GGA_XC_B97_D", "XC_GGA_XC_B97_GGA1", "XC_GGA_XC_B97_K", "XC_GGA_XC_EDF1",
                                 "XC_GGA_XC_HCTH_120", "XC_GGA_XC_HCTH_147", "XC_GGA_XC_HCTH_407",
                                 "XC_GGA_XC_HCTH_407P", "XC_GGA_XC_HCTH_93", "XC_GGA_XC_HCTH_A", "XC_GGA_XC_HCTH_P14",
                                 "XC_GGA_XC_HCTH_P76", "XC_GGA_XC_KT2", "XC_GGA_XC_MOHLYP", "XC_GGA_XC_MOHLYP2",
                                 "XC_GGA_XC_MPWLYP1W", "XC_GGA_XC_PBE1W", "XC_GGA_XC_PBELYP1W", "XC_GGA_XC_SB98_1a",
                                 "XC_GGA_XC_SB98_1b", "XC_GGA_XC_SB98_1c", "XC_GGA_XC_SB98_2a", "XC_GGA_XC_SB98_2b",
                                 "XC_GGA_XC_SB98_2c", "XC_GGA_XC_TH1", "XC_GGA_XC_TH2", "XC_GGA_XC_TH3",
                                 "XC_GGA_XC_TH4", "XC_GGA_XC_TH_FC", "XC_GGA_XC_TH_FCFO", "XC_GGA_XC_TH_FCO",
                                 "XC_GGA_XC_TH_FL", "XC_GGA_XC_XLYP", "XC_HYB_GGA_XC_B1LYP", "XC_HYB_GGA_XC_B1PW91",
                                 "XC_HYB_GGA_XC_B1WC", "XC_HYB_GGA_XC_B3LYP", "XC_HYB_GGA_XC_B3P86",
                                 "XC_HYB_GGA_XC_B3PW91", "XC_HYB_GGA_XC_B97", "XC_HYB_GGA_XC_B97_1",
                                 "XC_HYB_GGA_XC_B97_2", "XC_HYB_GGA_XC_B97_3", "XC_HYB_GGA_XC_B97_K",
                                 "XC_HYB_GGA_XC_BHANDH", "XC_HYB_GGA_XC_BHANDHLYP", "XC_HYB_GGA_XC_MB3LYP_RC04",
                                 "XC_HYB_GGA_XC_MPW3LYP", "XC_HYB_GGA_XC_MPW3PW", "XC_HYB_GGA_XC_O3LYP",
                                 "XC_HYB_GGA_XC_PBEH", "XC_HYB_GGA_XC_SB98_1a", "XC_HYB_GGA_XC_SB98_1b",
                                 "XC_HYB_GGA_XC_SB98_1c", "XC_HYB_GGA_XC_SB98_2a", "XC_HYB_GGA_XC_SB98_2b",
                                 "XC_HYB_GGA_XC_SB98_2c", "XC_HYB_GGA_XC_X3LYP", "XC_HYB_GGA_XC_mPW1K",
                                 "XC_HYB_GGA_XC_mPW1PW", "XC_LDA_XC_TETER93", "none"])} 

xsLO_attribute_types = {"emax": (float, 1), "lmax": (int, 1), "maxnodes": (int, 1)} 

lorecommendation_attribute_types = {"lmaxlo": (int, 1), "nodesmaxlo": (int, 1)} 


# relax information 
relax_attribute_types = {"addtohistory": (bool, 1),
                         "endbfgs": (str, 1),
                         "epsforce": (float, 1),
                         "history": (bool, 1),
                         "historyformat": (str, 1),
                         "maxbfgs": (int, 1),
                         "maxsteps": (int, 1),
                         "method": (str, 1),
                         "outputlevel": (str, ["high", "low", "normal"]),
                         "printtorque": (bool, 1),
                         "taubfgs": (float, 1),
                         "taunewton": (float, 1)} 


# phonons information 
phonons_attribute_types = {"canonical": (bool, 1),
                           "delete_eigensystem_response": (bool, 1),
                           "deltaph": (float, 1),
                           "do": (str, ["dry", "fromfile", "fromscratch", "skip"]),
                           "drynumprocs": (int, 1),
                           "gamma": (str, ["onestep", "standard", "twostep"]),
                           "maxprocsperpart": (int, 1),
                           "method": (str, ["dfpt", "sc"]),
                           "minprocsperpart": (int, 1),
                           "ngridq": (int, 3),
                           "polar": (bool, 1),
                           "reduceq": (bool, 1),
                           "sumrule": (bool, 1),
                           "write_schedule": (bool, 1)} 
phonons_valid_subtrees = ["qpointset", "phonondos", "phonondispplot", "reformatdynmat", "interpolate", "parts"] 

phonondos_attribute_types = {"inttype": (str, ["sum", "tetra"]),
                             "ngrdos": (int, 1),
                             "ngridqint": (int, 3),
                             "nsmdos": (int, 1),
                             "ntemp": (int, 1),
                             "nwdos": (int, 1)} 

phonondispplot_valid_subtrees = ["plot1d"] 
phonondispplot_mandatory_attributes = ["plot1d"] 

interpolate_attribute_types = {"ngridq": (int, 3), "vqloff": (float, 3), "writeeigenvectors": (bool, 1)} 
interpolate_mandatory_attributes = ["ngridq"] 


# properties information 
properties_valid_subtrees = ["spintext", "coreoverlap", "bandstructure", "stm", "wfplot", "dos", "LSJ", "masstensor",
                             "chargedensityplot", "TSvdW", "DFTD2", "exccplot", "elfplot", "mvecfield", "xcmvecfield",
                             "electricfield", "gradmvecfield", "fermisurfaceplot", "EFG", "mossbauer", "expiqr",
                             "elnes", "eliashberg", "momentummatrix", "dielmat", "boltzequ", "raman", "moke", "shg",
                             "wannier", "wannierplot", "wanniergap", "ldos", "polarization"] 

spintext_attribute_types = {"bands": (int, 2)} 
spintext_valid_subtrees = ["plot2d"] 
spintext_mandatory_attributes = ["plot2d"] 

coreoverlap_attribute_types = {"coreatom": (int, 1), "corespecies": (int, 1)} 

bandstructure_attribute_types = {"character": (bool, 1),
                                 "deriv": (bool, 1),
                                 "scissor": (float, 1),
                                 "wannier": (bool, 1)} 
bandstructure_valid_subtrees = ["plot1d"] 
bandstructure_mandatory_attributes = ["plot1d"] 

stm_attribute_types = {"bias": (float, 1),
                       "stmmode": (str, ["constantHeight", "topographic"]),
                       "stmtype": (str, ["differentialConductance", "integratedLDOS"])} 
stm_valid_subtrees = ["plot2d", "region"] 

region_attribute_types = {"grid2d": (int, 2), "grid3d": (int, 3), "height": (float, 1), "zrange": (float, 2)} 

wfplot_attribute_types = {"version": (str, 1)} 
wfplot_valid_subtrees = ["kstlist", "plot1d", "plot2d", "plot3d"] 
wfplot_mandatory_attributes = ["kstlist"] 

dos_attribute_types = {"inttype": (str, ["tetra", "trilin", "trilin+"]),
                       "jdos": (bool, 1),
                       "linkpt": (float, 1),
                       "lmirep": (bool, 1),
                       "lonly": (bool, 1),
                       "ngrdos": (int, 1),
                       "ngridkint": (int, 3),
                       "nsmdos": (int, 1),
                       "nwdos": (int, 1),
                       "scissor": (float, 1),
                       "sqados": (float, 3),
                       "wannier": (bool, 1),
                       "winddos": (float, 2)} 

LSJ_valid_subtrees = ["kstlist"] 

masstensor_attribute_types = {"deltaem": (float, 1), "ndspem": (int, 1), "vklem": (float, 3)} 

chargedensityplot_attribute_types = {"nocore": (bool, 1)} 
chargedensityplot_valid_subtrees = ["plot1d", "plot2d", "plot3d"] 

exccplot_valid_subtrees = ["plot1d", "plot2d", "plot3d"] 

elfplot_valid_subtrees = ["plot1d", "plot2d", "plot3d"] 

mvecfield_valid_subtrees = ["plot2d", "plot3d"] 

xcmvecfield_valid_subtrees = ["plot2d", "plot3d"] 

electricfield_valid_subtrees = ["plot2d", "plot3d"] 

gradmvecfield_valid_subtrees = ["plot1d", "plot2d", "plot3d"] 

fermisurfaceplot_attribute_types = {"nstfsp": (int, 1)} 
fermisurfaceplot_valid_subtrees = ["plot2d", "plot3d"] 

expiqr_valid_subtrees = ["kstlist"] 

elnes_attribute_types = {"ngrid": (int, 1),
                         "vecql": (float, 3),
                         "wgrid": (int, 1),
                         "wmax": (float, 1),
                         "wmin": (float, 1)} 

eliashberg_attribute_types = {"mustar": (float, 1)} 

momentummatrix_attribute_types = {"fastpmat": (bool, 1)} 

dielmat_attribute_types = {"drude": (float, 2),
                           "intraband": (bool, 1),
                           "scissor": (float, 1),
                           "swidth": (float, 1),
                           "tevout": (bool, 1),
                           "wgrid": (int, 1),
                           "wmax": (float, 1)} 
dielmat_valid_subtrees = ["epscomp"] 
dielmat_multiple_children = ["epscomp"] 

boltzequ_attribute_types = {"chemicalPotentialRange": (float, 2),
                            "chemicalPotentialSpacing": (float, 1),
                            "dopingConcentration": (float, 1),
                            "energyReference": (str, ["efermi", "none"]),
                            "evOutputEnergies": (bool, 1),
                            "siOutputUnits": (bool, 1),
                            "temperatureRange": (float, 2),
                            "temperatureSpacing": (float, 1),
                            "transportDfBroadening": (float, 1),
                            "transportDfRange": (float, 2),
                            "transportDfSpacing": (float, 1),
                            "useDopingConcentration": (bool, 1),
                            "useTransportDf": (bool, 1)} 
boltzequ_valid_subtrees = ["etCoeffComponents"] 
boltzequ_multiple_children = ["etCoeffComponents"] 

raman_attribute_types = {"broad": (float, 1),
                         "degree": (int, 1),
                         "displ": (float, 1),
                         "doequilibrium": (bool, 1),
                         "elaser": (float, 1),
                         "elaserunit": (str, ["Ha", "cm-1", "eV", "nm"]),
                         "getphonon": (str, ["fromfile", "fromscratch", "readinput", "symvec", "symveccheck"]),
                         "mode": (int, 1),
                         "molecule": (bool, 1),
                         "ninter": (int, 1),
                         "nstate": (int, 1),
                         "nstep": (int, 1),
                         "temp": (float, 1),
                         "useforces": (bool, 1),
                         "usesym": (bool, 1),
                         "writefunc": (bool, 1),
                         "xmax": (float, 1),
                         "xmin": (float, 1)} 
raman_valid_subtrees = ["eigvec", "energywindow"] 
raman_mandatory_attributes = ["energywindow"] 
raman_multiple_children = ["eigvec"] 

eigvec_attribute_types = {"comp": (float, 2)} 
eigvec_mandatory_attributes = ["comp"] 

moke_attribute_types = {"drude": (float, 2),
                        "intraband": (bool, 1),
                        "scissor": (float, 1),
                        "swidth": (float, 1),
                        "tevout": (bool, 1),
                        "wgrid": (int, 1),
                        "wmax": (float, 1)} 

shg_attribute_types = {"etol": (float, 1),
                       "scissor": (float, 1),
                       "swidth": (float, 1),
                       "tevout": (bool, 1),
                       "wgrid": (int, 1),
                       "wmax": (float, 1)} 
shg_valid_subtrees = ["chicomp"] 
shg_mandatory_attributes = ["chicomp"] 
shg_multiple_children = ["chicomp"] 

wannier_attribute_types = {"cutshell": (bool, 1),
                           "do": (str, ["fromfile", "fromscratch", "maxfromfile", "skip"]),
                           "fermizero": (bool, 1),
                           "input": (str, ["gs", "gw", "hybrid", "qsgw"]),
                           "mindist": (bool, 1),
                           "minshell": (int, 1),
                           "nbzshell": (int, 3),
                           "printproj": (bool, 1)} 
wannier_valid_subtrees = ["projection", "group"] 
wannier_multiple_children = ["group"] 

projection_attribute_types = {"dordmax": (int, 1), "epsld": (float, 1), "nprojtot": (int, 1), "nunocc": (int, 1)} 

group_attribute_types = {"epsdis": (float, 1),
                         "epsmax": (float, 1),
                         "epsopf": (float, 1),
                         "epsproj": (float, 1),
                         "fst": (int, 1),
                         "innerwindow": (float, 2),
                         "lst": (int, 1),
                         "maxitdis": (int, 1),
                         "maxitmax": (int, 1),
                         "maxitopf": (int, 1),
                         "memlendis": (int, 1),
                         "memlenmax": (int, 1),
                         "memlenopf": (int, 1),
                         "method": (str, ["auto", "disFull", "disSMV", "opf", "opfmax", "pro", "promax"]),
                         "minitdis": (int, 1),
                         "minitmax": (int, 1),
                         "minitopf": (int, 1),
                         "minstepdis": (float, 1),
                         "minstepmax": (float, 1),
                         "minstepopf": (float, 1),
                         "neighcells": (bool, 1),
                         "nproj": (int, 1),
                         "nwf": (int, 1),
                         "nwrite": (int, 1),
                         "optim": (str, ["cg", "lbfgs"]),
                         "outerwindow": (float, 2),
                         "writeconv": (bool, 1)} 
group_valid_subtrees = ["projector"] 
group_multiple_children = ["projector"] 

projector_attribute_types = {"nr": (int, 1)} 
projector_mandatory_attributes = ["nr"] 

wannierplot_attribute_types = {"cell": (int, 3), "fst": (int, 1), "lst": (int, 1)} 
wannierplot_valid_subtrees = ["plot1d", "plot2d", "plot3d"] 

wanniergap_attribute_types = {"auto": (bool, 1), "ngridkint": (int, 3)} 
wanniergap_valid_subtrees = ["pointband"] 
wanniergap_multiple_children = ["pointband"] 

pointband_attribute_types = {"band": (int, 1), "extremal": (bool, 1), "vkl": (float, 3)} 
pointband_mandatory_attributes = ["band", "vkl"] 

ldos_attribute_types = {"delta": (float, 1),
                        "grid": (int, 3),
                        "newint": (bool, 1),
                        "ngrdos": (int, 1),
                        "nsmdos": (int, 1),
                        "nwdos": (int, 1),
                        "scissor": (float, 1),
                        "tol": (float, 1),
                        "winddos": (float, 2)} 


# xs information 
xs_attribute_types = {"bfieldc": (float, 3),
                      "broad": (float, 1),
                      "dbglev": (int, 1),
                      "dfoffdiag": (bool, 1),
                      "dogroundstate": (str, ["fromfile", "fromscratch"]),
                      "emattype": (int, 1),
                      "emaxdf": (float, 1),
                      "epsdfde": (float, 1),
                      "fastpmat": (bool, 1),
                      "gqmax": (float, 1),
                      "gqmaxtype": (str, ["|G+q|", "|G|"]),
                      "h5fname": (str, 1),
                      "h5gname": (str, 1),
                      "lmaxapwwf": (int, 1),
                      "lmaxemat": (int, 1),
                      "maxscl": (int, 1),
                      "nempty": (int, 1),
                      "ngridk": (int, 3),
                      "ngridq": (int, 3),
                      "nosym": (bool, 1),
                      "pwmat": (str, ["fft", "mm"]),
                      "reducek": (bool, 1),
                      "reduceq": (bool, 1),
                      "rgkmax": (float, 1),
                      "scissor": (float, 1),
                      "skipgnd": (bool, 1),
                      "swidth": (float, 1),
                      "tappinfo": (bool, 1),
                      "tevout": (bool, 1),
                      "vkloff": (float, 3),
                      "writexsgrids": (bool, 1),
                      "xstype": (str, ["BSE", "RT-TDDFT", "TDDFT", "fastBSE"])} 
xs_valid_subtrees = ["storeexcitons", "pwelements", "writeexcitons", "writekpathweights", "excitonPlot",
                     "realTimeTDDFT", "tddft", "screening", "phonon_screening", "expand_eps", "BSE", "fastBSE",
                     "transitions", "qpointset", "tetra", "energywindow", "plan"] 
xs_mandatory_attributes = ["xstype"] 

storeexcitons_attribute_types = {"MaxEnergyExcitons": (float, 1),
                                 "MaxNumberExcitons": (int, 1),
                                 "MinEnergyExcitons": (float, 1),
                                 "MinNumberExcitons": (int, 1),
                                 "selectenergy": (bool, 1),
                                 "useev": (bool, 1)} 

pwelements_attribute_types = {"band_combinations": (str, ["all", "occ_occ", "occ_unocc", "unocc_occ", "unocc_unocc"])} 
pwelements_mandatory_attributes = ["band_combinations"] 

writeexcitons_attribute_types = {"MaxEnergyExcitons": (float, 1),
                                 "MaxNumberExcitons": (int, 1),
                                 "MinEnergyExcitons": (float, 1),
                                 "MinNumberExcitons": (int, 1),
                                 "abscutares": (float, 2),
                                 "abscutres": (float, 2),
                                 "selectenergy": (bool, 1),
                                 "useev": (bool, 1)} 

writekpathweights_attribute_types = {"MaxEnergyExcitons": (float, 1),
                                     "MaxNumberExcitons": (int, 1),
                                     "MinEnergyExcitons": (float, 1),
                                     "MinNumberExcitons": (int, 1),
                                     "intorder": (int, 1),
                                     "printgridweights": (bool, 1),
                                     "selectenergy": (bool, 1),
                                     "useev": (bool, 1)} 

excitonPlot_attribute_types = {"epstol": (float, 1)} 
excitonPlot_valid_subtrees = ["exciton", "hole", "electron"] 
excitonPlot_mandatory_attributes = ["electron", "hole"] 
excitonPlot_multiple_children = ["exciton"] 

exciton_attribute_types = {"fix": (str, 1), "lambda": (int, 1)} 

hole_valid_subtrees = ["plot1d", "plot2d", "plot3d"] 

electron_valid_subtrees = ["plot1d", "plot2d", "plot3d"] 

realTimeTDDFT_attribute_types = {"TaylorOrder": (int, 1),
                                 "calcNonlocalCurrentDensity": (bool, 1),
                                 "calculateNExcitedElectrons": (bool, 1),
                                 "calculateTotalEnergy": (bool, 1),
                                 "endTime": (float, 1),
                                 "normalizeWF": (bool, 1),
                                 "printAfterIterations": (int, 1),
                                 "printTimingDetailed": (bool, 1),
                                 "printTimingGeneral": (bool, 1),
                                 "propagator": (str, ["AETRS", "CFM4", "EH", "EHM", "EMR", "RK4", "SE"]),
                                 "subtractJ0": (bool, 1),
                                 "timeStep": (float, 1),
                                 "vectorPotentialSolver": (str, ["euler", "improvedeuler", "midpoint", "rk4"])} 
realTimeTDDFT_valid_subtrees = ["predictorCorrector", "screenshots", "laser", "pmat"] 
realTimeTDDFT_mandatory_attributes = ["pmat"] 

predictorCorrector_attribute_types = {"maxIterations": (int, 1), "tol": (float, 1)} 

screenshots_attribute_types = {"niter": (int, 1)} 
screenshots_valid_subtrees = ["eigenvalues", "projectionCoefficients", "occupations"] 

eigenvalues_attribute_types = {"nEigenvalues": (int, 1), "tolerance": (float, 1)} 

projectionCoefficients_attribute_types = {"format": (str, 1), "printAbsoluteValue": (bool, 1)} 

occupations_attribute_types = {"format": (str, 1)} 

laser_attribute_types = {"fieldType": (str, ["external", "total"])} 
laser_valid_subtrees = ["kick", "trapCos", "sinSq"] 
laser_multiple_children = ["kick", "sinSq", "trapCos"] 

kick_attribute_types = {"amplitude": (float, 1),
                        "direction": (str, ["x", "y", "z"]),
                        "t0": (float, 1),
                        "width": (float, 1)} 

trapCos_attribute_types = {"amplitude": (float, 1),
                           "direction": (str, ["x", "y", "z"]),
                           "omega": (float, 1),
                           "phase": (float, 1),
                           "riseTime": (float, 1),
                           "t0": (float, 1),
                           "width": (float, 1)} 

sinSq_attribute_types = {"amplitude": (float, 1),
                         "direction": (str, ["x", "y", "z"]),
                         "omega": (float, 1),
                         "phase": (float, 1),
                         "pulseLength": (float, 1),
                         "t0": (float, 1)} 

pmat_attribute_types = {"forceHermitian": (bool, 1), "readFromFile": (bool, 1), "writeToFile": (bool, 1)} 

tddft_attribute_types = {"acont": (bool, 1),
                         "ahc": (bool, 1),
                         "alphalrc": (float, 1),
                         "alphalrcdyn": (float, 1),
                         "aresdf": (bool, 1),
                         "aresfxc": (bool, 1),
                         "betalrcdyn": (float, 1),
                         "do": (str, ["fromkernel", "fromscratch"]),
                         "drude": (float, 2),
                         "fxcbsesplit": (float, 1),
                         "fxctype": (str,
                                     ["ALDA", "BO", "BO_SCALAR", "LRCdyn", "LRCdyn_NLF", "LRCstatic", "LRCstatic_NLF",
                                      "MB1", "MB1_NLF", "RBO", "RPA"]),
                         "intraband": (bool, 1),
                         "kerndiag": (bool, 1),
                         "lindhard": (bool, 1),
                         "lmaxalda": (int, 1),
                         "mdfqtype": (int, 1),
                         "nwacont": (int, 1),
                         "torddf": (bool, 1),
                         "tordfxc": (bool, 1)} 

screening_attribute_types = {"do": (str, ["fromscratch", "skip"]),
                             "intraband": (bool, 1),
                             "nempty": (int, 1),
                             "ngridk": (int, 3),
                             "ngridq_interpolation": (int, 3),
                             "nosym": (bool, 1),
                             "nqpt_unique": (int, 1),
                             "qpointsgamma": (bool, 1),
                             "quasiparticle_correction": (str, ["none", "scissor"]),
                             "reducek": (bool, 1),
                             "rgkmax": (float, 1),
                             "screentype": (str, ["diag", "full", "longrange", "noinvdiag"]),
                             "tr": (bool, 1),
                             "vkloff": (float, 3)} 

phonon_screening_attribute_types = {"alat_qe": (float, 1),
                                    "excitation_energy": (float, 1),
                                    "file_type": (str, ["exciting", "quantum_espresso"]),
                                    "phonon_file": (str, 1),
                                    "zstar_file": (str, 1)} 
phonon_screening_mandatory_attributes = ["alat_qe", "excitation_energy", "phonon_file", "zstar_file"] 

expand_eps_attribute_types = {"supercell_1": (int, 3), "supercell_2": (int, 3)} 

BSE_attribute_types = {"aresbse": (bool, 1),
                       "blocks": (str, ["both", "ra", "rr"]),
                       "brixshdf5": (bool, 1),
                       "bsedirsing": (bool, 1),
                       "bsetype": (str, ["IP", "RPA", "singlet", "triplet"]),
                       "checkposdef": (bool, 1),
                       "chibar0": (bool, 1),
                       "chibar0comp": (int, 1),
                       "chibarq": (bool, 1),
                       "coupling": (bool, 1),
                       "cuttype": (str, ["0d", "2d", "none"]),
                       "dichroism": (bool, 1),
                       "distribute": (bool, 1),
                       "econv": (float, 2),
                       "eecs": (int, 1),
                       "efind": (bool, 1),
                       "fbzq": (bool, 1),
                       "iqmtrange": (int, 2),
                       "lmaxdielt": (int, 1),
                       "measure": (bool, 1),
                       "nexc": (int, 1),
                       "ngridksub": (int, 3),
                       "nleblaik": (int, 1),
                       "nosym": (bool, 1),
                       "nosymspec": (bool, 1),
                       "nstlbse": (int, 4),
                       "nstlxas": (int, 2),
                       "outputlevel": (str, ["expert", "normal"]),
                       "readstatetask446": (bool, 1),
                       "reducek": (bool, 1),
                       "rgkmax": (float, 1),
                       "sciavbd": (bool, 1),
                       "sciavqbd": (bool, 1),
                       "sciavqhd": (bool, 1),
                       "sciavqwg": (bool, 1),
                       "sciavtype": (str, ["invscreendiag", "screendiag", "spherical"]),
                       "scrherm": (int, 1),
                       "solver": (str, ["direct", "fastBSE"]),
                       "vkloff": (float, 3),
                       "xas": (bool, 1),
                       "xasatom": (int, 1),
                       "xasedge": (str,
                                   ["K", "L1", "L2", "L23", "L3", "M1", "M2", "M23", "M3", "M4", "M45", "M5", "N1",
                                    "N2", "N23", "N3", "N4", "N45", "N5"]),
                       "xasspecies": (int, 1),
                       "xes": (bool, 1)} 

fastBSE_attribute_types = {"clanczos": (float, 1),
                           "cvtsteplim": (int, 1),
                           "cvttol": (float, 1),
                           "ngridr": (int, 3),
                           "nisdf": (int, 3),
                           "nlanczos": (int, 1),
                           "saveQ": (bool, 1),
                           "seed": (str, ["clock", "fixed"])} 

transitions_valid_subtrees = ["individual", "ranges", "lists"] 

individual_valid_subtrees = ["trans"] 
individual_multiple_children = ["trans"] 

trans_attribute_types = {"action": (str, ["exclude", "include"]),
                         "final": (int, 1),
                         "initial": (int, 1),
                         "kpointnumber": (int, 1)} 

ranges_valid_subtrees = ["range"] 
ranges_multiple_children = ["range"] 

range_attribute_types = {"action": (str, ["exclude", "include"]),
                         "kpointnumber": (int, 1),
                         "start": (int, 1),
                         "statestype": (str, ["finalstates", "initialstates"]),
                         "stop": (int, 1)} 
range_mandatory_attributes = ["statestype"] 

lists_valid_subtrees = ["istate"] 
lists_multiple_children = ["istate"] 

istate_attribute_types = {"action": (str, ["exclude", "include"]),
                          "kpointnumber": (int, 1),
                          "state": (int, 1),
                          "statestype": (str, ["finalstates", "initialstates"])} 
istate_mandatory_attributes = ["statestype"] 

tetra_attribute_types = {"cw1k": (bool, 1),
                         "kordexc": (bool, 1),
                         "qweights": (int, 1),
                         "tetradf": (bool, 1),
                         "tetraocc": (bool, 1)} 

plan_valid_subtrees = ["doonly"] 
plan_multiple_children = ["doonly"] 

doonly_attribute_types = {"task": (str,
                                   ["bse", "bsegenspec", "bsesurvey", "df", "df2", "dielectric", "emattest",
                                    "exccoulint", "excitonWavefunction", "expand_add_eps",
                                    "fastBSE_groundstate_properties", "fastBSE_human_readable_output",
                                    "fastBSE_isdf_cvt", "fastBSE_main", "fxc_alda_check", "idf", "kernxc_bse",
                                    "kernxc_bse3", "phonon_screening", "planewave_elements", "pmatxs2orig",
                                    "portstate(-1)", "portstate(-2)", "portstate(1)", "portstate(2)", "scrcoulint",
                                    "screen", "scrgeneigvec", "scrtetcalccw", "scrwritepmat", "testmain", "testxs",
                                    "tetcalccw", "write_dielectric_matrix", "write_pmat_hdf5_xs", "write_screen",
                                    "write_screened_coulomb", "writebandgapgrid", "writebevec", "writeemat",
                                    "writeematasc", "writekpathweights", "writeoverlapxs", "writepmat", "writepmatasc",
                                    "writepmatxs", "writepwmat", "x0toasc", "x0tobin", "xsestimate", "xsgeneigvec"])} 
doonly_mandatory_attributes = ["task"] 


# gw information 
gw_attribute_types = {"at1": (int, 1),
                      "at2": (int, 1),
                      "coreflag": (str, 1),
                      "debug": (bool, 1),
                      "degeneracyAbsoluteTolerance": (float, 1),
                      "degeneracyRelativeTolerance": (float, 1),
                      "enforceDegeneracy": (bool, 1),
                      "eph": (str, 1),
                      "ibgw": (int, 1),
                      "ibmax": (int, 1),
                      "ibmax2": (int, 1),
                      "ibmin": (int, 1),
                      "ibmin2": (int, 1),
                      "igmax": (int, 1),
                      "igmin": (int, 1),
                      "iik": (int, 1),
                      "jjk": (int, 1),
                      "mblksiz": (int, 1),
                      "nbgw": (int, 1),
                      "nempty": (int, 1),
                      "ngridq": (int, 3),
                      "printSelfC": (bool, 1),
                      "printSpectralFunction": (bool, 1),
                      "qdepw": (str, 1),
                      "reduceq": (bool, 1),
                      "rmax": (float, 1),
                      "rpath": (str, 1),
                      "rpmat": (bool, 1),
                      "skipgnd": (bool, 1),
                      "taskname": (str, 1),
                      "vqloff": (float, 3),
                      "wlo": (float, 1),
                      "wto": (float, 1)} 
gw_valid_subtrees = ["plot1d", "freqgrid", "selfenergy", "mixbasis", "barecoul", "scrcoul", "taskGroup"] 

freqgrid_attribute_types = {"eta": (float, 1),
                            "fconv": (str, 1),
                            "fgrid": (str, 1),
                            "freqmax": (float, 1),
                            "freqmin": (float, 1),
                            "nomeg": (int, 1)} 

selfenergy_attribute_types = {"actype": (str, 1),
                              "eqpsolver": (int, 1),
                              "eshift": (int, 1),
                              "method": (str, 1),
                              "nempty": (int, 1),
                              "singularity": (str, ["crg", "low_dim", "mpb", "none"]),
                              "swidth": (float, 1),
                              "tol": (float, 1)} 
selfenergy_valid_subtrees = ["wgrid"] 

wgrid_attribute_types = {"size": (int, 1), "type": (str, 1), "wmax": (float, 1), "wmin": (float, 1)} 

mixbasis_attribute_types = {"epsmb": (float, 1), "gmb": (float, 1), "lmaxmb": (int, 1)} 

barecoul_attribute_types = {"barcevtol": (float, 1),
                            "basis": (str, 1),
                            "cutofftype": (str, 1),
                            "pwm": (float, 1),
                            "stctol": (float, 1)} 

scrcoul_attribute_types = {"averaging": (str, ["2d", "isotropic"]),
                           "omegap": (float, 1),
                           "q0eps": (float, 3),
                           "scrtype": (str, 1),
                           "subgrid_q0": (int, 3)} 

taskGroup_attribute_types = {"outputFormat": (str, ["binary", "text"])} 
taskGroup_valid_subtrees = ["Coulomb", "epsilon", "invertEpsilon", "sigmac", "sigmax", "vxc"] 

Coulomb_valid_subtrees = ["qpoints"] 
Coulomb_mandatory_attributes = ["qpoints"] 
Coulomb_multiple_children = ["qpoints"] 

epsilon_attribute_types = {"printPolarizabilityFactor": (bool, 1)} 
epsilon_valid_subtrees = ["qpoints"] 
epsilon_mandatory_attributes = ["qpoints"] 
epsilon_multiple_children = ["qpoints"] 

invertEpsilon_valid_subtrees = ["qpoints"] 
invertEpsilon_mandatory_attributes = ["qpoints"] 
invertEpsilon_multiple_children = ["qpoints"] 

sigmac_valid_subtrees = ["kpoints"] 
sigmac_mandatory_attributes = ["kpoints"] 
sigmac_multiple_children = ["kpoints"] 

sigmax_valid_subtrees = ["kpoints"] 
sigmax_mandatory_attributes = ["kpoints"] 
sigmax_multiple_children = ["kpoints"] 

vxc_valid_subtrees = ["kpoints"] 
vxc_mandatory_attributes = ["kpoints"] 
vxc_multiple_children = ["kpoints"] 


# MD information 
MD_attribute_types = {"basisDerivative": (bool, 1),
                      "coreCorrections": (bool, 1),
                      "integrationAlgorithm": (str, ["HeunSimplified"]),
                      "printAllForces": (bool, 1),
                      "timeStep": (float, 1),
                      "type": (str, ["Ehrenfest"]),
                      "updateOverlap": (bool, 1),
                      "updatePmat": (bool, 1),
                      "valenceCorrections": (bool, 1)} 


# eph information 
eph_attribute_types = {"debugeph": (bool, 1),
                       "ibeph": (int, 1),
                       "ibsumeph": (int, 1),
                       "nbeph": (int, 1),
                       "nbsumeph": (int, 1),
                       "nemptyeph": (int, 1),
                       "ngridqeph": (int, 3),
                       "tasknameeph": (str, 1),
                       "vqloffeph": (float, 3)} 
eph_valid_subtrees = ["freqgrideph", "selfenergyeph"] 

freqgrideph_attribute_types = {"freqmaxeph": (float, 1), "nomegeph": (int, 1)} 

selfenergyeph_valid_subtrees = ["SpectralFunctionPloteph"] 

SpectralFunctionPloteph_attribute_types = {"axis": (str, 1),
                                           "eta": (float, 1),
                                           "nwgrid": (int, 1),
                                           "wmax": (float, 1),
                                           "wmin": (float, 1)} 


# valid entries for the xs subtree 'plan'
valid_plan_entries = ["bse", "bsegenspec", "bsesurvey", "df", "df2", "dielectric", "emattest", "exccoulint",
                      "excitonWavefunction", "expand_add_eps", "fastBSE_groundstate_properties",
                      "fastBSE_human_readable_output", "fastBSE_isdf_cvt", "fastBSE_main", "fxc_alda_check", "idf",
                      "kernxc_bse", "kernxc_bse3", "phonon_screening", "planewave_elements", "pmatxs2orig",
                      "portstate(-1)", "portstate(-2)", "portstate(1)", "portstate(2)", "scrcoulint", "screen",
                      "scrgeneigvec", "scrtetcalccw", "scrwritepmat", "testmain", "testxs", "tetcalccw",
                      "write_dielectric_matrix", "write_pmat_hdf5_xs", "write_screen", "write_screened_coulomb",
                      "writebandgapgrid", "writebevec", "writeemat", "writeematasc", "writekpathweights",
                      "writeoverlapxs", "writepmat", "writepmatasc", "writepmatxs", "writepwmat", "x0toasc", "x0tobin",
                      "xsestimate", "xsgeneigvec"] 
