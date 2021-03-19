import os

import groundStateParser
import propertiesParser
import BSEParser
import GWParser
import RT_TDDFTParser

def parserChooser(path):
    fileName = os.path.split(path)[1]
    if fileName == 'INFO.OUT':
        return groundStateParser.excitingINFO(path)
    elif fileName == 'info.xml':
        return groundStateParser.excitingInfo(path)
    elif fileName == 'atoms.xml':
        return groundStateParser.excitingAtoms(path)
    elif fileName == 'evalcore.xml':
        return groundStateParser.excitingEvalcore(path)
    elif fileName == 'eigval.xml':
        return groundStateParser.excitingEigval(path)
    elif fileName == 'geometry.xml':
        return groundStateParser.excitingGeometry(path)
    elif fileName == 'RHO3D.xml':
        return propertiesParser.excitingPlot3D(path)
    elif fileName == 'VCL3D.xml':
        return propertiesParser.excitingPlot3D(path)
    elif fileName == 'VXC3D.xml':
        return propertiesParser.excitingPlot3D(path)
    elif fileName == 'WF3D.xml':
        return propertiesParser.excitingPlot3D(path)
    elif fileName == 'ELF3D.xml':
        return propertiesParser.excitingPlot3D(path)
    elif fileName == 'EF3D.xml':
        return propertiesParser.excitingPlot3D(path)
    elif fileName == 'LSJ.xml':
        return propertiesParser.excitingLSJ(path)
    elif fileName == 'EFG.xml':
        return propertiesParser.excitingEFG(path)
    elif fileName == 'mossbauer.xml':
        return propertiesParser.excitingMossbauer(path)
    elif fileName == 'expiqr.xml':
        return propertiesParser.excitingExpiqr(path)
    elif fileName == 'effmass.xml':
        return propertiesParser.excitingEffmass(path)
    elif fileName == 'bandstructure.xml':
        return propertiesParser.excitingBandstructure(path)
    elif fileName == 'dos.xml':
        return propertiesParser.excitingDos(path)
    elif fileName == 'KERR.OUT':
        return propertiesParser.excitingKerr(path)
    elif fileName == 'EPSILON_11.OUT':
        return propertiesParser.excitingEpsilon(path)
    elif fileName == 'EPSILON_12.OUT':
        return propertiesParser.excitingEpsilon(path)
    elif fileName == 'EPSILON_33.OUT':
        return propertiesParser.excitingEpsilon(path)
    elif fileName == 'CHI_111.OUT':
        return propertiesParser.excitingChi(path)
    elif fileName == 'ELNES.OUT':
        return propertiesParser.excitingElnes(path)
    elif fileName == 'SEEBECK_11.OUT':
        return propertiesParser.excitingSeebeck(path)
    elif fileName == 'ELECTCOND_11.OUT':
        return propertiesParser.excitingSeebeck(path)
    elif fileName == 'THERMALCOND_11.OUT':
        return propertiesParser.excitingSeebeck(path)
    elif fileName == 'Z_11.OUT':
        return propertiesParser.excitingSeebeck(path)
    elif fileName == 'ldos.out':
        return propertiesParser.excitingLdos(path)
    elif fileName == 'band_edges.out':
        return propertiesParser.excitingBand_edges(path)
    elif fileName == 'spintext.xml':
        return propertiesParser.excitingSpintext(path)
    elif fileName == 'POLARIZATION.OUT':
        return propertiesParser.excitingPolarization(path)
    elif fileName == 'TDOS_WANNIER.OUT':
        return propertiesParser.excitingTdosWannier(path)
    elif fileName == 'WANNIER_INFO.OUT':
        return propertiesParser.excitingWannierInfo(path)
    elif fileName == 'coreoverlap.xml':
        return propertiesParser.coreOverlap(path)
    elif fileName == 'EPSILON_NAR_BSE-singlet-TDA-BAR_SCR-full_OC11.OUT':
        return BSEParser.excitingEPSILON_NAR(path)
    elif fileName == 'EPSILON_NAR_BSE-singlet-TDA-BAR_SCR-full_OC22.OUT':
        return BSEParser.excitingEPSILON_NAR(path)
    elif fileName == 'EPSILON_NAR_BSE-singlet-TDA-BAR_SCR-full_OC33.OUT':
        return BSEParser.excitingEPSILON_NAR(path)
    elif fileName == 'EPSILON_NAR_FXCMB1_OC11_QMT001.OUT':
        return BSEParser.excitingEPSILON_NAR(path)
    elif fileName == 'EPSILON_NAR_FXCMB1_OC22_QMT001.OUT':
        return BSEParser.excitingEPSILON_NAR(path)
    elif fileName == 'EPSILON_NAR_FXCMB1_OC33_QMT001.OUT':
        return BSEParser.excitingEPSILON_NAR(path)
    elif fileName == 'EPSILON_NAR_NLF_FXCMB1_OC11_QMT001.OUT':
        return BSEParser.excitingEPSILON_NAR(path)
    elif fileName == 'EPSILON_NAR_NLF_FXCMB1_OC22_QMT001.OUT':
        return BSEParser.excitingEPSILON_NAR(path)
    elif fileName == 'EPSILON_NAR_NLF_FXCMB1_OC33_QMT001.OUT':
        return BSEParser.excitingEPSILON_NAR(path)
    elif fileName == 'EPSILON_BSE-singlet-TDA-BAR_SCR-full_OC11.OUT':
        return BSEParser.excitingEPSILON_NAR(path)
    elif fileName == 'EPSILON_BSE-singlet-TDA-BAR_SCR-full_OC22.OUT':
        return BSEParser.excitingEPSILON_NAR(path)
    elif fileName == 'EPSILON_BSE-singlet-TDA-BAR_SCR-full_OC33.OUT':
        return BSEParser.excitingEPSILON_NAR(path)
    elif fileName == 'LOSS_NAR_FXCMB1_OC11_QMT001.OUT':
        return BSEParser.excitingLOSS_NAR(path)
    elif fileName == 'LOSS_NAR_FXCMB1_OC22_QMT001.OUT':
        return BSEParser.excitingLOSS_NAR(path)
    elif fileName == 'LOSS_NAR_FXCMB1_OC33_QMT001.OUT':
        return BSEParser.excitingLOSS_NAR(path)
    elif fileName == 'LOSS_NAR_NLF_FXCMB1_OC11_QMT001.OUT':
        return BSEParser.excitingLOSS_NAR(path)
    elif fileName == 'LOSS_NAR_NLF_FXCMB1_OC22_QMT001.OUT':
        return BSEParser.excitingLOSS_NAR(path)
    elif fileName == 'LOSS_NAR_NLF_FXCMB1_OC33_QMT001.OUT':
        return BSEParser.excitingLOSS_NAR(path)
    elif fileName == 'LOSS_BSE-singlet-TDA-BAR_SCR-full_OC11.OUT':
        return BSEParser.excitingLOSS_NAR(path)
    elif fileName == 'LOSS_BSE-singlet-TDA-BAR_SCR-full_OC22.OUT':
        return BSEParser.excitingLOSS_NAR(path)
    elif fileName == 'LOSS_BSE-singlet-TDA-BAR_SCR-full_OC33.out':
        return BSEParser.excitingLOSS_NAR(path)
    elif fileName == 'EXCITON_NAR_BSE-singlet-TDA-BAR_SCR-full_OC11.OUT':
        return BSEParser.excitingEXCITON_NAR_BSE(path)
    elif fileName == 'EXCITON_NAR_BSE-singlet-TDA-BAR_SCR-full_OC22.OUT':
        return BSEParser.excitingEXCITON_NAR_BSE(path)
    elif fileName == 'EXCITON_NAR_BSE-singlet-TDA-BAR_SCR-full_OC33.OUT':
        return BSEParser.excitingEXCITON_NAR_BSE(path)
    elif fileName == 'EXCITON_BSE-singlet-TDA-BAR_SCR-full_OC11.OUT':
        return BSEParser.excitingEXCITON_NAR_BSE(path)
    elif fileName == 'EXCITON_BSE-singlet-TDA-BAR_SCR-full_OC22.OUT':
        return BSEParser.excitingEXCITON_NAR_BSE(path)
    elif fileName == 'EXCITON_BSE-singlet-TDA-BAR_SCR-full_OC33.OUT':
        return BSEParser.excitingEXCITON_NAR_BSE(path)
    elif fileName == 'EFERMI_GW.OUT':
        return GWParser.excitingEFERMI_GW(path)
    elif fileName == 'EVALQP.DAT':
        return GWParser.excitingEVALQP(path)
    elif fileName == 'VXCNN.DAT':
        return GWParser.excitingVXCNN(path)
    elif fileName == 'EPS00_GW.OUT':
        return GWParser.excitingEPS00_GW(path)
    elif fileName == 'JIND.OUT':
        return RT_TDDFTParser.excitingJIND(path)
    elif fileName == 'NEXC.OUT':
        return RT_TDDFTParser.excitingNEXC(path,skiprows=1)
    elif fileName == 'ETOT_RTTDDFT.OUT':
        return RT_TDDFTParser.excitingETOT(path)
    elif 'EIGVAL_' in fileName:
        return RT_TDDFTParser.excitingScreenshots(path)
    elif 'PROJ_' in fileName:
        return RT_TDDFTParser.excitingScreenshots(path,convertTo1D=True)
    else:
        raise NameError("%s does not exist."%fileName)

#   elif(fileName == 'someother.xml'):
#       return groundStateParser.exctingSomeOther(path) and so on
