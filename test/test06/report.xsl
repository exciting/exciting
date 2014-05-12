<?xml version="1.0" encoding="UTF-8" ?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:xs="http://www.w3.org/2001/XMLSchema" xmlns:math="http://exslt.org/math">
  <xsl:output method="xml" indent='yes'/> 
 <xsl:template match="/">
 <report>

  <test>
    <status>
    <xsl:choose>
    <xsl:when test="document('runPbTiO3/atoms.xml')/atomlist/atom[1]/@chemicalSymbol='Pb'"><xsl:text>passed</xsl:text></xsl:when>
    <xsl:otherwise><xsl:text>failed</xsl:text></xsl:otherwise>
    </xsl:choose>
    </status>
    <name>  atom solver works for Pb </name>
    <description>passes if the first atom in atoms.xml is Pb</description>
    <directory>test06/runPbTiO3 </directory>
  </test>
  <test>
    <status>
    <xsl:choose>
    <xsl:when test="math:abs(document('runPbTiO3/atoms.xml')/atomlist/atom[1]/spectrum/state[1]/@energy+ 3.230713949860e3)&lt;1e-8">
     <xsl:text>passed</xsl:text></xsl:when>
     <xsl:otherwise> <xsl:text>failed</xsl:text></xsl:otherwise>
    </xsl:choose>
    </status>
    <name>  Pb 1s energy from the atom solver </name>
    <description>passes if Pb 1s energy is correct:
     <xsl:value-of select="document('runPbTiO3/atoms.xml')/atomlist/atom[1]/spectrum/state[1]/@energy"/> vs. -3.230713949860e3 
    </description>
    <directory>test06/runPbTiO3 </directory>
  </test>
   <test>
    <status>
    <xsl:choose>
    <xsl:when test="math:abs(sum(document('runPbTiO3/atoms.xml')/atomlist/atom[1]/spectrum//state/@energy)+5531.19958104999)&lt;1e-7">
     <xsl:text>passed</xsl:text></xsl:when>
     <xsl:otherwise> <xsl:text>failed</xsl:text></xsl:otherwise>
    </xsl:choose>
    </status>
    <name>  Pb spectrum with the Dirac Hamiltonian </name>
    <description>passes if the simple sum of eigenenergies is correct:
       <xsl:value-of select="sum(document('runPbTiO3/atoms.xml')/atomlist/atom[1]/spectrum//state/@energy)"/> vs. -5531.19958104999 
    </description>
    <directory>test06/runPbTiO3 </directory>
  </test>
  <test>
    <status>
    <xsl:choose>
    <xsl:when test="math:abs(document('runPbTiO3/atoms.xml')/atomlist/atom[1]/NumericalSetup/@rmax - 2.755859979645e1)&lt;1e-6">
     <xsl:text>passed</xsl:text></xsl:when>
     <xsl:otherwise> <xsl:text>failed</xsl:text></xsl:otherwise>
    </xsl:choose>
    </status>
    <name>  correct He species is used </name>
    <description>passes if rmax is the same as in reference:
     <xsl:value-of select="document('runPbTiO3/atoms.xml')/atomlist/atom[1]/NumericalSetup/@rmax"/> vs. 27.55859979645
    </description>
    <directory>test01/runHe </directory>
  </test>

     <test>
    <status>
    <xsl:choose>
    <xsl:when test="math:abs(document('runPbTiO3/info-nonsym.xml')/info/groundstate/scl/iter[last()]/energies/@totalEnergy - document('runPbTiO3/info-sym.xml')/info/groundstate/scl/iter[last()]/energies/@totalEnergy)&lt;1e-5">
     <xsl:text>passed</xsl:text></xsl:when>
     <xsl:otherwise> <xsl:text>failed</xsl:text></xsl:otherwise>
    </xsl:choose>
    </status>
    <name>  Symmetric kinetic energy vs. non-symmetric kinetic energy </name>
    <description>passes if the total energies match:
<xsl:value-of select="document('runPbTiO3/info-nonsym.xml')/info/groundstate/scl/iter[last()]/energies/@totalEnergy"/> vs. <xsl:value-of select="document('runPbTiO3/info-sym.xml')/info/groundstate/scl/iter[last()]/energies/@totalEnergy"/> 
    </description>
    <directory>test06/runPbTiO3 </directory>
  </test>
      <test>
    <status>
    <xsl:choose>
    <xsl:when test="(document('runPbTiO3/info-nonsym.xml')/info/groundstate/scl/iter[last()]/energies/@totalEnergy - document('runPbTiO3/info-sym.xml')/info/groundstate/scl/iter[last()]/energies/@totalEnergy)&lt;1e-5">
     <xsl:text>passed</xsl:text></xsl:when>
     <xsl:otherwise> <xsl:text>failed</xsl:text></xsl:otherwise>
    </xsl:choose>
    </status>
    <name>  Symmetric kinetic energy vs. non-symmetric kinetic energy </name>
    <description>passes if the total energies match:
<xsl:value-of select="document('runPbTiO3/info-nonsym.xml')/info/groundstate/scl/iter[last()]/energies/@totalEnergy"/> vs. <xsl:value-of select="document('runPbTiO3/info-sym.xml')/info/groundstate/scl/iter[last()]/energies/@totalEnergy"/> 
    </description>
    <directory>test06/runPbTiO3 </directory>
  </test>
      <test>
    <status>
    <xsl:choose>
    <xsl:when test="document('runPbTiO3/info-nonsym.xml')/info/groundstate/scl/iter[last()]/@iteration &lt; document('runPbTiO3/input-nonsym.xml')/input/groundstate/@maxscl">
     <xsl:text>passed</xsl:text></xsl:when>
     <xsl:otherwise> <xsl:text>failed</xsl:text></xsl:otherwise>
    </xsl:choose>
    </status>
    <name>  Convergence of the calculation with non-symmetric kinetic energy </name>
    <description>passes scf converges before the threshold given in the input:
<xsl:value-of select="document('runPbTiO3/info-nonsym.xml')/info/groundstate/scl/iter[last()]/@iteration"/> vs. <xsl:value-of select="document('runPbTiO3/input-nonsym.xml')/input/groundstate/@maxscl"/>
    </description>
    <directory>test06/runPbTiO3 </directory>
  </test> 
  <test>
    <status>
    <xsl:choose>
    <xsl:when test="document('runPbTiO3/info-sym.xml')/info/groundstate/scl/iter[last()]/@iteration &lt; document('runPbTiO3/input-sym.xml')/input/groundstate/@maxscl">
     <xsl:text>passed</xsl:text></xsl:when>
     <xsl:otherwise> <xsl:text>failed</xsl:text></xsl:otherwise>
    </xsl:choose>
    </status>
    <name>  Convergence of the calculation with non-symmetric kinetic energy </name>
    <description>passes scf converges before the threshold given in the input:
<xsl:value-of select="document('runPbTiO3/info-sym.xml')/info/groundstate/scl/iter[last()]/@iteration"/> vs. <xsl:value-of select="document('runPbTiO3/input-sym.xml')/input/groundstate/@maxscl"/>
    </description>
    <directory>test06/runPbTiO3 </directory>
  </test>
</report>
</xsl:template>
</xsl:stylesheet> 
