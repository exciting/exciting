<?xml version="1.0" encoding="UTF-8" ?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:xs="http://www.w3.org/2001/XMLSchema" xmlns:math="http://exslt.org/math">
  <xsl:output method="xml" indent='yes'/> 
 <xsl:template match="/">
 <report>

  <test>
    <status>
    <xsl:choose>
    <xsl:when test="document('runHe/atoms.xml')/atomlist/atom/@chemicalSymbol='He'"><xsl:text>passed</xsl:text></xsl:when>
    <xsl:otherwise><xsl:text>failed</xsl:text></xsl:otherwise>
    </xsl:choose>
    </status>
    <name>  atom solver works for He </name>
    <description>passes if atoms.xml exists</description>
    <directory>test01/runHe </directory>
  </test>
  <test>
    <status>
    <xsl:choose>
    <xsl:when test="document('runHe/atoms.xml')/atomlist/atom/@chemicalSymbol='He'"><xsl:text>passed</xsl:text></xsl:when>
    <xsl:otherwise><xsl:text>failed</xsl:text></xsl:otherwise>
    </xsl:choose>
    </status>
    <name>  atom solver works for He </name>
    <description>passes if atoms.xml exists</description>
    <directory>test01/runHe </directory>
  </test>
  <test>
    <status>
    <xsl:choose>
    <xsl:when test="math:abs(document('runHe/atoms.xml')/atomlist/atom/NumericalSetup/@rmax - 10.05929533894)&lt;1e-8">
     <xsl:text>passed</xsl:text></xsl:when>
     <xsl:otherwise> <xsl:text>failed</xsl:text></xsl:otherwise>
    </xsl:choose>
    </status>
    <name>  correct He species is used </name>
    <description>passes if rmax is the same as in reference:
     <xsl:value-of select="document('runHe/atoms.xml')/atomlist/atom/NumericalSetup/@rmax"/> vs. 10.05929533894
    </description>
    <directory>test01/runHe </directory>
  </test>
  <test>
    <status>
    <xsl:choose>
    <xsl:when test="math:abs(document('runHe/atoms.xml')/atomlist/atom/spectrum/state[1]/@energy+ 0.570255918361270)&lt;1e-8">
     <xsl:text>passed</xsl:text></xsl:when>
     <xsl:otherwise> <xsl:text>failed</xsl:text></xsl:otherwise>
    </xsl:choose>
    </status>
    <name>  He 1s energy from the atom solver </name>
    <description>passes if He 1s energy is correct:
     <xsl:value-of select="document('runHe/atoms.xml')/atomlist/atom/spectrum/state[1]/@energy"/> vs. -0.570255918361270
    </description>
    <directory>test01/runHe </directory>
  </test>
  <test>
    <status>
    <xsl:choose>
    <xsl:when test="document('runHe/info.xml')/info/groundstate/@status='finished'"><xsl:text>passed</xsl:text></xsl:when>
    <xsl:otherwise><xsl:text>failed</xsl:text></xsl:otherwise>
    </xsl:choose>
    </status>
    <name>  He groundstate </name>
    <description>passes if groundstate works</description>
    <directory>test01/runHe </directory>
  </test>
  <test>
    <status>
    <xsl:choose>
    <xsl:when test="math:abs(document('runHe/info.xml')/info/groundstate/scl/iter[last()]/energies/@totalEnergy+ 2.834836)&lt;0.000001">
      <xsl:text>passed</xsl:text></xsl:when>
      <xsl:otherwise><xsl:text>failed</xsl:text></xsl:otherwise>
    </xsl:choose>
    </status>
    <name>  SCF energy  </name>
    <description>passes if the total energy is correct within 1d-6:
    <xsl:value-of select="document('runHe/info.xml')/info/groundstate/scl/iter[last()]/energies/@totalEnergy"/> vs. -2.834836
    </description>
    <directory>test01/runHe </directory>
  </test>
   <test>
    <status>
    <xsl:choose>
    <xsl:when test="document('runAr/atoms-default.xml')/atomlist/atom/@chemicalSymbol='Ar'"><xsl:text>passed</xsl:text></xsl:when>
    <xsl:otherwise><xsl:text>failed</xsl:text></xsl:otherwise>
    </xsl:choose>
    </status>
    <name>  atom solver works for Ar with the Dirac Hamiltonian</name>
    <description>passes if atoms.xml exists</description>
    <directory>test01/runAr </directory>
  </test>
  <test>
    <status>
    <xsl:choose>
    <xsl:when test="document('runAr/atoms-nr.xml')/atomlist/atom/@chemicalSymbol='Ar'"><xsl:text>passed</xsl:text></xsl:when>
    <xsl:otherwise><xsl:text>failed</xsl:text></xsl:otherwise>
    </xsl:choose>
    </status>
    <name>  atom solver works for Ar with the non-relativistic Hamiltonian</name>
    <description>passes if atoms-nr.xml exists</description>
    <directory>test01/runAr </directory>
  </test>
   <test>
    <status>
    <xsl:choose>
    <xsl:when test="document('runAr/atoms-zora.xml')/atomlist/atom/@chemicalSymbol='Ar'"><xsl:text>passed</xsl:text></xsl:when>
    <xsl:otherwise><xsl:text>failed</xsl:text></xsl:otherwise>
    </xsl:choose>
    </status>
    <name>  atom solver works for Ar with ZORA </name>
    <description>passes if atoms-zora.xml exists</description>
    <directory>test01/runAr </directory>
  </test>
   <test>
    <status>
    <xsl:choose>
    <xsl:when test="document('runAr/atoms-iora.xml')/atomlist/atom/@chemicalSymbol='Ar'"><xsl:text>passed</xsl:text></xsl:when>
    <xsl:otherwise><xsl:text>failed</xsl:text></xsl:otherwise>
    </xsl:choose>
    </status>
    <name>  atom solver works for Ar with IORA</name>
    <description>passes if atoms-iora.xml exists</description>
    <directory>test01/runAr </directory>
  </test>
rmax="1.745072029189e1"
   <test>
    <status>
    <xsl:choose>
    <xsl:when test="math:abs(document('runAr/atoms-default.xml')/atomlist/atom/NumericalSetup/@rmax - 17.45072029189)&lt;1e-6">
     <xsl:text>passed</xsl:text></xsl:when>
     <xsl:otherwise> <xsl:text>failed</xsl:text></xsl:otherwise>
    </xsl:choose>
    </status>
    <name>  correct Ar species is used </name>
    <description>passes if rmax is the same as in reference:
       <xsl:value-of select="document('runAr/atoms-default.xml')/atomlist/atom/NumericalSetup/@rmax"/> vs. 17.45072029189
    </description>
    <directory>test01/runAr </directory>
  </test>
   <test>
    <status>
    <xsl:choose>
    <xsl:when test="math:abs(sum(document('runAr/atoms-default.xml')/atomlist/atom/spectrum//state/@energy)+143.741067299331)&lt;1e-8">
     <xsl:text>passed</xsl:text></xsl:when>
     <xsl:otherwise> <xsl:text>failed</xsl:text></xsl:otherwise>
    </xsl:choose>
    </status>
    <name>  Ar spectrum with the Dirac Hamiltonian </name>
    <description>passes if the simple sum of eigenenergies is correct:
       <xsl:value-of select="sum(document('runAr/atoms-default.xml')/atomlist/atom/spectrum//state/@energy)"/> vs. -143.741067299331
    </description>
    <directory>test01/runAr </directory>
  </test>
    <test>
    <status>
    <xsl:choose>
    <xsl:when test="math:abs(sum(document('runAr/atoms-nr.xml')/atomlist/atom/spectrum//state/@energy)+ 143.128392441857)&lt;1e-8">
     <xsl:text>passed</xsl:text></xsl:when>
     <xsl:otherwise> <xsl:text>failed</xsl:text></xsl:otherwise>
    </xsl:choose>
    </status>
    <name>  Ar spectrum with the non-relativistic Hamiltonian </name>
    <description>passes if the simple sum of eigenenergies is correct:
       <xsl:value-of select="sum(document('runAr/atoms-nr.xml')/atomlist/atom/spectrum//state/@energy)"/> vs. -143.128392441857
    </description>
    <directory>test01/runAr </directory>
  </test>
    <test>
    <status>
    <xsl:choose>
    <xsl:when test="math:abs(sum(document('runAr/atoms-zora.xml')/atomlist/atom/spectrum//state/@energy)+144.17745033111)&lt;1e-8">
     <xsl:text>passed</xsl:text></xsl:when>
     <xsl:otherwise> <xsl:text>failed</xsl:text></xsl:otherwise>
    </xsl:choose>
    </status>
    <name>  Ar spectrum with ZORA </name>
    <description>passes if the simple sum of eigenenergies is correct:
       <xsl:value-of select="sum(document('runAr/atoms-zora.xml')/atomlist/atom/spectrum//state/@energy)"/> vs. -144.177450331111
    </description>
    <directory>test01/runAr </directory>
  </test>
   <test>
    <status>
     <xsl:choose>
     <xsl:when test="math:abs(sum(document('runAr/atoms-iora.xml')/atomlist/atom/spectrum//state/@energy)+143.73151277489)&lt;1e-8">
     <xsl:text>passed</xsl:text></xsl:when>
     <xsl:otherwise> <xsl:text>failed</xsl:text></xsl:otherwise>
    </xsl:choose>
    </status>
    <name>  Ar spectrum with IORA </name>
    <description>passes if the simple sum of eigenenergies is correct:
       <xsl:value-of select="sum(document('runAr/atoms-iora.xml')/atomlist/atom/spectrum//state/@energy)"/> vs. -143.73151277489
    </description>
    <directory>test01/runAr </directory>
  </test>
    <test>
    <status>
    <xsl:choose>
    <xsl:when test="math:abs(document('runAr/info-default.xml')/info/groundstate/scl/iter[last()]/energies/@totalEnergy+527.817961012)&lt;1e-7">
     <xsl:text>passed</xsl:text></xsl:when>
     <xsl:otherwise> <xsl:text>failed</xsl:text></xsl:otherwise>
    </xsl:choose>
    </status>
    <name>  Total energy of Ar with the default (Dirac+ZORA) Hamiltonian </name>
    <description>passes if the total energy is correct:
<xsl:value-of select="document('runAr/info-default.xml')/info/groundstate/scl/iter[last()]/energies/@totalEnergy"/> vs. -527.817961012
    </description>
    <directory>test01/runAr </directory>
  </test>
    <test>
    <status>
    <xsl:choose>
    <xsl:when test="math:abs(document('runAr/info-nr.xml')/info/groundstate/scl/iter[last()]/energies/@totalEnergy+525.946157347)&lt;1e-7">
     <xsl:text>passed</xsl:text></xsl:when>
     <xsl:otherwise> <xsl:text>failed</xsl:text></xsl:otherwise>
    </xsl:choose>
    </status>
    <name>  Total energy of Ar with the non-relativistic Hamiltonian </name>
    <description>passes if the total energy is correct:
<xsl:value-of select="document('runAr/info-nr.xml')/info/groundstate/scl/iter[last()]/energies/@totalEnergy"/> vs. -525.946157347
    </description>
    <directory>test01/runAr </directory>
  </test>
    <test>
    <status>
    <xsl:choose>
    <xsl:when test="math:abs(document('runAr/info-zora.xml')/info/groundstate/scl/iter[last()]/energies/@totalEnergy+528.802397729)&lt;1e-7">
     <xsl:text>passed</xsl:text></xsl:when>
     <xsl:otherwise> <xsl:text>failed</xsl:text></xsl:otherwise>
    </xsl:choose>
    </status>
    <name>  Total energy of Ar with ZORA </name>
    <description>passes if the total energy is correct:
<xsl:value-of select="document('runAr/info-zora.xml')/info/groundstate/scl/iter[last()]/energies/@totalEnergy"/> vs. -528.802397729
    </description>
    <directory>test01/runAr </directory>
  </test>
    <test>
    <status>
    <xsl:choose>
    <xsl:when test="math:abs(document('runAr/info-iora.xml')/info/groundstate/scl/iter[last()]/energies/@totalEnergy+527.818391513)&lt;1e-7">
     <xsl:text>passed</xsl:text></xsl:when>
     <xsl:otherwise> <xsl:text>failed</xsl:text></xsl:otherwise>
    </xsl:choose>
    </status>
    <name>  Total energy of Ar with the IORA </name>
    <description>passes if the total energy is correct:
<xsl:value-of select="document('runAr/info-iora.xml')/info/groundstate/scl/iter[last()]/energies/@totalEnergy"/> vs. -527.818391513
    </description>
    <directory>test01/runAr </directory>
  </test>
  
</report>
</xsl:template>
</xsl:stylesheet> 
