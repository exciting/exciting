<?xml version="1.0" encoding="UTF-8" ?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:xs="http://www.w3.org/2001/XMLSchema" xmlns:math="http://exslt.org/math">
  <xsl:output method="xml" indent='yes' />
  <xsl:template match="/">
    <report>
      
      <test>
        <status>
          <xsl:choose>
            <xsl:when test="document('runBSE/LOSS_NAR_BSEsinglet_SCRfull_OC11.OUT.xml')/loss">
              <xsl:text>passed</xsl:text>
            </xsl:when>
            <xsl:otherwise>
              <xsl:text>failed</xsl:text>
            </xsl:otherwise>
          </xsl:choose>
        </status>
        <name> BSE loss funktion works</name>
        <description>pass if loss funktion is written</description>
        <directory>test07/runBSE</directory>
      </test>
          
      <test>
        <status>
          <xsl:choose>
            <xsl:when test="document('runBSE/EPSILON_NAR_BSEsinglet_SCRfull_OC11.OUT.xml')/dielectric">
              <xsl:text>passed</xsl:text>
            </xsl:when>
            <xsl:otherwise>
              <xsl:text>failed</xsl:text>
            </xsl:otherwise>
          </xsl:choose>
        </status>
        <name> BSE dielectric funktion works</name>
        <description>pass if dielectric funktion is written</description>
        <directory>test07/runBSE</directory>
      </test>
       <test>
        <status>
          <xsl:choose>
            <xsl:when test="document('runtddft/LOSS_NAR_FXCMB1_OC11_QMT001.OUT.xml')/loss">
              <xsl:text>passed</xsl:text>
            </xsl:when>
            <xsl:otherwise>
              <xsl:text>failed</xsl:text>
            </xsl:otherwise>
          </xsl:choose>
        </status>
        <name> BSE loss funktion works</name>
        <description>pass if loss funktion is written</description>
        <directory>test07/runBSE</directory>
      </test>
          
      <test>
        <status>
          <xsl:choose>
            <xsl:when test="document('runtddft/EPSILON_NAR_FXCMB1_OC11_QMT001.OUT.xml')/dielectric">
              <xsl:text>passed</xsl:text>
            </xsl:when>
            <xsl:otherwise>
              <xsl:text>failed</xsl:text>
            </xsl:otherwise>
          </xsl:choose>
        </status>
        <name> BSE dielectric funktion works</name>
        <description>pass if dielectric funktion is written</description>
        <directory>test07/runBSE</directory>
      </test>
    </report>
  </xsl:template>
</xsl:stylesheet>