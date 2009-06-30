<?xml version="1.0"?>
<xsl:transform xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
  version="1.0" xmlns:xs="http://www.w3.org/2001/XMLSchema">

  <xsl:output method="text"/>

  <xsl:strip-space elements="*"/>
  <xsl:template name="elementtodtd">
  <xsl:text> 
  &lt;!ELEMENT </xsl:text><xsl:value-of select="@name"></xsl:value-of>
  <xsl:if test="*/*/xs:element">
   <xsl:text> (</xsl:text>
   </xsl:if>
    <xsl:if test="not(*/*/xs:element) and not(@type)">
   <xsl:text> EMPTY </xsl:text>
   </xsl:if>
  <xsl:for-each select="*/*/xs:element">
  <xsl:text> </xsl:text>
<xsl:value-of select="@name|@ref"/>
<xsl:choose>
<xsl:when test="@minOccurs=1 and @maxOccurs=1"><xsl:text></xsl:text></xsl:when>
<xsl:when test="@minOccurs >=1"><xsl:text>+</xsl:text></xsl:when>
<xsl:when test="not (./@maxOccurs) or ./@maxOccurs=1"><xsl:text>?</xsl:text></xsl:when>
<xsl:when test="@maxOccurs='unbounded'"><xsl:text>*</xsl:text></xsl:when>
</xsl:choose>
<xsl:if test="not(position()=count(../xs:element))">
   <xsl:text>|</xsl:text>
 </xsl:if>
<xsl:text> </xsl:text>
  </xsl:for-each>
   <xsl:if test="@type">
    <xsl:text> (#PCDATA) </xsl:text>
   </xsl:if>
  <xsl:if test="*/*/xs:element">
   <xsl:text> )*</xsl:text></xsl:if>
    <xsl:text>&gt;</xsl:text>
    <xsl:if test ="*/xs:attribute">
        <xsl:text>
        &lt;!ATTLIST  </xsl:text><xsl:value-of select="@name"/><xsl:text> 
    </xsl:text> 
    <xsl:for-each select="*/xs:attribute">
    <xsl:value-of select="@name|@ref"/><xsl:text> </xsl:text>
        <xsl:value-of select="'CDATA'"/><xsl:text> </xsl:text>
      <xsl:choose>
      <xsl:when test="@default">
     <xsl:value-of select="'#IMPLIED'"/><xsl:text>
    </xsl:text>
      </xsl:when>
      <xsl:when test="@use='required'">
       <xsl:text>#REQUIRED
    </xsl:text>
      </xsl:when>
      <xsl:otherwise>
      <xsl:text>#IMPLIED
    </xsl:text>
      </xsl:otherwise>
      </xsl:choose>

    </xsl:for-each>   
    <xsl:if test="@name='input' or @name='inputset' ">
<xsl:text>xsi:noNamespaceSchemaLocation  CDATA "../../../xml/excitinginput.xsd" 
    xmlns:xsi  CDATA       "http://www.w3.org/2001/XMLSchema-instance"
</xsl:text>
    </xsl:if>
    <xsl:text> &gt;</xsl:text>
    </xsl:if>
  
  </xsl:template>
  
  <xsl:template match="/">
  
  <xsl:for-each select="//xs:element[@name]">
  <xsl:call-template name="elementtodtd"/>
  </xsl:for-each>
  </xsl:template>
  
</xsl:transform>