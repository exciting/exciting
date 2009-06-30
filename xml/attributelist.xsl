<?xml version="1.0" encoding="UTF-8" ?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:xs="http://www.w3.org/2001/XMLSchema">
  <xsl:output method="xml"/>
  <xsl:variable name="newline">
    <xsl:text>
</xsl:text>

  </xsl:variable>
  <xsl:template name="getparents">
    <xsl:param name="separator"/>
    <xsl:for-each select="parent::node()">
      <xsl:call-template name="getparents">
        <xsl:with-param name="separator" select="$separator"/>
      </xsl:call-template>
    </xsl:for-each>
    <xsl:if test="(@maxOccurs='unbounded' or @maxOccurs&gt;1) and $separator='%' and (not (@type))">
      <if test="(parent::node())">
        <xsl:value-of select="$separator"/>
      </if>
      <xsl:value-of select="@name"/>
      <xsl:text>array(</xsl:text>
      <xsl:value-of select="@name"/>
      <xsl:text>index)</xsl:text>
    </xsl:if>
    <xsl:if test="name()='xs:element' and not(@type)">
      <xsl:if test="@name!='input'">
        <xsl:value-of select="$separator"/>
      </xsl:if>
      <xsl:value-of select="@name"/>
    </xsl:if>
  </xsl:template>
  <xsl:template match="/">
    <xsl:element name="attributelist">
      <xsl:value-of select="$newline"/>

      <xsl:for-each select="//xs:attribute|//xs:element[@name!='input']">
        <xsl:element name="attribute">
          <xsl:attribute name="vname">
					<xsl:choose>
					<xsl:when test="xs:annotation/xs:appinfo/oldname">
					<xsl:value-of select="xs:annotation/xs:appinfo/oldname"/>
                    </xsl:when>
					<xsl:otherwise>
					<xsl:value-of select="@name|@ref"/>
					</xsl:otherwise>
					</xsl:choose>
					</xsl:attribute>
          <xsl:attribute name="xpath">
 <xsl:call-template name="getparents"> 
 <xsl:with-param name="separator" select="'/'"/>
  </xsl:call-template>
  <xsl:value-of select="@name|@ref"/>
  </xsl:attribute>
          <xsl:attribute name="fortranname">
					
					<xsl:if test="name(.)='xs:element' and (not( ./@type))">
					<xsl:text>associated(</xsl:text>
					</xsl:if>
 <xsl:call-template name="getparents"> 
 <xsl:with-param name="separator" select="'%'"/>
  </xsl:call-template>
   <xsl:if test="not(name(.)='xs:element'and not(@type))">
   <xsl:text>%</xsl:text>
  <xsl:value-of select="@name|@ref"/>
  </xsl:if>
  <xsl:if test="name(.)='xs:element' and (not( ./@type))">
					<xsl:text>)</xsl:text>
					</xsl:if>
          <xsl:if test="xs:simpleType/xs:restriction/xs:enumeration">
          <xsl:text>number</xsl:text>
          </xsl:if>
                    </xsl:attribute>
          <xsl:choose>
            <xsl:when test="@ref">
              <xsl:attribute name="isref">true</xsl:attribute>
            </xsl:when>
          </xsl:choose>
        </xsl:element>
        <xsl:value-of select="$newline"/>
      </xsl:for-each>
 
    </xsl:element>
  </xsl:template>
</xsl:stylesheet>