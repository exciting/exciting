<?xml version="1.0" encoding="UTF-8" ?>
<xsl:stylesheet version="1.0"
	xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:xs="http://www.w3.org/2001/XMLSchema">
	<xsl:output method="text" />
	<xsl:variable name="newline">
<xsl:text>
</xsl:text>
	</xsl:variable>
 <xsl:template match="displaymath">
<xsl:text> 
[[math label]] 
</xsl:text> <xsl:value-of select="."></xsl:value-of><xsl:text>
[[/math]]
</xsl:text>
  </xsl:template>
   <xsl:template match="inlinemath">
 <xsl:text> [[$ </xsl:text> <xsl:value-of select="normalize-space(.)"></xsl:value-of><xsl:text> $]]</xsl:text>
  </xsl:template>
    <xsl:template match="pre">
 <xsl:text> {{ </xsl:text> <xsl:value-of select="normalize-space(.)"></xsl:value-of><xsl:text> }}</xsl:text>
  </xsl:template>
  
    <xsl:template match="it">
 <xsl:text> // </xsl:text> <xsl:value-of select="normalize-space(.)"></xsl:value-of><xsl:text> //</xsl:text>
  </xsl:template>
   <xsl:template match="bf">
 <xsl:text> **</xsl:text> <xsl:value-of select="normalize-space(.)"></xsl:value-of><xsl:text>** </xsl:text>
  </xsl:template>
   <xsl:template match="text()">
 <xsl:value-of select="normalize-space(.)"/>
  </xsl:template>
  
   <xsl:template match="xs:documentation">
<xsl:apply-templates select="text()|inlinemath|displaymath|pre|it"/>
  </xsl:template>
	<xsl:template name="atributdescriptions" match="none">
	
            
		<!--  <xsl:value-of select="./xs:complexType/xs:annotation/xs:documentation" />-->

  <xsl:apply-templates select="xs:annotation/xs:documentation"/>
		<xsl:value-of select="$newline" />
		<xsl:if test="./xs:complexType/xs:attribute|./xs:attribute">
      <xsl:text>++++* Attributes:
</xsl:text>
			<xsl:for-each select="./xs:complexType/xs:attribute|./xs:attribute">
				<xsl:text>+++++* @</xsl:text>
				<xsl:value-of select="./@name" />
				<xsl:value-of select="./@ref" />
				<xsl:text>
</xsl:text>
				<xsl:if test="./@ref">
				<xsl:apply-templates select="//*[@name=./@ref]/*/xs:documentation"/>
					<xsl:text> [[#</xsl:text>
					<xsl:value-of select="./@ref" />
					<xsl:text> see]] </xsl:text>
					
				</xsl:if>
				<xsl:if test="./@name">
					<xsl:text>
[[# </xsl:text>
					<xsl:value-of select="./@name" />
					<xsl:text>]]</xsl:text>
				</xsl:if>

				<xsl:apply-templates select="xs:annotation/xs:documentation"/>
				
				<xsl:value-of select="$newline" />
						<xsl:if test="./xs:simpleType/xs:restriction/xs:enumeration/@value|./@type">	
					<xsl:text>
</xsl:text>

	<xsl:text>**type:** </xsl:text>

<xsl:if test="./xs:simpleType/xs:restriction/xs:enumeration/@value">
<xsl:text> select:
</xsl:text>
<xsl:for-each select="./xs:simpleType/xs:restriction/xs:enumeration/@value  ">
<xsl:text>
* </xsl:text><xsl:value-of select="."/>
</xsl:for-each>
<xsl:text> </xsl:text>
</xsl:if>
					
					<xsl:value-of select="./@type" />

					<xsl:if test="./@use">
						<xsl:text> **use:** </xsl:text>
						<xsl:value-of select="./@use" />
					</xsl:if>
					<xsl:if test="./@default">
						<xsl:text> **default-value:** </xsl:text>
						<xsl:value-of select="./@default" />
					</xsl:if>
					<xsl:text> 
</xsl:text>
			</xsl:if>
			</xsl:for-each>
			<xsl:for-each select="./*/xs:attributeGroup|./xs:attributeGroup">
				<xsl:if test="./@ref">
					<xsl:text>
+++++* attribute group 
</xsl:text>
					<xsl:value-of select="./@ref" />
					<xsl:apply-templates select="//*[@name=./@ref]/*/xs:documentation"/>
					<xsl:text>
[[#</xsl:text><xsl:value-of select="./@ref" />
					<xsl:text> See]]</xsl:text>
				 
          <xsl:apply-templates select="xs:annotation/xs:documentation"/>
				</xsl:if>
			</xsl:for-each>
			<xsl:text>
</xsl:text>
		</xsl:if>
	</xsl:template>


	<xsl:template name="section">
		<xsl:param name="depth" />
		<xsl:choose>
			<xsl:when test="$depth=0">
				<xsl:text>+ </xsl:text>
			</xsl:when>
			<xsl:when test="$depth=1">
				<xsl:text>++ </xsl:text>
			</xsl:when>
			<xsl:when test="$depth=2">
				<xsl:text>+++ </xsl:text>
			</xsl:when>
			<xsl:when test="$depth=3">
				<xsl:text>++++* </xsl:text>
			</xsl:when>
			<xsl:when test="$depth=4">
				<xsl:text>+++++* </xsl:text>
			</xsl:when>
		</xsl:choose>
    <xsl:text>&lt;</xsl:text>
		<xsl:value-of select="./@name" />
		<xsl:value-of select="./@ref" />
		<xsl:choose>
			<xsl:when test="name(.)='xs:attributeGroup'">
				<xsl:text> attribute group</xsl:text>
			</xsl:when>
			<xsl:otherwise>
				<xsl:text>&gt;</xsl:text>
			</xsl:otherwise>
		</xsl:choose>
		<xsl:if test="./@ref">
			<xsl:text>
[[#</xsl:text>
			<xsl:value-of select="./@ref" />
			<xsl:text> See]]</xsl:text>
			<xsl:apply-templates select="//*[@name=./@ref]/*/xs:documentation"/>
			
		</xsl:if>
		<xsl:if test="./@name">
			<xsl:text>
[[# </xsl:text>
			<xsl:value-of select="./@name" />
			<xsl:text>]]</xsl:text>
		</xsl:if>
		<xsl:call-template name="atributdescriptions" />
		<xsl:for-each select="./*/*/xs:element|./*/*/xs:group|./*/xs:group">
			<xsl:call-template name="section">
				<xsl:with-param name="depth" select="$depth+1" />
			</xsl:call-template>
		</xsl:for-each>
	</xsl:template>
	<xsl:template match="/">
		<xsl:text>
</xsl:text>
		<xsl:text>
[[f>toc]]
+ &lt;Root&gt;  </xsl:text>

		<xsl:value-of select="/xs:schema/xs:annotation/xs:appinfo/root"/>
		<xsl:text>
[[# input]]</xsl:text>

		<xsl:for-each select="/*/xs:element[@name=/xs:schema/xs:annotation/xs:appinfo/root]">
			<xsl:call-template name="atributdescriptions" />
		</xsl:for-each>
		<xsl:for-each select="/*/xs:element[@name=/xs:schema/xs:annotation/xs:appinfo/root]/*/*/xs:element">
			<xsl:call-template name="section">
				<xsl:with-param name="depth" select="0" />
			</xsl:call-template>
		</xsl:for-each>
+ reused Elements
		These elements make sense in more than only one context. In this documentation there are references placed if one of these applies.
		<xsl:for-each select="/*/xs:element[@name!=/xs:schema/xs:annotation/xs:appinfo/root]|/*/xs:group">
			<xsl:call-template name="section">
				<xsl:with-param name="depth" select="1" />
			</xsl:call-template>
		</xsl:for-each>
		<xsl:text>
+ reused attributes</xsl:text>
These attributes make sense in more than only one context.  In this documentation there are references placed if one of these applies.
	<xsl:for-each select="/*">
			<xsl:call-template name="atributdescriptions" />
		</xsl:for-each>

		<xsl:for-each select="/*/xs:attributeGroup">
			<xsl:call-template name="section">
				<xsl:with-param name="depth" select="1" />
			</xsl:call-template>
		</xsl:for-each>
		<xsl:text>
</xsl:text>
	</xsl:template>
 
</xsl:stylesheet>