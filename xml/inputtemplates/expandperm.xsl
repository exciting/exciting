<?xml version="1.0" encoding="UTF-8" ?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:exsl="http://exslt.org/common"
 xmlns:str="http://exslt.org/strings">
	<!--
		documentation :
		http://exciting-code.org/expand-all-parameter-permutations
	-->
	<xsl:output method="xml" encoding="UTF-8" indent="yes" />
	<xsl:template name="expand">
		<xsl:param name="depth" select="1" />
		<xsl:param name="values" select="''" />
		<xsl:param name="path" select="''" />
		<xsl:for-each select="/setup/param[$depth]/val">
			<xsl:variable name="newvalue">
				<xsl:element name="x">
     <xsl:element name="param">
      <xsl:element name="name">
       <xsl:value-of select="normalize-space(/setup/param[$depth]/@name)" />
      </xsl:element>
      <xsl:element name="val">
							<xsl:value-of select="normalize-space(.)" />
      </xsl:element>
     </xsl:element>
					<xsl:for-each select="dep">
      <xsl:element name="param">
       <xsl:element name="name">
        <xsl:value-of select="normalize-space(@name)" />
       </xsl:element>
       <xsl:element name="val">
								<xsl:value-of select="normalize-space(@val)" />
       </xsl:element>
      </xsl:element>
					</xsl:for-each>
				</xsl:element>
			</xsl:variable>
			<xsl:variable name="newpath">
				<xsl:value-of select="normalize-space(/setup/param[$depth]/@name)" />
				<xsl:text>_</xsl:text>
				<xsl:value-of select="normalize-space(.)" />
				<xsl:text>/</xsl:text>
			</xsl:variable>
			<xsl:choose>
				<xsl:when test="$depth &lt; (count(/setup/param))">
					<xsl:call-template name="expand">
						<xsl:with-param name="values">
							<xsl:copy-of select="$values" />
							<xsl:copy-of select="$newvalue" />
						</xsl:with-param>
						<xsl:with-param name="depth" select="$depth+1" />
						<xsl:with-param name="path">
							<xsl:value-of select="$path" />
							<xsl:value-of select="$newpath" />
						</xsl:with-param>
					</xsl:call-template>
				</xsl:when>
				<xsl:otherwise>
					<xsl:element name="set">
      <xsl:for-each select="exsl:node-set($values)/*/param">
       <xsl:attribute name="{name}">
       <xsl:value-of select="val" />
       </xsl:attribute>
						</xsl:for-each>
      <xsl:for-each select="exsl:node-set($newvalue)/*/param">
       <xsl:attribute name="{name}">
       <xsl:value-of select="val" />
       </xsl:attribute>
      </xsl:for-each>
						<xsl:attribute name="path">
          <xsl:value-of select="$path" />
          <xsl:value-of select="$newpath" />
           </xsl:attribute>
					</xsl:element>
				</xsl:otherwise>
			</xsl:choose>
		</xsl:for-each>
	</xsl:template>
	<xsl:template match="/">
		<xsl:element name="experiment">
			<xsl:call-template name="expand" />
		</xsl:element>
	</xsl:template>
</xsl:stylesheet>
