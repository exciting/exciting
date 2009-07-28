<?xml version="1.0" encoding="UTF-8" ?>
<xsl:stylesheet version="1.0"
	xmlns:xsl="http://www.w3.org/1999/XSL/Transform">	
	<xsl:output method="text" />
	<xsl:template match="/">
		<xsl:for-each select = "//file[@name='input.xml']">
		<xsl:text>fu</xsl:text>
		</xsl:for-each>
		<xsl:for-each select="dirtree/directory/directory/file [@name='input.xml']">
			<xsl:text>cd </xsl:text>
			<xsl:value-of select="path" />
			<xsl:text>/</xsl:text>
			<xsl:variable name="file">
				<xsl:text>./</xsl:text>
				<xsl:value-of select="/dirtree/directory/@name" />
				<xsl:text>/input.xml</xsl:text>
			</xsl:variable>
			<xsl:value-of select="document($file)/info/groundstate/@status" />
		</xsl:for-each>
		<xsl:choose>
		<xsl:when test="//info/groundstate[@status='finished']">
			<xsl:text>calculation finished</xsl:text>
		</xsl:when>
		<xsl:otherwise>
			<xsl:text>failed</xsl:text>
		</xsl:otherwise>
		</xsl:choose>
	</xsl:template>
</xsl:stylesheet>