<?xml version="1.0" encoding="UTF-8" ?>
<xsl:stylesheet version="1.0"
xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
 <xsl:output method="text" cdata-section-elements="at"/>
<xsl:template match="//point">
<xsl:for-each select="@*|band/@eval|bc/@character">
<xsl:value-of select="."/><xsl:text>  </xsl:text>
</xsl:for-each>
<xsl:text>
</xsl:text>
</xsl:template>
<xsl:template match="//text()">
<!-- remove all whitespace -->
</xsl:template>
<xsl:template match="//band">
<xsl:apply-templates/>
<xsl:text>
</xsl:text>
</xsl:template>
</xsl:stylesheet>
	