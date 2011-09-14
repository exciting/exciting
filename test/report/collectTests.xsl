<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet 
    xmlns:dt="http://exslt.org/dates-and-times" 
    extension-element-prefixes="dt"
    xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
    version="1.0">
    <xsl:template match="/">
        <reports>
        <report>
            <xsl:attribute name="githash">
                <xsl:value-of select="//githash/@hash"/>
            </xsl:attribute>
            <xsl:attribute name="date">
                <xsl:value-of select="dt:date-time()"/>
            </xsl:attribute>
            <xsl:for-each select="//dir">
                <xsl:variable name="doc" >
                    <xsl:text>../</xsl:text>
                    <xsl:value-of select="."/>
                    <xsl:text>/report.xml</xsl:text>
                    
                </xsl:variable>
                <xsl:copy-of select="document($doc)//test"/>
            </xsl:for-each>
        </report>
             <xsl:copy-of select="document('oldreports.xml')//report"/>
        </reports>
        
    </xsl:template>
</xsl:stylesheet>