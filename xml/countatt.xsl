<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
    version="1.0"
    xmlns:sets="http://exslt.org/sets"
    xmlns:xs="http://www.w3.org/2001/XMLSchema"
    >
    <xsl:output  method="xml" indent="yes" />
    <xsl:template match="/">
        <out>
            <xsl:for-each select="sets:distinct(//xs:attribute/@name)">
                <xsl:sort select="."/>
                <xsl:variable name="name" select="."></xsl:variable>
                
                
           
                <att>
                    <xsl:attribute name="occurs">
                        <xsl:value-of select="count(//xs:attribute[@name=$name or @ref=$name])"/>
                    </xsl:attribute>
                    <xsl:attribute name="defs">
                        <xsl:value-of select="count(//xs:attribute[@name=$name ])"/>
                        
                    </xsl:attribute>
                 <xsl:value-of select="$name"/>
               
                
                 
             </att>
            </xsl:for-each>
           
        </out>
    </xsl:template>
  
</xsl:stylesheet>