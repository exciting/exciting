<?xml version="1.0" encoding="UTF-8" ?>
    <xsl:stylesheet version="1.0"
      xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
      <xsl:output method="xml" />
      <xsl:template match="/">
        <xsl:comment>
          created from template
        </xsl:comment>
    
    
    <!-- ############# -->
        <input>
          <title><xsl:value-of select="//params/@title"/></title>
          <structure speciespath="./">
            <xsl:copy-of select="/input/structure/crystal" />
            <xsl:copy-of select="/input/structure/species" />
          </structure>
          <groundstate tforce="true" >
          <xsl:copy-of select= "//params/@ngridk" />
          <xsl:copy-of select= "//params/@maxscl" />
          </groundstate>
        </input>
    <!-- ############# -->
    
    
    
      </xsl:template>
    </xsl:stylesheet>
