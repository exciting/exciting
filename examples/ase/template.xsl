<?xml version="1.0" encoding="UTF-8" ?>
<xsl:stylesheet version="1.0"
  xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
  <xsl:output method="xml" />
  <xsl:template match="/">
    <xsl:comment>
      created from template.xsl
    </xsl:comment>


<!-- ############# -->
    <input>
      <title></title>
      <structure speciespath="../../species">
        <xsl:copy-of select="/input/structure/crystal" />
        <xsl:copy-of select="/input/structure/species" />
      </structure>
      <groundstate ngridk="4  4  4" maxscl="3"  vkloff="0.5  0.5  0.5" tforce="true" />
    </input>
<!-- ############# -->



  </xsl:template>
</xsl:stylesheet>