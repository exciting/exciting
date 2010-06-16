<?xml version="1.0" encoding="UTF-8" ?>
<xsl:stylesheet version="1.0"
xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
  <xsl:output method="text" cdata-section-elements="at"/>

<xsl:template match="/">
 # Comments
   object 1 class gridpositions counts <xsl:value-of select="/plot3d/grid/@gridticks"/>
   origin <xsl:value-of select="/plot3d/grid/@origin"/>
   delta <xsl:value-of select="/plot3d/grid/axis[1]/@delta"/>
   delta <xsl:value-of select="/plot3d/grid/axis[2]/@delta"/>
   delta <xsl:value-of select="/plot3d/grid/axis[3]/@delta"/>
   object 2 class gridconnections counts <xsl:value-of select="/plot3d/grid/@gridticks"/>
<xsl:for-each select="/plot3d/function">
   object 3 class array type double rank 0 items  <xsl:value-of select="/plot3d/function/@n"/> data follows
<xsl:for-each select="*/*">
<xsl:value-of select="."/>
</xsl:for-each>
</xsl:for-each>
   attribute "dep" string "positions"
   object "regular positions regular connections" class field
   component "positions" value 1
   component "connections" value 2
   component "data" value 3
   </xsl:template>
   </xsl:stylesheet>