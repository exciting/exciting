<?xml version="1.0" encoding="UTF-8" ?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:str="http://exslt.org/strings" xmlns:math="http://exslt.org/math"
  xmlns:sets="http://exslt.org/sets">
  <xsl:output method="xml" indent="yes" />
  <xsl:template match="/">
  
    <xsl:element name="input">
 <xsl:text>
</xsl:text>  
      <xsl:element name="structure">
      <xsl:attribute name="speciespath"></xsl:attribute>
 <xsl:text>
 </xsl:text>       
      <xsl:element name="crystal">
      <xsl:text>
  </xsl:text>       
       <xsl:element name="basevect">
     </xsl:element>
     <xsl:text>
  </xsl:text>  
       <xsl:element name="basevect">
     </xsl:element>
     <xsl:text>
  </xsl:text>  
       <xsl:element name="basevect">
     </xsl:element>
     <xsl:text>
 </xsl:text>  
</xsl:element>
 <xsl:text>
 </xsl:text> 
        <xsl:for-each select="sets:distinct(/system/atom/@elementType)">
          <xsl:variable name="v">
            <xsl:value-of select="." />
          </xsl:variable>
          <xsl:element name="species">
            <xsl:attribute name="speciesfile">
              <xsl:value-of select="$v" />
              <xsl:text>.xml</xsl:text>
              </xsl:attribute>
         <xsl:text>
  </xsl:text>              
            <xsl:for-each select="/system/atom[@elementType=$v]">
              <xsl:element name="atom">
                <xsl:attribute name="coord">
<xsl:value-of select="@x3" />
<xsl:text> </xsl:text>
<xsl:value-of select="@z3" />
<xsl:text> </xsl:text>
<xsl:value-of select="@y3" />
<xsl:text> </xsl:text>
</xsl:attribute>
              </xsl:element>
    <xsl:text>
  </xsl:text>          
            </xsl:for-each>
          </xsl:element>
          <xsl:text> 
</xsl:text>
        </xsl:for-each>
      </xsl:element>
  <xsl:text>
</xsl:text>
    </xsl:element>

  </xsl:template>
</xsl:stylesheet>