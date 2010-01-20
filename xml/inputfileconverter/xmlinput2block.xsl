<?xml version="1.0" encoding="UTF-8" ?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
  <xsl:output method="text" />
  <!-- usage: xsltproc xmltobloc.xsl input.xml >exciting.in -->
  <xsl:template match="/">
    <xsl:variable name="newline">
      <xsl:text>
</xsl:text>
    </xsl:variable>
    <xsl:text>tasks
</xsl:text>
    <xsl:for-each select="/input/tasks/*">
      <xsl:choose>
        <xsl:when test="name(.)='groundstate'">
          <xsl:choose>
            <xsl:when test="@fromscratch">
              <xsl:text> 0
</xsl:text>
            </xsl:when>
            <xsl:otherwise>
              <xsl:text> 1
</xsl:text>
            </xsl:otherwise>
          </xsl:choose>
        </xsl:when>
        <xsl:when test="name(.)='bandstructure'">
          <xsl:text> 20
</xsl:text>
        </xsl:when>
      </xsl:choose>
    </xsl:for-each>
    <xsl:text>
</xsl:text>
    <xsl:text>avec
</xsl:text>
    <xsl:for-each select="input/structure/crystal/basevect">
      <xsl:text> </xsl:text>
      <xsl:value-of select="." />
      <xsl:value-of select="$newline" />
    </xsl:for-each>
    <xsl:value-of select="$newline" />
    <xsl:if test="/input/structure/species">
      <xsl:text>atoms
 </xsl:text>
      <xsl:value-of select="count(/input/structure/species)" />
      <xsl:value-of select="$newline" />
      <xsl:for-each select="/input/structure/species">
        <xsl:text> '</xsl:text>
        <xsl:value-of select="/input/structure/species/@speciesfile" />
        <xsl:text>'
  </xsl:text>
        <xsl:value-of select="count(atom)" />
        <xsl:value-of select="$newline" />
        <xsl:for-each select="atom">
          <xsl:text>   </xsl:text>
          <xsl:value-of select="@coord" />
          <xsl:text>  </xsl:text>
          <xsl:choose>
            <xsl:when test="@bfcmt">
              <xsl:value-of select="@bfcmt" />
            </xsl:when>
            <xsl:otherwise>
              0.0 0.0 0.0
            </xsl:otherwise>
          </xsl:choose>
          <xsl:value-of select="$newline" />
        </xsl:for-each>
      </xsl:for-each>
      <xsl:value-of select="$newline" />
    </xsl:if>
    <xsl:for-each select="//@*">
      <xsl:choose>
        <xsl:when test="contains(name(.),'xsi:')" />
        <xsl:when test="name(.)='abrev'" />
        <xsl:when test="name(.)='fromscratch'" />
        <xsl:when test="name(.)='atomnr'" />
        <xsl:when test="name(.)='speciesnr'" />
        <xsl:when test="name(.)='grid'" />
        <xsl:when test="name(.)='steps'" />
        <xsl:when test="name(.)='coord'" />
        <xsl:when test="name(.)='bfcmt'" />
        <xsl:when test="name(.)='stretch'">
          <xsl:text>scale1
</xsl:text>
          <xsl:variable name="scale1" select="substring-before(normalize-space(.),' ')" />
          <xsl:value-of select="$scale1" />
          <xsl:text>
					
</xsl:text>
          <xsl:text>scale2
</xsl:text>
          <xsl:variable name="scale2" select="substring-before(substring-after(normalize-space(.),' '),' ')" />
          <xsl:value-of select="$scale2" />
          <xsl:text>
					
</xsl:text>
          <xsl:text>scale3
</xsl:text>
          <xsl:variable name="scale3" select="substring-after(normalize-space(.),concat($scale1,' ',$scale2,' '))" />
          <xsl:value-of select="$scale3" />
          <xsl:text>
					
</xsl:text>
        </xsl:when>
        <xsl:otherwise>
          <xsl:value-of select="name(.)" />
          <xsl:text>
 </xsl:text>
          <xsl:value-of select="." />
          <xsl:value-of select="$newline" />
          <xsl:value-of select="$newline" />
        </xsl:otherwise>
      </xsl:choose>
    </xsl:for-each>
    <xsl:for-each select="//plot1d | //plot2d | //plot3d">
      <xsl:value-of select="name()" />
      <xsl:text>
 </xsl:text>
      <xsl:if test="name()='plot1d'">
        <xsl:value-of select="count(./path/point)" />
        <xsl:text> </xsl:text>
      </xsl:if>
      <xsl:value-of select="./*/@grid | ./*/@steps" />
      <xsl:text>
 </xsl:text>
      <xsl:for-each select="./*/point | ./*/origin">
        <xsl:value-of select="." />
        <xsl:text>
 </xsl:text>
      </xsl:for-each>
      <xsl:text>
</xsl:text>
    </xsl:for-each>
  </xsl:template>
</xsl:stylesheet>