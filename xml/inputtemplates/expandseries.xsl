<?xml version="1.0" encoding="UTF-8" ?>
 
<xsl:stylesheet version="1.0"
  xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:exsl="http://exslt.org/common">
  <!--
    documentation :
    http://exciting-code.org/series-expansion-template
  -->
  <xsl:output method="xml" encoding="UTF-8" indent="yes" />
  
  <xsl:template match="/">
    <xsl:element name="setup">
      <xsl:for-each select="/setup/param">
        <xsl:choose>
          <xsl:when test="val">
            <xsl:copy-of select="." />
          </xsl:when>
          <xsl:when test="series">
            <xsl:element name="param">
              <xsl:attribute name="name">
           <xsl:value-of select="@name" />
           </xsl:attribute>
              <xsl:call-template name="expandseries">
                <xsl:with-param name="stop">
                  <xsl:value-of select="series/@stop" />
                </xsl:with-param>
                <xsl:with-param name="increment">
                  <xsl:value-of select="series/@increment" />
                </xsl:with-param>
                <xsl:with-param name="value">
                  <xsl:value-of select="series/@start" />
                </xsl:with-param>
              </xsl:call-template>
            </xsl:element>
          </xsl:when>
          <xsl:when test="geomseries">
            <xsl:element name="param">
              <xsl:attribute name="name">
           <xsl:value-of select="@name" />
           </xsl:attribute>
              <xsl:call-template name="expandgeometricseries">
                <xsl:with-param name="stop">
                  <xsl:value-of select="geomseries/@stop" />
                </xsl:with-param>
                <xsl:with-param name="factor">
                  <xsl:value-of select="geomseries/@factor" />
                </xsl:with-param>
                <xsl:with-param name="value">
                  <xsl:value-of select="geomseries/@start" />
                </xsl:with-param>
              </xsl:call-template>
            </xsl:element>
          </xsl:when>
        </xsl:choose>
 
      </xsl:for-each>
    </xsl:element>
 
  </xsl:template>
  <xsl:template name="expandseries">
    <xsl:param name="stop" />
    <xsl:param name="increment" />
    <xsl:param name="value" />
 
    <xsl:element name="val">
      <xsl:value-of select="$value" />
    </xsl:element>
    <xsl:if test="($increment &gt; 0 and $value &lt; $stop) or
($increment &lt; 0 and $value &gt; $stop)     ">
      <xsl:call-template name="expandseries">
        <xsl:with-param name="stop">
          <xsl:value-of select="$stop" />
        </xsl:with-param>
        <xsl:with-param name="increment">
          <xsl:value-of select="$increment" />
        </xsl:with-param>
        <xsl:with-param name="value">
          <xsl:value-of select="$value + $increment" />
        </xsl:with-param>
      </xsl:call-template>
    </xsl:if>
  </xsl:template>
  <xsl:template name="expandgeometricseries">
    <xsl:param name="stop" />
    <xsl:param name="factor" />
    <xsl:param name="value" />
 
    <xsl:element name="val">
      <xsl:value-of select="$value" />
    </xsl:element>
    <xsl:if test="($factor &gt; 1 and $value &lt; $stop) or
($factor &lt; 1 and $value &gt; $stop)     ">
      <xsl:call-template name="expandgeometricseries">
        <xsl:with-param name="stop">
          <xsl:value-of select="$stop" />
        </xsl:with-param>
        <xsl:with-param name="factor">
          <xsl:value-of select="$factor" />
        </xsl:with-param>
        <xsl:with-param name="value">
          <xsl:value-of select="$value * $factor" />
        </xsl:with-param>
      </xsl:call-template>
    </xsl:if>
  </xsl:template>
</xsl:stylesheet>