<?xml version="1.0" encoding="UTF-8" ?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
  <xsl:template match="/">
    <html>
      <head>
        <title>
          <xsl:value-of select="//title" />
        </title>
      </head>
      <body>
        <table border="0">
          <tr>
            <td>
              <h1>
                <xsl:value-of select="//title" />
              </h1>
            </td>
          </tr>
          <xsl:for-each select="/input/*[name()!='structure' and name()!='title']">
            <tr>
              <td>
                <xsl:call-template name="attributestotable" />
              </td>
            </tr>
          </xsl:for-each>
        </table>
        
        <a href="?_xsl=no">get xml</a>
            <a href="info.xml?_xsl=http://xml.exciting-code.org/info.xsl">Results (infoxml)</a>
      </body>
    </html>
  </xsl:template>
  <xsl:template name="attributestotable">
    <h2>
    <xsl:choose>
    <xsl:when test="name(..)='properties'">
    <xsl:element name="a">
    <xsl:attribute name="href">
    <xsl:text>./</xsl:text>      <xsl:value-of select="name(.)" />
    <xsl:text>.xml</xsl:text>
    </xsl:attribute>
      <xsl:value-of select="name(.)" />
    </xsl:element>
    </xsl:when>
    <xsl:otherwise>
          <xsl:value-of select="name(.)" />
    </xsl:otherwise>
    </xsl:choose>
    </h2>
    <table>
      <xsl:for-each select="@*">
        <tr>
          <td>
            <xsl:element name="a">
            <xsl:attribute name="target">_blank</xsl:attribute>
              <xsl:attribute name="href">
              <xsl:text>http://exciting-code.org/input-reference-essential-expert#att</xsl:text>
              <xsl:value-of select="name(.)" />
	</xsl:attribute>
              <xsl:value-of select="name(.)" />
            </xsl:element>
          </td>
          <td>
            <xsl:value-of select="." />
          </td>
        </tr>
      </xsl:for-each>
      <xsl:for-each select="*">
        <tr>
          <td>
            <xsl:call-template name="attributestotable" />
          </td>
        </tr>
      </xsl:for-each>
    </table>
  </xsl:template>
</xsl:stylesheet>
