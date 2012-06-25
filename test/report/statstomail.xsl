<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:transform xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
  version="1.0" xmlns:xs="http://www.w3.org/2001/XMLSchema">

  <xsl:import href="reports2html.xsl"/>

   <xsl:output  method="xml" omit-xml-declaration="yes"/>
 <xsl:template match ="/">From:tester@g44228@unileoben.ac.at
To:christian.meisenbichler@mu-leoben.at,stephan.sagmeister@unileoben.ac.at
Subject:exciting master <xsl:value-of select="statistics/run[1]/timestamp/@timestring "/>
<xsl:text> failed:</xsl:text>
<xsl:value-of select="/statistics/run[1]/failed/@count "/>
<xsl:text>
MIME-Version: 1.0
Content-Type: text/html; charset=ISO-8859-1
</xsl:text>
  <xsl:call-template name="reports2html">
  <xsl:with-param name="urlpf" >http://g44222/tests/master/</xsl:with-param></xsl:call-template>
  </xsl:template>


</xsl:transform>
