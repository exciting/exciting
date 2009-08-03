<?xml version="1.0" encoding="UTF-8" ?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" >
  <xsl:output method="text"/>
  <xsl:template match="/">
   <xsl:text>#!/bin/sh</xsl:text>
  <xsl:for-each select="//directory">
  <xsl:if test="file/@name='exciting.in'">
  <xsl:text>perl ../xml/blocktoxml.pl .</xsl:text>
  <xsl:value-of select="path"/>
  <xsl:text>/exciting.in > .</xsl:text>
 <xsl:value-of select="path"/>
 <xsl:text>/input.xml
</xsl:text>
  
  </xsl:if>
  </xsl:for-each>
  </xsl:template>
</xsl:stylesheet>