<?xml version="1.0" encoding="ISO-8859-1"?>
<xsl:transform xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
  version="1.0" xmlns:xs="http://www.w3.org/2001/XMLSchema">

<xsl:template match="/">
<xsl:call-template name="stats2html"></xsl:call-template>
</xsl:template>
<xsl:template name="stats2html">
<xsl:param name="urlpf"></xsl:param>
<html xsl:version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
	xmlns="http://www.w3.org/1999/xhtml">
	<body style="font-family:Arial;font-size:12pt;background-color:#ffffff">
<h1>Statistics </h1>
		<xsl:for-each select="statistics/run">
			<xsl:sort select="timestamp/time" />
			<div style="margin-bottom:1em;margin-top:2em;font-size:12pt;background-color:#dddddd">
			<span style="font-size:18pt;margin-right:1em">Test run Nr:<xsl:value-of select="@name"/> </span>
<span style="margin-right:2em"><b>Time: </b> 
				<xsl:value-of select="timestamp/@timestring" />
</span>
<span style="margin-right:2em"><b>githash: </b>
 
				<xsl:value-of select="githash/@hash" />
</span>
		</div>
			<div
				style="background-color:#ffffff;margin-bottom:1em;font-size:10pt">
				
			
				<img>
				<xsl:attribute name="src">
				<xsl:text>http://chart.apis.google.com/chart?cht=bhs&amp;chd=t:</xsl:text>
				<xsl:value-of select="passed/@percent" />
				<xsl:text>|</xsl:text>
				<xsl:value-of select="unspecified/@percent" />
				<xsl:text>|</xsl:text><xsl:value-of select="failed/@percent" />
				<xsl:text>&amp;chs=</xsl:text><xsl:value-of select="all/@count*20"/>
                <xsl:text>x80&amp;chdl=passed|unspecified|failed&amp;chdlp=t&amp;chl=passed:</xsl:text>
				<xsl:value-of select="passed/@count" />
				<xsl:text>|unspecified:</xsl:text>
				<xsl:value-of select="unspecified/@count" />
				<xsl:text>|failed:</xsl:text>
				<xsl:value-of select="failed/@count" />
				<xsl:text>&amp;chco=006600,f0f000,cc0033</xsl:text>
				</xsl:attribute>
				 </img>
		
			</div>	
		<div>
		    <span style="margin-right:1em"> <a><xsl:attribute name="href"><xsl:value-of select="$urlpf"/><xsl:value-of select="@name"/>passed.xml</xsl:attribute>passed   <xsl:value-of select="passed/@count" /></a>  </span>
		    <span style="margin-right:1em"> <a><xsl:attribute name="href"><xsl:value-of select="$urlpf"/><xsl:value-of select="@name"/>unspecified.xml</xsl:attribute>unspecified   <xsl:value-of select="unspecified/@count" /></a></span>
		    <span style="margin-right:1em"> <a><xsl:attribute name="href"><xsl:value-of select="$urlpf"/><xsl:value-of select="@name"/>failed.xml</xsl:attribute>failed   <xsl:value-of select="failed/@count" /></a></span>
		
		
		</div>
		</xsl:for-each>
		</body>
</html>
</xsl:template>
</xsl:transform>
