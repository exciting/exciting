<?xml version="1.0" encoding="ISO-8859-1"?>
	<!-- Edited by XMLSpy® -->
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
				http://chart.apis.google.com/
				chart?cht=bhs&amp;chd=t:
				<xsl:value-of select="passed/@percent" />
				|
				<xsl:value-of select="unspecified/@percent" />
				|<xsl:value-of select="failed/@percent" />
				&amp;chs=<xsl:value-of select="all/@count*20"/>x80&amp;chdl=passed|unspecified|failed&amp;chdlp=t&amp;chl=passed:
				<xsl:value-of select="passed/@count" />
				|unspecified:
				<xsl:value-of select="unspecified/@count" />
				|failed:
				<xsl:value-of select="failed/@count" />
				&amp;chco=006600,f0f000,cc0033
				</xsl:attribute>
				 </img>
		
			</div>	
		<div>
		    <span style="margin-right:1em"> <a><xsl:attribute name="href"><xsl:value-of select="@name"/>passed.xml</xsl:attribute>passed</a>  </span>
		    <span style="margin-right:1em"> <a><xsl:attribute name="href"><xsl:value-of select="@name"/>unspecified.xml</xsl:attribute>unspecified</a></span>
		    <span style="margin-right:1em"> <a><xsl:attribute name="href"><xsl:value-of select="@name"/>failed.xml</xsl:attribute>failed</a></span>
		
		
		</div>
		</xsl:for-each>
		</body>
</html>
