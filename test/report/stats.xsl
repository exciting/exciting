<?xml version="1.0" encoding="ISO-8859-1"?>
	<!-- Edited by XMLSpy® -->
<html xsl:version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
	xmlns="http://www.w3.org/1999/xhtml">
	<body style="font-family:Arial;font-size:12pt;background-color:#ffffff">
<h1>stats</h1>
		<xsl:for-each select="statistics/run">
			<xsl:sort select="timestamp/time" />
			<div style="margin-bottom:1em;font-size:12pt;background-color:#00ffff">
			<p>
Time:
				<xsl:value-of select="timestamp/@timestring" />
</p><p>
githash: 
				<xsl:value-of select="githash/@hash" />
		</p>	</div>
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
				&amp;chs=500x80&amp;chdl=passed|unspecified|failed&amp;chdlp=t&amp;chl=passed:
				<xsl:value-of select="passed/@count" />
				|unspecified:
				<xsl:value-of select="unspecified/@count" />
				|failed:
				<xsl:value-of select="failed/@count" />
				&amp;chco=006600,f0f000,cc0033
				</xsl:attribute>
				 </img>
		
			</div>	
		</xsl:for-each>
		</body>
</html>
