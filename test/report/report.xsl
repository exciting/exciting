<?xml version="1.0" encoding="ISO-8859-1"?>
	<!-- Edited by XMLSpy® -->
<html xsl:version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
	xmlns="http://www.w3.org/1999/xhtml">
	<body style="font-family:Arial;font-size:12pt;background-color:#EEEEEE">
		<xsl:for-each select="report/test">
			<xsl:sort select="status"/>
			<xsl:sort select="directory"/>
			<xsl:sort select="name"/>
			<div>
			<xsl:variable name="mystatus">
		<xsl:value-of select="status" /> 
			</xsl:variable>
			
				<xsl:choose>
					<xsl:when test="$mystatus='failed'">
						<xsl:attribute name="style">
				 <xsl:text>background-color:red;color:white;padding:4px</xsl:text>
			   			</xsl:attribute>
					</xsl:when>
					<xsl:when test="$mystatus='passed'">
						<xsl:attribute name="style">
				 <xsl:text>background-color:green;color:white;padding:4px</xsl:text>
			   			</xsl:attribute>
					</xsl:when>
					<xsl:otherwise >
					<xsl:attribute name="style">
				 <xsl:text>background-color:yellow;color:black;padding:4px</xsl:text>
			   			</xsl:attribute>
					</xsl:otherwise>
				</xsl:choose>

				<span style="font-weight:bold">
					<xsl:value-of select="name" />
				</span>
			</div>
			<div style="margin-left:20px;margin-bottom:1em;font-size:10pt">

				<xsl:value-of select="description" />
			</div>
			<div
				style="background-color:#dddddd;margin-left:20px;margin-bottom:1em;font-size:10pt">
				Directory:
				<xsl:value-of select="directory" />
			</div>
			<div
				style="background-color:#dddddd;margin-left:20px;margin-bottom:1em;font-size:10pt">
				<xsl:value-of select="status" /></div>
		</xsl:for-each>
	</body>
</html>
