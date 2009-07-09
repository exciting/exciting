<?xml version="1.0" encoding="UTF-8" ?>

<xsl:stylesheet version="1.0"
	xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

	<xsl:template name="agraph">
		<xsl:param name="data" />
		<xsl:param name="title" />
		<xsl:element name="img">
			<xsl:attribute name="src"><xsl:text>http://chart.apis.google.com/chart?cht=lc&amp;chs=600x300&amp;chd=t:</xsl:text>
		<xsl:for-each select="$data">
			<xsl:value-of select="." />
			<xsl:if test="not(position()=count($data))">
				<xsl:text>,</xsl:text>
			</xsl:if>
		</xsl:for-each>
		<xsl:variable name="ymin">
			<xsl:for-each select="$data">
				<xsl:sort select="." order="ascending" data-type="number" />
				<xsl:if test="position()=1">
					<xsl:value-of select="." />
				</xsl:if>
			</xsl:for-each>
		</xsl:variable>
		<xsl:variable name="ymax">
			<xsl:for-each select="$data">
				<xsl:sort select="." order="descending" data-type="number" />
				<xsl:if test="position()=1">
					<xsl:value-of select="." />
				</xsl:if>
			</xsl:for-each>
		</xsl:variable>
		<xsl:text>&amp;chds=</xsl:text>
<xsl:value-of select="$ymin"/><xsl:text>,</xsl:text><xsl:value-of select="$ymax"/>
<xsl:text>&amp;chtt=</xsl:text><xsl:value-of select="$title"/>
<xsl:text>&amp;chxt=x,y&amp;chxr=0,1,</xsl:text>
<xsl:value-of select="count($data)"/><xsl:text>|1,</xsl:text>
<xsl:value-of select="$ymin"/><xsl:text>,</xsl:text><xsl:value-of select="$ymax"/>
</xsl:attribute></xsl:element>
</xsl:template>
	<xsl:template match="/">
		<html>
			<body>
				<xsl:call-template name="agraph" select ="/scl/iter/@totalEnergy">
				<xsl:with-param name="data" select="/scl/iter/@totEnergy" />
				<xsl:with-param name="title" select="'Total Energy'" />
				</xsl:call-template>
				<xsl:call-template name="agraph" select ="/scl/iter/@totalEnergy">
				<xsl:with-param name="data" select="/scl/iter/@log10rms" />
				<xsl:with-param name="title" select="'convergence+(powers+of+10)'" />
				</xsl:call-template>

	</body>
	</html>
	
	</xsl:template>
</xsl:stylesheet>