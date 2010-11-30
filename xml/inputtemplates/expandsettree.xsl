<?xml version="1.0" encoding="UTF-8" ?>
<xsl:stylesheet version="1.0"
	xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:exsl="http://exslt.org/common">
	<!-- documentation : http://exciting-code.org/expand-all-parameter-permutations -->
	<xsl:output method="xml" encoding="UTF-8" indent="yes" />
	<xsl:variable name="newline">
		<xsl:text>
</xsl:text>
	</xsl:variable>
	<!-- //////////////////////////////////////////////////////////////////////////////// -->
	<xsl:template name="crawltree">
		<xsl:param name="depth" select="1" />
		<xsl:param name="values" select="''" />
		<xsl:param name="path" select="''" />
		<xsl:for-each select="/setup/param[$depth]/val">
			<xsl:variable name="newpath">
				<xsl:value-of select="normalize-space(/setup/param[$depth]/@name)" />
				<xsl:text>_</xsl:text>
				<xsl:value-of select="normalize-space(.)" />
				<xsl:text>/</xsl:text>
			</xsl:variable>
			<xsl:choose>
				<xsl:when test="$depth &lt; (count(/setup/param))">
					<!-- write subsetup.xml files on each node of the tree but the innermost -->
					<xsl:document href="{concat($path,$newpath)}/subexperiment.xml"
						method="xml" encoding="UTF-8" indent="yes">
						<xsl:element name="experiment">
							<xsl:value-of select="$newline" />
							<xsl:call-template name="expand">
								<xsl:with-param name="depth" select="$depth+1" />
							</xsl:call-template>
						</xsl:element>
					</xsl:document>
					<xsl:call-template name="crawltree">
						<xsl:with-param name="depth" select="$depth+1" />
						<xsl:with-param name="path">
							<xsl:value-of select="$path" />
							<xsl:value-of select="$newpath" />
						</xsl:with-param>
					</xsl:call-template>
				</xsl:when>
			</xsl:choose>
		</xsl:for-each>
	</xsl:template>
	<!-- //////////////////////////////////////////////////////////////////////////////// -->
	<xsl:template name="expand">
		<xsl:param name="depth" select="1" />
		<xsl:param name="values" select="''" />
		<xsl:param name="path" select="''" />
		<xsl:for-each select="/setup/param[$depth]/val">
			<xsl:variable name="newvalue">
				<xsl:element name="x">
					<xsl:attribute name="{normalize-space(/setup/param[$depth]/@name)}">
        <xsl:value-of select="normalize-space(.)" />
      </xsl:attribute>
					<xsl:for-each select="dep">
						<xsl:attribute name="{normalize-space(@name)}">
        <xsl:value-of select="normalize-space(@val)" />
      </xsl:attribute>
					</xsl:for-each>
				</xsl:element>
			</xsl:variable>
			<xsl:variable name="newpath">
				<xsl:value-of select="normalize-space(/setup/param[$depth]/@name)" />
				<xsl:text>_</xsl:text>
				<xsl:value-of select="normalize-space(.)" />
				<xsl:text>/</xsl:text>
			</xsl:variable>
			<xsl:choose>
				<xsl:when test="$depth &lt; (count(/setup/param))">
					<xsl:call-template name="expand">
						<xsl:with-param name="values">
							<xsl:copy-of select="$values" />
							<xsl:copy-of select="$newvalue" />
						</xsl:with-param>
						<xsl:with-param name="depth" select="$depth+1" />
						<xsl:with-param name="path">
							<xsl:value-of select="$path" />
							<xsl:value-of select="$newpath" />
						</xsl:with-param>
					</xsl:call-template>
					<!-- write subsetup.xml files on each node of the tree but the innermost -->
					<xsl:document href="{concat($path,$newpath)}/subsetup.xml"
						method="xml" encoding="UTF-8" indent="yes">
						<xsl:element name="setup">
							<xsl:value-of select="$newline" />
							<xsl:for-each select="/setup/param[position()&gt;$depth]">
								<xsl:copy-of select="." />
							</xsl:for-each>
							<xsl:value-of select="$newline" />
						</xsl:element>
					</xsl:document>
				</xsl:when>
				<xsl:otherwise>
					<xsl:element name="set">
						<xsl:for-each
							select="exsl:node-set($values)/*/@*|exsl:node-set($newvalue)/*/@*">
							<xsl:attribute name="{name(.)}">
                       <xsl:value-of select="." />
                       </xsl:attribute>
						</xsl:for-each>
						<xsl:attribute name="path">
                     <xsl:value-of select="$path" />
                     <xsl:value-of select="$newpath" />
                  </xsl:attribute>
					</xsl:element>
					<xsl:value-of select="$newline" />
				</xsl:otherwise>
			</xsl:choose>
		</xsl:for-each>
	</xsl:template>
	<!-- //////////////////////////////////////////////////////////////////////////////// -->
	<xsl:template match="/">
		<!-- main experiment.xml -->
		<xsl:document href="experiment.xml" method="xml" indent="yes">
			<xsl:element name="experiment">
				<xsl:value-of select="$newline" />
				<xsl:call-template name="expand" />
			</xsl:element>
		</xsl:document>
		<!-- subexperiment.xml files on every node of the tree but the innermost -->
		<xsl:call-template name="crawltree" />
	</xsl:template>
</xsl:stylesheet>
