<?xml version="1.0" encoding="UTF-8" ?>
<xsl:stylesheet version="1.0"
	xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:xs="http://www.w3.org/2001/XMLSchema">
	<xsl:output method="xml" encoding="UTF-8" indent="yes" />
	<xsl:variable name="newline">
		<xsl:text />
	</xsl:variable>
	<xsl:template match="/">

		<input>
			<xsl:value-of select="$newline" />
			<xsl:for-each
				select="//xs:element[*/xs:attribute/xs:annotation/xs:appinfo/olddefault]">
				<xsl:element name="{@name}">
					<xsl:for-each
						select=".//xs:attribute[xs:annotation/xs:appinfo/olddefault]">
						<xsl:attribute name="{@name}">
              <xsl:value-of select="./xs:annotation/xs:appinfo/olddefault/@default" />
              </xsl:attribute>
					</xsl:for-each>
				</xsl:element>
				<xsl:value-of select="$newline" />
			</xsl:for-each>
		</input>

	</xsl:template>
</xsl:stylesheet>