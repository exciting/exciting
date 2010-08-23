<?xml version="1.0" encoding="UTF-8" ?>
<xsl:stylesheet version="1.0"
	xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:xs="http://www.w3.org/2001/XMLSchema">
	<xsl:output method="xml" encoding="UTF-8" indent="yes" />
	<xsl:variable name="newline">
		<xsl:text>		
</xsl:text>
	</xsl:variable>
	<xsl:template match="/">

		<input>
         <xsl:value-of select="$newline" />
			<xsl:comment>
				Use the "species.input.APW+lo-Legacy" input file, located in the
				"species" directory, for
				generating the species files related to the old default settings.
			</xsl:comment>
			<xsl:value-of select="$newline" />
			<xsl:for-each
				select="//xs:element[*/xs:attribute/xs:annotation/xs:appinfo/olddefault]">
				<xsl:element name="{@name}">
					<xsl:for-each select=".//xs:attribute[xs:annotation/xs:appinfo/olddefault]">
					   <xsl:if test="@name">
						   <xsl:attribute name="{@name}">
                        <xsl:value-of select="./xs:annotation/xs:appinfo/olddefault/@default" />
                     </xsl:attribute>
                  </xsl:if>
                  <xsl:if test="@ref">
                     <xsl:attribute name="{@ref}">
                        <xsl:value-of select="./xs:annotation/xs:appinfo/olddefault/@default" />
                     </xsl:attribute>
                  </xsl:if>
					</xsl:for-each>
				</xsl:element>
            <xsl:value-of select="$newline" />
			</xsl:for-each>
		</input>

	</xsl:template>
</xsl:stylesheet>