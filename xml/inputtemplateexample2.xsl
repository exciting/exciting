<?xml version="1.0" encoding="UTF-8" ?>
<xsl:stylesheet version="1.0"
	xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
	<xsl:template match="/">
		<inputset xsi:noNamespaceSchemaLocation="./excitinginput.xsd"
			xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
			<xsl:for-each select="/parameterset/scale">
				<xsl:variable name="scale" select="." />
				
				<xsl:for-each select="/parameterset/steps">
					<xsl:variable name="steps" select="." />
						<xsl:variable name="dogndstate" select="@dogndstate" />
				
					<input xsi:noNamespaceSchemaLocation="./excitinginput.xsd"
						xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
						<xsl:attribute name="id">
				<xsl:value-of select="concat('sc',$scale,'_steps',$steps)" />
				</xsl:attribute>
				
						   <xsl:attribute name="depends">
		<xsl:value-of select="concat('sc',$scale,'_steps',../steps[@dogndstate='true'])" />
				</xsl:attribute>

						<structure speciespath="../species">
							<crystal primcell="false">
								<xsl:attribute name="scale">
				<xsl:value-of select="$scale" />
				</xsl:attribute>
								<basevect>1 0 0</basevect>
								<basevect>0 1 0</basevect>
								<basevect>0 0 1</basevect>
							</crystal>
							<species abrev="Cu">
								<position>0 0 0</position>

							</species>
							<species abrev="O">
								<position>0 0 .5</position>
								<position>0 0 0</position>
							</species>
						</structure>

						<tasks>
							<xsl:if test="$dogndstate">
								<groundstate fromscratch="true" maxscl="30" solver="Arpack"
									mixer="msec" tforce="true" nempty="3" epsocc="1e-8">
									

									<ngkgrit autokpt="true">0 1 3</ngkgrit>
									<spin spinpol="true">
										<mommtfix>
											<atommoment speciesnr="1" atomnr="1">2 4.23 3</atommoment>
										</mommtfix>
									</spin>
								</groundstate>
							</xsl:if>
							<bandstructure>
								<plot1d>
									<path steps="">
										<xsl:attribute name="steps"><xsl:value-of
											select="$steps"></xsl:value-of></xsl:attribute>
										<point>0 0 0</point>
										<point>0 1 0</point>
										<point>8 6 5</point>
									</path>
								</plot1d>
							</bandstructure>
						</tasks>
					</input>
				</xsl:for-each>
				<xsl:variable name="level1" select="position()" />
			</xsl:for-each>
		</inputset>
	</xsl:template>
</xsl:stylesheet>