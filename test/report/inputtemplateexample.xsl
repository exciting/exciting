<?xml version="1.0" encoding="UTF-8" ?>
<xsl:stylesheet version="1.0"
	xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
	<xsl:template match="/">

		<inputset xsi:noNamespaceSchemaLocation="./excitinginput.xsd"
			xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
			<xsl:for-each select="/parameterset/scale">
				<xsl:variable name="scale" select="." />
				<xsl:for-each select="/parameterset/rgkmax">
					<xsl:variable name="rgkmax" select="." />

					<input xsi:noNamespaceSchemaLocation="./excitinginput.xsd"
						xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
						<xsl:attribute name="id">
				<xsl:value-of select="concat('sc',$scale,'_rgk',$rgkmax)" />
				</xsl:attribute>
					<xsl:attribute name="depends">
			<xsl:value-of select="concat('sc',$scale,'_rgk',$rgkmax)" />
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
							<groundstate fromscratch="true" maxscl="30" solver="Arpack" mixer="msec" tforce="true" nempty="3" epsocc="1e-8">
								<xsl:attribute name="rgkmax"><xsl:value-of select="$rgkmax" /></xsl:attribute>

								<ngkgrit autokpt="true">0 1 3</ngkgrit>
								<spin spinpol="true">
									<mommtfix>
										<atommoment speciesnr="1" atomnr="1">2 4.23 3</atommoment>
									</mommtfix>
								</spin>
							</groundstate><bandstructure>
								<plot1d>
									<path steps="200">
										<point>0 0 0</point>
										<point>0 1 0</point>
										<point>8 6 5</point>
									</path>
								</plot1d>
							</bandstructure>
							<wfplot>
								<plot1d>
									<path steps="200">
										<point>0 0 0</point>
										<point>0 1 0</point>
									</path>
								</plot1d>
								<plot2d>
									<parallelogram grid="40 40">
										<origin>2 2 1e12</origin>
										<point>0 0 0</point>
										<point>3 4 3</point>
									</parallelogram>
								</plot2d>
							</wfplot>
						</tasks>
					</input>
				</xsl:for-each>
				<xsl:variable name="level1" select="position()" />
			</xsl:for-each>
		</inputset>
	</xsl:template>
</xsl:stylesheet>