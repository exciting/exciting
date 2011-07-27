<?xml version="1.0" encoding="UTF-8" ?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
   <xsl:output method="text" />
<xsl:param name="unit"/>

<xsl:template match="/">
   <xsl:variable name="C">
      <xsl:choose>
         <xsl:when test="$unit='Ha' or $unit='Hartree'">
            <xsl:value-of select="1"/>
         </xsl:when>
         <xsl:otherwise>
            <xsl:value-of select="27.211384"/>
         </xsl:otherwise>
      </xsl:choose>
   </xsl:variable>

   <xsl:if test="/bandstructure/species/atom/band">
      <xsl:for-each select="/bandstructure/species[1]/atom[1]/band[1]/point[1]/bc">
         <xsl:variable name="l" select="@l"/>
         <xsl:for-each select="/bandstructure/species">
            <xsl:for-each select="./atom">
               <xsl:variable name="filename">
                  <xsl:value-of select="/bandstructure/title"/>
                  <xsl:text>_band_</xsl:text><xsl:value-of select="../@name"/>
                  <xsl:text>_</xsl:text><xsl:value-of select="position()"/>
                  <xsl:text>_l</xsl:text><xsl:value-of select="$l"/><xsl:text>.agr</xsl:text>
               </xsl:variable>
               <xsl:document href="{$filename}"  method="text">
                  <xsl:call-template name="bandstograce">
                     <xsl:with-param name="data" select="./band"/>
                     <xsl:with-param name="C" select="$C"/>
                     <xsl:with-param name="l" select="$l"/>
                  </xsl:call-template>
               </xsl:document>
            </xsl:for-each>
         </xsl:for-each>
      </xsl:for-each>
   </xsl:if>

   <xsl:if test="/bandstructure/band">
      <xsl:variable name="filename">
         <xsl:value-of select="/bandstructure/title"/><xsl:text>_bandstructure.agr</xsl:text>
      </xsl:variable>
      <xsl:document href="{$filename}" method="text">
         <xsl:call-template name="bandstograce">
            <xsl:with-param name="data" select="bandstructure/band"/>
            <xsl:with-param name="C" select="$C"/>
         </xsl:call-template>
      </xsl:document>
   </xsl:if>
</xsl:template>

<xsl:template name="bandstograce">
   <xsl:param name="data"/>
   <xsl:param name="C"/>
   <xsl:param name="l" select="false"/>
 
   <xsl:variable name="min_x_data">
      <xsl:for-each select="$data/point">
         <xsl:sort select="@distance" data-type="number" order="ascending"/>
         <xsl:if test="position()=1">
            <xsl:value-of select="@distance"/>
         </xsl:if>
      </xsl:for-each>
   </xsl:variable>

   <xsl:variable name="max_x_data">
      <xsl:for-each select="$data/point">
         <xsl:sort select="@distance" data-type="number" order="descending"/>
         <xsl:if test="position()=1">
            <xsl:value-of select="@distance"/>
         </xsl:if>
      </xsl:for-each>
   </xsl:variable>

   <xsl:variable name="min_y_data">
      <xsl:for-each select="$data/point">
         <xsl:sort select="@eval" data-type="number" order="ascending"/>
         <xsl:if test="position()=1">
            <xsl:value-of select="@eval*$C"/>
         </xsl:if>
      </xsl:for-each>
    </xsl:variable>

    <xsl:variable name="max_y_data">
       <xsl:for-each select="$data/point">
          <xsl:sort select="@eval" data-type="number" order="descending"/>
          <xsl:if test="position()=1">
             <xsl:value-of select="@eval*$C"/>
          </xsl:if>
       </xsl:for-each>
    </xsl:variable>

    <xsl:variable name="min_yaxis" select="$min_y_data - ($max_y_data - $min_y_data)*0.05"/>
    <xsl:variable name="max_yaxis" select="$max_y_data + ($max_y_data - $min_y_data)*0.05"/>

<xsl:text>
# Grace project file
#
@version 50122
@page size 600,600 
@page scroll 5%
@page inout 5%
@link page off
@map font 8 to "Courier", "Courier"
@map font 10 to "Courier-Bold", "Courier-Bold"
@map font 11 to "Courier-BoldOblique", "Courier-BoldOblique"
@map font 9 to "Courier-Oblique", "Courier-Oblique"
@map font 4 to "Helvetica", "Helvetica"
@map font 6 to "Helvetica-Bold", "Helvetica-Bold"
@map font 7 to "Helvetica-BoldOblique", "Helvetica-BoldOblique"
@map font 15 to "Helvetica-Narrow", "Helvetica-Narrow"
@map font 16 to "Helvetica-Narrow-Bold", "Helvetica-Narrow-Bold"
@map font 17 to "Helvetica-Narrow-BoldOblique", "Helvetica-Narrow-BoldOblique"
@map font 18 to "Helvetica-Narrow-Oblique", "Helvetica-Narrow-Oblique"
@map font 5 to "Helvetica-Oblique", "Helvetica-Oblique"
@map font 20 to "NewCenturySchlbk-Bold", "NewCenturySchlbk-Bold"
@map font 21 to "NewCenturySchlbk-BoldItalic", "NewCenturySchlbk-BoldItalic"
@map font 22 to "NewCenturySchlbk-Italic", "NewCenturySchlbk-Italic"
@map font 23 to "NewCenturySchlbk-Roman", "NewCenturySchlbk-Roman"
@map font 24 to "Palatino-Bold", "Palatino-Bold"
@map font 25 to "Palatino-BoldItalic", "Palatino-BoldItalic"
@map font 26 to "Palatino-Italic", "Palatino-Italic"
@map font 27 to "Palatino-Roman", "Palatino-Roman"
@map font 12 to "Symbol", "Symbol"
@map font 2 to "Times-Bold", "Times-Bold"
@map font 3 to "Times-BoldItalic", "Times-BoldItalic"
@map font 1 to "Times-Italic", "Times-Italic"
@map font 0 to "Times-Roman", "Times-Roman"
@map font 33 to "ZapfChancery-MediumItalic", "ZapfChancery-MediumItalic"
@map font 13 to "ZapfDingbats", "ZapfDingbats"
@map font 35 to "LMRoman10-Bold", "LMRoman10-Bold"
@map font 36 to "LMRoman10-BoldItalic", "LMRoman10-BoldItalic"
@map font 37 to "LMRoman10-BoldOblique", "LMRoman10-BoldOblique"
@map font 38 to "LMRoman10-CapsOblique", "LMRoman10-CapsOblique"
@map font 39 to "LMRoman10-CapsRegular", "LMRoman10-CapsRegular"
@map font 40 to "LMRoman10-Demi", "LMRoman10-Demi"
@map font 41 to "LMRoman10-DemiOblique", "LMRoman10-DemiOblique"
@map font 42 to "LMRoman10-Dunhill", "LMRoman10-Dunhill"
@map font 43 to "LMRoman10-DunhillOblique", "LMRoman10-DunhillOblique"
@map font 44 to "LMRoman10-Italic", "LMRoman10-Italic"
@map font 45 to "LMRoman10-Oblique", "LMRoman10-Oblique"
@map font 46 to "LMRoman10-Regular", "LMRoman10-Regular"
@map font 47 to "LMRoman10-Unslanted", "LMRoman10-Unslanted"
@map font 48 to "LMSans10-Bold", "LMSans10-Bold"
@map font 49 to "LMSans10-BoldOblique", "LMSans10-BoldOblique"
@map font 50 to "LMSans10-DemiCondensed", "LMSans10-DemiCondensed"
@map font 51 to "LMSans10-DemiCondensedOblique", "LMSans10-DemiCondensedOblique"
@map font 52 to "LMSans10-Oblique", "LMSans10-Oblique"
@map font 53 to "LMSans10-Regular", "LMSans10-Regular"
@map font 54 to "LMSansQuotation8-Bold", "LMSansQuotation8-Bold"
@map font 55 to "LMSansQuotation8-BoldOblique", "LMSansQuotation8-BoldOblique"
@map font 56 to "LMSansQuotation8-Oblique", "LMSansQuotation8-Oblique"
@map font 57 to "LMSansQuotation8-Regular", "LMSansQuotation8-Regular"
@map font 58 to "LMTypewriter10-CapsOblique", "LMTypewriter10-CapsOblique"
@map font 59 to "LMTypewriter10-CapsRegular", "LMTypewriter10-CapsRegular"
@map font 60 to "LMTypewriter10-Dark", "LMTypewriter10-Dark"
@map font 61 to "LMTypewriter10-DarkOblique", "LMTypewriter10-DarkOblique"
@map font 62 to "LMTypewriter10-Italic", "LMTypewriter10-Italic"
@map font 63 to "LMTypewriter10-Light", "LMTypewriter10-Light"
@map font 64 to "LMTypewriter10-LightCondensed", "LMTypewriter10-LightCondensed"
@map font 65 to "LMTypewriter10-LightCondensedOblique", "LMTypewriter10-LightCondensedOblique"
@map font 66 to "LMTypewriter10-LightOblique", "LMTypewriter10-LightOblique"
@map font 67 to "LMTypewriter10-Oblique", "LMTypewriter10-Oblique"
@map font 68 to "LMTypewriter10-Regular", "LMTypewriter10-Regular"
@map font 69 to "LMTypewriterVarWd10-Dark", "LMTypewriterVarWd10-Dark"
@map font 70 to "LMTypewriterVarWd10-DarkOblique", "LMTypewriterVarWd10-DarkOblique"
@map font 71 to "LMTypewriterVarWd10-Light", "LMTypewriterVarWd10-Light"
@map font 72 to "LMTypewriterVarWd10-LightOblique", "LMTypewriterVarWd10-LightOblique"
@map font 73 to "LMTypewriterVarWd10-Oblique", "LMTypewriterVarWd10-Oblique"
@map font 74 to "LMTypewriterVarWd10-Regular", "LMTypewriterVarWd10-Regular"
@map font 75 to "TeXGyreAdventor-Bold", "TeXGyreAdventor-Bold"
@map font 76 to "TeXGyreAdventor-BoldItalic", "TeXGyreAdventor-BoldItalic"
@map font 77 to "TeXGyreAdventor-Italic", "TeXGyreAdventor-Italic"
@map font 78 to "TeXGyreAdventor-Regular", "TeXGyreAdventor-Regular"
@map font 79 to "TeXGyreBonum-Bold", "TeXGyreBonum-Bold"
@map font 80 to "TeXGyreBonum-BoldItalic", "TeXGyreBonum-BoldItalic"
@map font 81 to "TeXGyreBonum-Italic", "TeXGyreBonum-Italic"
@map font 82 to "TeXGyreBonum-Regular", "TeXGyreBonum-Regular"
@map font 83 to "TeXGyreChorus-Italic", "TeXGyreChorus-Italic"
@map font 84 to "TeXGyreCursor-Bold", "TeXGyreCursor-Bold"
@map font 85 to "TeXGyreCursor-BoldItalic", "TeXGyreCursor-BoldItalic"
@map font 86 to "TeXGyreCursor-Italic", "TeXGyreCursor-Italic"
@map font 87 to "TeXGyreCursor-Regular", "TeXGyreCursor-Regular"
@map font 88 to "TeXGyreHeros-Bold", "TeXGyreHeros-Bold"
@map font 89 to "TeXGyreHeros-BoldItalic", "TeXGyreHeros-BoldItalic"
@map font 90 to "TeXGyreHeros-Condensed", "TeXGyreHeros-Condensed"
@map font 91 to "TeXGyreHeros-CondensedBold", "TeXGyreHeros-CondensedBold"
@map font 92 to "TeXGyreHeros-CondensedBoldItalic", "TeXGyreHeros-CondensedBoldItalic"
@map font 93 to "TeXGyreHeros-CondensedItalic", "TeXGyreHeros-CondensedItalic"
@map font 94 to "TeXGyreHeros-Italic", "TeXGyreHeros-Italic"
@map font 95 to "TeXGyreHeros-Regular", "TeXGyreHeros-Regular"
@map font 96 to "TeXGyrePagella-Bold", "TeXGyrePagella-Bold"
@map font 97 to "TeXGyrePagella-BoldItalic", "TeXGyrePagella-BoldItalic"
@map font 98 to "TeXGyrePagella-Italic", "TeXGyrePagella-Italic"
@map font 99 to "TeXGyrePagella-Regular", "TeXGyrePagella-Regular"
@map font 100 to "TeXGyreSchola-Bold", "TeXGyreSchola-Bold"
@map font 101 to "TeXGyreSchola-BoldItalic", "TeXGyreSchola-BoldItalic"
@map font 102 to "TeXGyreSchola-Italic", "TeXGyreSchola-Italic"
@map font 103 to "TeXGyreSchola-Regular", "TeXGyreSchola-Regular"
@map font 104 to "TeXGyreTermes-Bold", "TeXGyreTermes-Bold"
@map font 105 to "TeXGyreTermes-BoldItalic", "TeXGyreTermes-BoldItalic"
@map font 106 to "TeXGyreTermes-Italic", "TeXGyreTermes-Italic"
@map color 0 to (255, 255, 255), "white"
@map color 1 to (0, 0, 0), "black"
@map color 2 to (255, 0, 0), "red"
@map color 3 to (0, 0, 255), "blue"
@map color 4 to (0, 255, 0), "green"
@map color 5 to (255, 255, 0), "yellow"
@map color 6 to (188, 143, 143), "brown"
@map color 7 to (100, 100, 100), "grey"
@map color 8 to (148, 0, 211), "violet"
@map color 9 to (0, 255, 255), "cyan"
@map color 10 to (255, 0, 255), "magenta"
@map color 11 to (255, 165, 0), "orange"
@map color 12 to (114, 33, 188), "indigo"
@map color 13 to (103, 7, 72), "maroon"
@map color 14 to (64, 224, 208), "turquoise"
@map color 15 to (0, 139, 0), "green4"
@default linewidth 2.5 
@default linestyle 1
@default pattern 1
@default font 0
@default char size 1.000000
@default symbol size 1.000000
@default sformat "%.15g"
@background color 0
@page background fill on
@ with line
@     line on
@     line loctype world
@     line g0
@     line 0, 0,</xsl:text> <xsl:value-of select="//point[last()]/@distance"></xsl:value-of> <xsl:text>, 0
@     line linewidth 2.5
@     line linestyle 3
@ line def
@ with string
@     string on
@     string loctype world
@     string g0
@     string </xsl:text> <xsl:value-of select="//point[last()]/@distance*1.01 "></xsl:value-of> <xsl:text>, </xsl:text><xsl:value-of select="($min_y_data - $max_y_data)*0.0070"/><xsl:text>
@     string font 4
@     string char size 1.650000
@     string def "E\sF"
@ r0 off
@ with g0
@     view 0.160000, 0.120000, 0.850000, 0.85
@     title "</xsl:text><xsl:value-of select="/bandstructure/title"/><xsl:text>"
@     title font 4
@     title size 1.800000
@     title color 1
@     subtitle size 1.100000
</xsl:text><xsl:if test="$l"><xsl:text>
@     subtitle "</xsl:text><xsl:value-of select="../@chemicalSymbol"/><xsl:text>, l= </xsl:text><xsl:value-of select="$l"/><xsl:text>"</xsl:text>
</xsl:if><xsl:text>
@     world </xsl:text><xsl:value-of select="$min_x_data"/><xsl:text>,</xsl:text><xsl:value-of select="$min_yaxis"/><xsl:text>,</xsl:text><xsl:value-of select="$max_x_data"/><xsl:text>,</xsl:text><xsl:value-of select="$max_yaxis"/><xsl:text>
@     yaxis  ticklabel font 4
@     yaxis  label font 4
@     yaxis  label "Energy [</xsl:text><xsl:if test="$C=1"><xsl:value-of select="'Ha'"/></xsl:if>
                            <xsl:if test="$C!=1"><xsl:value-of select="'eV'"/></xsl:if><xsl:text>]"
@     yaxis  label char size 1.800000
@     yaxis  ticklabel char size 1.50000
@     yaxis  bar linewidth 3.0
@     yaxis  tick place rounded true
@     yaxis  tick in
@     yaxis  tick major size 0.800000
@     yaxis  tick major color 1
@     yaxis  tick major linewidth 3.0
@     yaxis  tick major linestyle 1
@     yaxis  tick major grid off
@     yaxis  tick minor color 1
@     yaxis  tick minor linewidth 3.0
@     yaxis  tick minor linestyle 1
@     yaxis  tick minor grid off
@     yaxis  tick minor size 0.500000
@     xaxis  bar on
@     xaxis  bar color 1
@     xaxis  bar linestyle 1
@     xaxis  bar linewidth 3.0
@     xaxis  tick on
@     xaxis  tick major 0.5
@     xaxis  tick minor ticks 1
@     xaxis  tick default 6
@     xaxis  tick place rounded true
@     xaxis  tick in
@     xaxis  tick major size 1.000000
@     xaxis  tick major color 1
@     xaxis  tick major linewidth 2.0
@     xaxis  tick major linestyle 1
@     xaxis  tick major grid on
@     xaxis  tick minor color 1
@     xaxis  ticklabel font 0 
@     xaxis  label font 4
@     xaxis  label place normal
@     xaxis  ticklabel offset spec
@     xaxis  ticklabel offset 0.000000 , 0.016000
@     xaxis  ticklabel char size 1.65000
@     xaxis  tick place both
@     xaxis  tick spec type both
@     xaxis  tick major grid on
@     xaxis  tick spec   </xsl:text>
   <xsl:value-of select="count(/bandstructure/vertex)" />
<xsl:text>
</xsl:text>
   <xsl:for-each select="/bandstructure/vertex">
<xsl:text>@     xaxis  tick major  </xsl:text>
      <xsl:value-of select="position()-1" />
<xsl:text> , </xsl:text>
      <xsl:value-of select="@distance" />
<xsl:text>
</xsl:text>
<xsl:text>@     xaxis  ticklabel </xsl:text>
      <xsl:value-of select="position()-1" />
<xsl:text>, "</xsl:text>
      <xsl:choose>
         <xsl:when test="@label='GAMMA'">
<xsl:text>\xG\f{}</xsl:text>
         </xsl:when>
         <xsl:when test="@label='Gamma'">
<xsl:text>\xG\f{}</xsl:text>
         </xsl:when>
         <xsl:when test="@label='DELTA'">
<xsl:text>\xD\f{}</xsl:text>
         </xsl:when>
         <xsl:when test="@label='Delta'">
<xsl:text>\xD\f{}</xsl:text>
         </xsl:when>
         <xsl:when test="@label='THETA'">
<xsl:text>\xQ\f{}</xsl:text>
         </xsl:when>
         <xsl:when test="@label='Theta'">
<xsl:text>\xQ\f{}</xsl:text>
         </xsl:when>
         <xsl:when test="@label='XI'">
<xsl:text>\xX\f{}</xsl:text>
         </xsl:when>
         <xsl:when test="@label='Xi'">
<xsl:text>\xX\f{}</xsl:text>
         </xsl:when>
         <xsl:when test="@label='PI'">
<xsl:text>\xP\f{}</xsl:text>
         </xsl:when>
         <xsl:when test="@label='Pi'">
<xsl:text>\xP\f{}</xsl:text>
         </xsl:when>
         <xsl:when test="@label='SIGMA'">
<xsl:text>\xS\f{}</xsl:text>
         </xsl:when>
         <xsl:when test="@label='Sigma'">
<xsl:text>\xS\f{}</xsl:text>
         </xsl:when>
         <xsl:when test="@label='PHI'">
<xsl:text>\xF\f{}</xsl:text>
         </xsl:when>
         <xsl:when test="@label='Phi'">
<xsl:text>\xF\f{}</xsl:text>
         </xsl:when>
         <xsl:when test="@label='PSI'">
<xsl:text>\xY\f{}</xsl:text>
         </xsl:when>
         <xsl:when test="@label='Psi'">
<xsl:text>\xY\f{}</xsl:text>
         </xsl:when>
         <xsl:when test="@label='OMEGA'">
<xsl:text>\xW\f{}</xsl:text>
         </xsl:when>
         <xsl:when test="@label='Omega'">
<xsl:text>\xW\f{}</xsl:text>
         </xsl:when>
         <xsl:otherwise>
            <xsl:value-of select="@label" />
         </xsl:otherwise>
      </xsl:choose>
<xsl:text>"
</xsl:text>
   </xsl:for-each>
<xsl:text>@ autoticks
</xsl:text>

   <xsl:if test="not($l)">
      <xsl:for-each select="$data">
         <xsl:variable name="setnr" select="position()-1"/>
         <xsl:variable name="species" select="../../@chemicalSymbol"/>
         <xsl:if test="point/bc">
            <xsl:for-each select="point">
               <xsl:variable name="x" select="@distance"/>
               <xsl:variable name="y"> 
                  <xsl:value-of select="@eval*$C"/>
               </xsl:variable>
               <xsl:variable name="n" select="position()"/>
               <xsl:for-each select="bc">
                  <xsl:if test="@character>0.1 and $n mod 2=@l mod 2">
@with string
@    string on
@    string loctype world
@    string g0
@    string <xsl:value-of select="$x"/> , <xsl:value-of select="$y"/>
@    string color <xsl:value-of select="@l"/>
@    string rot 40
@    string font 4 
@    string just 0
@    string char size <xsl:value-of select="@character*2"/>
@    string def "<xsl:value-of select="$species"/> l=<xsl:value-of select="@l"/>"
                  </xsl:if>
               </xsl:for-each>
            </xsl:for-each>
         </xsl:if>
      </xsl:for-each>
   </xsl:if>

   <xsl:for-each select="$data">
<xsl:text>
@target G0.S</xsl:text>
      <xsl:value-of select="position()" />
      <xsl:choose>
         <xsl:when test="$l">
<xsl:text>
@    s</xsl:text><xsl:value-of select="position()"/><xsl:text> line color 2
@    s</xsl:text><xsl:value-of select="position()"/><xsl:text> symbol 11
@    s</xsl:text><xsl:value-of select="position()"/><xsl:text> symbol size 1.000000
@    s</xsl:text><xsl:value-of select="position()"/><xsl:text> symbol color 1
@    s</xsl:text><xsl:value-of select="position()"/><xsl:text> symbol pattern 1
@    s</xsl:text><xsl:value-of select="position()"/><xsl:text> symbol fill color 1
@    s</xsl:text><xsl:value-of select="position()"/><xsl:text> symbol fill pattern 0
@    s</xsl:text><xsl:value-of select="position()"/><xsl:text> symbol linewidth 1.5
@    s</xsl:text><xsl:value-of select="position()"/><xsl:text> symbol linestyle 1
@    s</xsl:text><xsl:value-of select="position()"/><xsl:text> symbol char 124
@    s</xsl:text><xsl:value-of select="position()"/><xsl:text> symbol char font 4
@    s</xsl:text><xsl:value-of select="position()"/><xsl:text> symbol skip 0
@type xysize
</xsl:text>
         </xsl:when>
         <xsl:otherwise>
<xsl:text>
@    s</xsl:text><xsl:value-of select="position()"/><xsl:text> line color 3
@type xy
</xsl:text>
         </xsl:otherwise>
      </xsl:choose>

      <xsl:for-each select="./point">
         <xsl:value-of select="@distance" />
<xsl:text>     </xsl:text>
         <xsl:value-of select="@eval*$C"/>

         <xsl:if test="$l">
<xsl:text>     </xsl:text>
            <xsl:value-of select="bc[@l=$l]/@character*2" />
         </xsl:if>
<xsl:text>
</xsl:text>
      </xsl:for-each>
<xsl:text>&amp;
</xsl:text>
   </xsl:for-each>
</xsl:template>

</xsl:stylesheet>
