<?xml version="1.0" encoding="UTF-8" ?>
<xsl:stylesheet version="1.0"
	xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
        xmlns:math="http://exslt.org/math">
<xsl:output method="text" encoding="UTF-8"/>

<xsl:template match="/">
   <xsl:if test="/thermodynamicproperties/vibrationalenergy">
      <xsl:for-each select="/thermodynamicproperties/vibrationalenergy">
         <xsl:document href="vibrationalenergy.agr" method="text">
            <xsl:variable name="C"><xsl:value-of select="math:power(10,3)"/></xsl:variable>
            <xsl:call-template name="insert">
               <xsl:with-param name="C"><xsl:value-of select="$C"/></xsl:with-param>
               <xsl:with-param name="y_unit"><xsl:value-of select="' [mHa]'"/></xsl:with-param>
               <xsl:with-param name="x_unit"><xsl:value-of select="' [K]'"/></xsl:with-param>
            </xsl:call-template>
         </xsl:document>
      </xsl:for-each>
   </xsl:if>

   <xsl:if test="/thermodynamicproperties/vibrationalfreeenergy">
      <xsl:for-each select="/thermodynamicproperties/vibrationalfreeenergy">
         <xsl:document href="vibrationalfreeenergy.agr" method="text">
           <xsl:variable name="C"><xsl:value-of select="math:power(10,3)"/></xsl:variable>
           <xsl:call-template name="insert">
              <xsl:with-param name="C"><xsl:value-of select="$C"/></xsl:with-param>
              <xsl:with-param name="y_unit"><xsl:value-of select="' [mHa]'"/></xsl:with-param>
              <xsl:with-param name="x_unit"><xsl:value-of select="' [K]'"/></xsl:with-param>
           </xsl:call-template>
         </xsl:document>
      </xsl:for-each>
   </xsl:if>

   <xsl:if test="/thermodynamicproperties/vibrationalentropy">
      <xsl:for-each select="/thermodynamicproperties/vibrationalentropy">
         <xsl:document href="vibrationalentropy.agr" method="text">
            <xsl:variable name="C"><xsl:value-of select="math:power(10,6)"/></xsl:variable>
            <xsl:call-template name="insert">
               <xsl:with-param name="C"><xsl:value-of select="$C"/></xsl:with-param>
               <xsl:with-param name="y_unit"><xsl:value-of select="' [\f{Symbol}m\f{}Ha/K]'"/></xsl:with-param>
               <xsl:with-param name="x_unit"><xsl:value-of select="' [K]'"/></xsl:with-param>
            </xsl:call-template>
         </xsl:document>
      </xsl:for-each>
   </xsl:if>

   <xsl:if test="/thermodynamicproperties/heatcapacity">
      <xsl:for-each select="/thermodynamicproperties/heatcapacity">
         <xsl:document href="heatcapacity.agr" method="text">
            <xsl:variable name="C"><xsl:value-of select="math:power(10,6)"/></xsl:variable>
            <xsl:call-template name="insert">
               <xsl:with-param name="C"><xsl:value-of select="$C"/></xsl:with-param>
               <xsl:with-param name="y_unit"><xsl:value-of select="' [\f{Symbol}m\f{}Ha/K]'"/></xsl:with-param>
               <xsl:with-param name="x_unit"><xsl:value-of select="' [K]'"/></xsl:with-param>
            </xsl:call-template>
         </xsl:document>
      </xsl:for-each>
   </xsl:if>
</xsl:template>

<xsl:template name="insert">
   <xsl:param name="C"/>
   <xsl:param name="y_unit"/>
   <xsl:param name="x_unit"/>

   <xsl:variable name="lowercase" select="'abcdefghijklmnopqrstuvwxyz'"/>
   <xsl:variable name="uppercase" select="'ABCDEFGHIJKLMNOPQRSTUVWXYZ'"/>

   <xsl:variable name="min_x_data">
      <xsl:for-each select="./map">
         <xsl:sort select="@variable1" data-type="number" order="ascending"/>
         <xsl:if test="position()=1">
            <xsl:value-of select="@variable1"/>
         </xsl:if>
      </xsl:for-each>
   </xsl:variable>

   <xsl:variable name="max_x_data">
      <xsl:for-each select="./map">
         <xsl:sort select="@variable1" data-type="number" order="descending"/>
         <xsl:if test="position()=1">
            <xsl:value-of select="@variable1"/>
         </xsl:if>
      </xsl:for-each>
   </xsl:variable>

   <xsl:variable name="min_y_data">
      <xsl:for-each select="./map">
         <xsl:sort select="@function1" data-type="number" order="ascending"/>
         <xsl:if test="position()=1">
            <xsl:value-of select="@function1 * $C"/>
         </xsl:if>
      </xsl:for-each>
    </xsl:variable>

    <xsl:variable name="max_y_data">
       <xsl:for-each select="./map">
          <xsl:sort select="@function1" data-type="number" order="descending"/>
          <xsl:if test="position()=1">
             <xsl:value-of select="@function1 * $C"/>
          </xsl:if>
       </xsl:for-each>
    </xsl:variable>

    <xsl:variable name="min_yaxis" select="$min_y_data - ($max_y_data - $min_y_data)*0.05"/>
    <xsl:variable name="max_yaxis" select="$max_y_data + ($max_y_data - $min_y_data)*0.05"/>

<xsl:text>
# Grace project file
#
@version 50122
@page size 792, 612
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
@reference date 0
@date wrap off
@date wrap year 1950
@default linewidth 2.5
@default linestyle 1
@default color 1
@default pattern 1
@default font 0
@default char size 1.000000
@default symbol size 1.000000
@default sformat "%.15g"
@background color 0
@page background fill on
@g0 on
@g0 hidden false
@g0 type XY
@g0 stacked false
@g0 bar hgap 0.000000
@with g0
@    world </xsl:text><xsl:value-of select="$min_x_data"/><xsl:text>, </xsl:text><xsl:value-of select="$min_yaxis"/><xsl:text>, </xsl:text><xsl:value-of select="$max_x_data"/><xsl:text>, </xsl:text><xsl:value-of select="$max_yaxis"/><xsl:text>
@    stack world 0, 0, 0, 0
@    znorm 1
@    view 0.230000, 0.150000, 1.200000, 0.850000
@    title ""
@    title font 4
@    title size 1.800000
@    title color 1
@    subtitle ""
@    subtitle font 4
@    subtitle size 1.200000
@    subtitle color 1
@    autoticks
@    xaxes scale Normal
@    yaxes scale Normal
@    xaxes invert off
@    yaxes invert off
@    xaxis  on
@    xaxis  type zero false
@    xaxis  offset 0.000000 , 0.000000
@    xaxis  bar on
@    xaxis  bar color 1
@    xaxis  bar linestyle 1
@    xaxis  bar linewidth 2.5
@    xaxis  label "</xsl:text><xsl:value-of select="concat(translate(substring(./mapdef/variable1/@name,1,1),$lowercase,$uppercase),substring(./mapdef/variable1/@name,2,string-length(./mapdef/variable1/@name)))"/> <xsl:value-of select="$x_unit"/><xsl:text>"
@    xaxis  label layout para
@    xaxis  label place auto
@    xaxis  label char size 1.800000
@    xaxis  label font 4
@    xaxis  label color 1
@    xaxis  label place normal
@    xaxis  tick on
@    xaxis  tick default 6
@    xaxis  tick place rounded true
@    xaxis  tick in
@    xaxis  tick major size 0.650000
@    xaxis  tick major color 1
@    xaxis  tick major linewidth 2.5
@    xaxis  tick major linestyle 1
@    xaxis  tick major grid off
@    xaxis  tick minor color 1
@    xaxis  tick minor linewidth 2.5
@    xaxis  tick minor linestyle 1
@    xaxis  tick minor grid off
@    xaxis  tick minor size 0.400000
@    xaxis  ticklabel on
@    xaxis  ticklabel formula ""
@    xaxis  ticklabel append ""
@    xaxis  ticklabel prepend ""
@    xaxis  ticklabel angle 0
@    xaxis  ticklabel skip 0
@    xaxis  ticklabel stagger 0
@    xaxis  ticklabel place normal
@    xaxis  ticklabel offset auto
@    xaxis  ticklabel offset 0.000000 , 0.010000
@    xaxis  ticklabel start type auto
@    xaxis  ticklabel start 0.000000
@    xaxis  ticklabel stop type auto
@    xaxis  ticklabel stop 0.000000
@    xaxis  ticklabel char size 1.600000
@    xaxis  ticklabel font 4
@    xaxis  ticklabel color 1
@    xaxis  tick place both
@    xaxis  tick spec type none
@    yaxis  on
@    yaxis  type zero false
@    yaxis  offset 0.000000 , 0.000000
@    yaxis  bar on
@    yaxis  bar color 1
@    yaxis  bar linestyle 1
@    yaxis  bar linewidth 3.0
@    yaxis  label "</xsl:text><xsl:value-of select="concat(translate(substring(./mapdef/function1/@name,1,1),$lowercase,$uppercase),substring(./mapdef/function1/@name,2,string-length(./mapdef/function1/@name)))"/> <xsl:value-of select="$y_unit"/><xsl:text>"
@    yaxis  label layout para
@    yaxis  label place auto
@    yaxis  label char size 1.800000
@    yaxis  label font 4
@    yaxis  label color 1
@    yaxis  label place normal
@    yaxis  tick on
@    yaxis  tick default 6
@    yaxis  tick place rounded true
@    yaxis  tick in
@    yaxis  tick major size 0.65000000
@    yaxis  tick major color 1
@    yaxis  tick major linewidth 2.5
@    yaxis  tick major linestyle 1
@    yaxis  tick major grid off
@    yaxis  tick minor color 1
@    yaxis  tick minor linewidth 2.5
@    yaxis  tick minor linestyle 1
@    yaxis  tick minor grid off
@    yaxis  tick minor size 0.400000
@    yaxis  ticklabel on
@    yaxis  ticklabel formula ""
@    yaxis  ticklabel append ""
@    yaxis  ticklabel prepend ""
@    yaxis  ticklabel angle 0
@    yaxis  ticklabel skip 0
@    yaxis  ticklabel stagger 0
@    yaxis  ticklabel place normal
@    yaxis  ticklabel offset auto
@    yaxis  ticklabel offset 0.000000 , 0.010000
@    yaxis  ticklabel start type auto
@    yaxis  ticklabel start 0.000000
@    yaxis  ticklabel stop type auto
@    yaxis  ticklabel stop 0.000000
@    yaxis  ticklabel char size 1.600000
@    yaxis  ticklabel font 4
@    yaxis  ticklabel color 1
@    yaxis  tick place both
@    yaxis  tick spec type none
@    s0 line color 2
</xsl:text>
   <xsl:if test="count(./map) &lt;= 20">
<xsl:text>
@    s0 symbol 3
@    s0 symbol size 1.000000
@    s0 symbol color 1
@    s0 symbol pattern 1
@    s0 symbol fill color 1
@    s0 symbol fill pattern 0
@    s0 symbol linewidth 2.5
@    s0 symbol linestyle 1
@    s0 symbol char 65
@    s0 symbol char font 0
@    s0 symbol skip 0
</xsl:text>
   </xsl:if>

<xsl:text>
@target G0.S0
@type xy
</xsl:text>
   <xsl:for-each select ="./map">
      <xsl:value-of select="@variable1"/><xsl:text> </xsl:text>
      <xsl:value-of select="@function1 * $C"/><xsl:text>
</xsl:text>
   </xsl:for-each>
<xsl:text>&amp;</xsl:text>
</xsl:template>

</xsl:stylesheet>

