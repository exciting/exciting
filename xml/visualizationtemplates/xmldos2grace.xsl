<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
                xmlns:exsl="http://exslt.org/common"
                extension-element-prefixes="exsl"
                version="1.0">
<xsl:import href="extract_diagram.xsl"/>
<xsl:output method="text" media-type="application/octet-stream"/>
<xsl:template match="/">

   <xsl:variable name="content">
      <xsl:apply-imports/>
   </xsl:variable>

   <xsl:variable name="min_x_data">
      <xsl:for-each select="exsl:node-set($content)/dos/*/diagram/point">
         <xsl:sort select="@e" data-type="number" order="ascending"/>
         <xsl:if test="position()=1">
            <xsl:value-of select="@e"/>
         </xsl:if>
      </xsl:for-each>
   </xsl:variable>

   <xsl:variable name="max_x_data">
      <xsl:for-each select="exsl:node-set($content)/dos/*/diagram/point">
         <xsl:sort select="@e" data-type="number" order="descending"/>
         <xsl:if test="position()=1">
            <xsl:value-of select="@e"/>
         </xsl:if>
      </xsl:for-each>
   </xsl:variable>

   <xsl:variable name="min_y_data">
      <xsl:for-each select="exsl:node-set($content)/dos/*/diagram/point">
         <xsl:sort select="@dos" data-type="number" order="ascending"/>
         <xsl:if test="position()=1">
            <xsl:value-of select="@dos"/>
         </xsl:if>
      </xsl:for-each>
    </xsl:variable>

    <xsl:variable name="max_y_data">
       <xsl:for-each select="exsl:node-set($content)/dos/*/diagram/point">
          <xsl:sort select="@dos" data-type="number" order="descending"/>
          <xsl:if test="position()=1">
             <xsl:value-of select="@dos"/>
          </xsl:if>
       </xsl:for-each>
    </xsl:variable>

    <xsl:variable name="min_xaxis" select="$min_x_data - ($max_x_data - $min_x_data)*0.05"/>
    <xsl:variable name="max_xaxis" select="$max_x_data + ($max_x_data - $min_x_data)*0.05"/>
    <xsl:variable name="min_yaxis" select="$min_y_data - ($max_y_data - $min_y_data)*0.05"/>
    <xsl:variable name="max_yaxis" select="$max_y_data + ($max_y_data - $min_y_data)*0.05"/>
<xsl:text>
# Grace project file
#
@version 50122
@page size 842, 595
@page scroll 5%
@page inout 5%
@link page off
@map font 0 to "Times-Roman", "Times-Roman"
@map font 1 to "Times-Italic", "Times-Italic"
@map font 2 to "Times-Bold", "Times-Bold"
@map font 3 to "Times-BoldItalic", "Times-BoldItalic"
@map font 4 to "Helvetica", "Helvetica"
@map font 5 to "Helvetica-Oblique", "Helvetica-Oblique"
@map font 6 to "Helvetica-Bold", "Helvetica-Bold"
@map font 7 to "Helvetica-BoldOblique", "Helvetica-BoldOblique"
@map font 8 to "Courier", "Courier"
@map font 9 to "Courier-Oblique", "Courier-Oblique"
@map font 10 to "Courier-Bold", "Courier-Bold"
@map font 11 to "Courier-BoldOblique", "Courier-BoldOblique"
@map font 13 to "ZapfDingbats", "ZapfDingbats"
@map font 15 to "Helvetica-Narrow", "Helvetica-Narrow"
@map font 16 to "Helvetica-Narrow-Bold", "Helvetica-Narrow-Bold"
@map font 17 to "Helvetica-Narrow-BoldOblique", "Helvetica-Narrow-BoldOblique"
@map font 18 to "Helvetica-Narrow-Oblique", "Helvetica-Narrow-Oblique"
@map font 20 to "NewCenturySchlbk-Bold", "NewCenturySchlbk-Bold"
@map font 21 to "NewCenturySchlbk-BoldItalic", "NewCenturySchlbk-BoldItalic"
@map font 22 to "NewCenturySchlbk-Italic", "NewCenturySchlbk-Italic"
@map font 23 to "NewCenturySchlbk-Roman", "NewCenturySchlbk-Roman"
@map font 24 to "Palatino-Bold", "Palatino-Bold"
@map font 25 to "Palatino-BoldItalic", "Palatino-BoldItalic"
@map font 26 to "Palatino-Italic", "Palatino-Italic"
@map font 27 to "Palatino-Roman", "Palatino-Roman"
@map color 0 to (255, 255, 255), "white"
@map color 1 to (0, 0, 0), "black"
@map color 2 to (255, 0, 0), "red"
@map color 3 to (0,255,0), "lime"
@map color 4 to (0, 0, 255), "blue"
@map color 5 to (0, 255, 255), "cyan"
@map color 6 to (140, 0, 140), "purple"
@map color 7 to (255, 140, 0), "darkorange"
@map color 8 to (255, 0, 255), "magneta"
@map color 9 to (25, 25, 112), "midnightblue"
@map color 10 to (0, 128, 0), "green"
@map color 11 to (140, 0, 0 ), "maroon"
@map color 12 to (255, 215, 0), "gold"
@map color 13 to (189, 183, 107), "darkkhaki"
@map color 14 to (199, 21, 133), "MediumVioletRed"
@map color 15 to (220, 220, 220), "grey"
@reference date 0
@date wrap off
@date wrap year 1950
@default linewidth 1.0
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
@    world </xsl:text><xsl:value-of select="$min_xaxis"/><xsl:text>,</xsl:text><xsl:value-of select="$min_yaxis"/><xsl:text>,</xsl:text><xsl:value-of select="$max_xaxis"/><xsl:text>,</xsl:text><xsl:value-of select="$max_yaxis"/><xsl:text>
@    stack world 0, 0, 0, 0
@    znorm 1
@    view 0.275000, 0.160000, 1.275000, 0.850000
@    title "</xsl:text><xsl:value-of select="/dos/title"/><xsl:text>\s "
@    title font 4
@    title size 1.500000
@    title color 1
@    subtitle "Density Of States"
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
@    xaxis  bar linewidth 3.0
@    xaxis  label "Energy (Ha)"
@    xaxis  label layout para
@    xaxis  label place spec
@    xaxis  label place 0.000000, 0.110000
@    xaxis  label char size 1.800000
@    xaxis  label font 4
@    xaxis  label color 1
@    xaxis  label place normal
@    xaxis  tick off
</xsl:text><!--
@    xaxis  tick major 0.05
--><xsl:text>
@    xaxis  tick minor ticks 4
@    xaxis  tick default 6
@    xaxis  tick place rounded true
@    xaxis  tick in
@    xaxis  tick major size 1.000000
@    xaxis  tick major color 15
@    xaxis  tick major linewidth 1.5
@    xaxis  tick major linestyle 1
@    xaxis  tick major grid on
@    xaxis  tick minor color 15 
@    xaxis  tick minor linewidth 1.0
@    xaxis  tick minor linestyle 1
@    xaxis  tick minor grid on
@    xaxis  tick minor size 0.500000
@    xaxis  ticklabel on
</xsl:text><!--
@    xaxis  ticklabel format decimal
@    xaxis  ticklabel prec 2
--><xsl:text>
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
@    yaxis  label "DOS (state/Ha)"
@    yaxis  label layout para
@    yaxis  label place spec
@    yaxis  label place 0.000000, 0.200000
@    yaxis  label char size 1.800000
@    yaxis  label font 4
@    yaxis  label color 1 
@    yaxis  label place normal 
@    yaxis  tick off
</xsl:text><!-- 
@    yaxis  tick major 500
--><xsl:text>  
@    yaxis  tick minor ticks 4
@    yaxis  tick default 6
@    yaxis  tick place rounded true
@    yaxis  tick in
@    yaxis  tick major size 1.000000
@    yaxis  tick major color 15 
@    yaxis  tick major linewidth 1.5
@    yaxis  tick major linestyle 1
@    yaxis  tick major grid on
@    yaxis  tick minor color 15 
@    yaxis  tick minor linewidth 1.0
@    yaxis  tick minor linestyle 1
@    yaxis  tick minor grid on
@    yaxis  tick minor size 0.500000
@    yaxis  ticklabel on
</xsl:text><!--
@    yaxis  ticklabel format decimal
@    yaxis  ticklabel prec 3
--><xsl:text>
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
@    legend on
@    legend loctype view
@    legend 0.30,0.82
@    legend box color 1
@    legend box pattern 1
@    legend box linewidth 2.0
@    legend box linestyle 1
@    legend box fill color 0
@    legend box fill pattern 1
@    legend font 4
@    legend char size 1.190000
@    legend color 1
@    legend length 4
@    legend vgap 1
@    legend hgap 1
@    legend invert false
@    frame type 0
@    frame linestyle 0
@    frame linewidth 1.0
@    frame color 1
@    frame pattern 1
@    frame background color 0
@    frame background pattern 0
</xsl:text>
<xsl:variable name="count_partial" select="count(exsl:node-set($content)/dos/partialdos/diagram)"/>
<xsl:variable name="count_total" select="count(exsl:node-set($content)/dos/totaldos/diagram)"/>
<xsl:variable name="count_interstitial" select="count(exsl:node-set($content)/dos/interstitialdos/diagram)"/>
<xsl:variable name="counter" select="$count_partial+$count_total+$count_interstitial"/>

<xsl:variable name="lll">
   <xsl:for-each select="exsl:node-set($content)/dos/*">
      <xsl:variable name="leg_p" select="concat(@speciessym,', ',@atom)"/>
      <xsl:for-each select="./diagram">
         <xsl:variable name="s0">
            <xsl:if test="@nspin=1"><xsl:value-of select="'up'"/></xsl:if>
            <xsl:if test="@nspin=2"><xsl:value-of select="'dn'"/></xsl:if>
         </xsl:variable>
         <xsl:variable name="legend">
            <xsl:choose>
            <xsl:when test="@type='totaldos'">            
               <xsl:value-of select="concat('t, s=',$s0)"/>
            </xsl:when>
            <xsl:when test="@type='interstitial'">
               <xsl:value-of select="concat('i, s=',$s0)"/>
            </xsl:when>
            <xsl:otherwise>
               <xsl:value-of select="concat($leg_p,', ', 's=',$s0,', ', 'l=',@l, ', ','m=',@m)"/>
            </xsl:otherwise>
            </xsl:choose>
         </xsl:variable>
         <xsl:value-of select="concat($legend,'!')"/>
      </xsl:for-each>
   </xsl:for-each>
</xsl:variable>

<xsl:call-template name="legends">
   <xsl:with-param name="lll" select="$lll"/>
   <xsl:with-param name="count" select="1"/>
   <xsl:with-param name="counter" select="$counter"/>
</xsl:call-template>

   <xsl:for-each select="exsl:node-set($content)/dos/*/diagram">
<xsl:text>
@target G0.S</xsl:text><xsl:value-of select="position()"/>
<xsl:text>
@type xy
</xsl:text>
      <xsl:for-each select="./point">
         <xsl:value-of select="@e"/>
<xsl:text>          </xsl:text>
         <xsl:value-of select="@dos"/>
<xsl:text>
</xsl:text>
      </xsl:for-each>
<xsl:text>&amp;
</xsl:text>
   </xsl:for-each>

</xsl:template>
<!-- ************************************************************************************************ -->
<xsl:template name="legends">
   <xsl:param name="lll"/>
   <xsl:param name="count"/>
   <xsl:param name="counter"/>
   <xsl:variable name="legend" select="substring-before($lll,'!')"/>
   <xsl:variable name="lll1" select="substring-after($lll,'!')"/>
<xsl:choose>
   <xsl:when test="$count&lt;$counter">
<xsl:text>
@s</xsl:text><xsl:value-of select="$count"/><xsl:text> legend "</xsl:text><xsl:value-of select="$legend"/><xsl:text>"</xsl:text>
   <xsl:call-template name="legends">
      <xsl:with-param name="lll" select="$lll1"/>
      <xsl:with-param name="count" select="$count+1"/>
      <xsl:with-param name="counter" select="$counter"/>
   </xsl:call-template>
   </xsl:when>
   <xsl:otherwise>
<xsl:text>
@s</xsl:text><xsl:value-of select="$count"/><xsl:text> legend "</xsl:text><xsl:value-of select="substring-before($lll,'!')"/><xsl:text>"</xsl:text>
   </xsl:otherwise>
</xsl:choose>
</xsl:template>
<!-- ************************************************************************************************ -->
</xsl:stylesheet>
