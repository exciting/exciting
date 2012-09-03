<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
                xmlns:exsl="http://exslt.org/common"
                extension-element-prefixes="exsl"
                version="1.0">
<xsl:output method="xml" media-type="application/octet-stream" indent="yes"/>
<xsl:param name="ID"/>
<!--*************************************************************************************-->
<xsl:template match="/">
  <xsl:element name="dos">
    <xsl:text></xsl:text>
    <title><xsl:value-of select="/dos/title"/></title>
    <xsl:text></xsl:text>
    <xsl:choose>
      <xsl:when test="contains($ID,',')='true'">
	<xsl:call-template name="ID_parser">
	  <xsl:with-param name="id1"><xsl:value-of select="substring-before($ID,',')"/></xsl:with-param>
	  <xsl:with-param name="id2"><xsl:value-of select="substring-after($ID,',')"/></xsl:with-param>
	</xsl:call-template>
      </xsl:when>
      <xsl:otherwise>                   <!--"contains($ID,',')='false'"-->
	<xsl:call-template name="main">
	  <xsl:with-param name="id1"><xsl:value-of select="$ID"/></xsl:with-param>
	</xsl:call-template>
      </xsl:otherwise>
    </xsl:choose>
   </xsl:element>
</xsl:template>
<!--*************************************************************************************-->
<xsl:template name="ID_parser">
   <xsl:param name="id1"/>
   <xsl:param name="id2"/>
   <xsl:call-template name="main">
      <xsl:with-param name="id1"> <xsl:value-of select="$id1"/> </xsl:with-param>
   </xsl:call-template>
   <xsl:choose>
      <xsl:when test="contains($id2,',')='true'">
         <xsl:call-template name="ID_parser">
            <xsl:with-param name="id1"><xsl:value-of select="substring-before($id2,',')"/></xsl:with-param>
            <xsl:with-param name="id2"><xsl:value-of select="substring-after($id2,',')"/></xsl:with-param>
         </xsl:call-template>
      </xsl:when>
      <xsl:otherwise>   <!-- test="contains($id2,',')='false'">-->
         <xsl:call-template name="main">
            <xsl:with-param name="id1"> <xsl:value-of select="$id2"/> </xsl:with-param>
         </xsl:call-template>
      </xsl:otherwise>  
   </xsl:choose>
</xsl:template>
<!--*************************************************************************************-->
<xsl:template name="main">
   <xsl:param name="id1"/>
   <xsl:variable name="sym.atom" select="substring-before($id1,'/')"/>
   <xsl:variable name="sym"      select="substring-before($sym.atom,'.')"/>
   <xsl:variable name="n1_n2"    select="substring-after($sym.atom,'.')"/>
   <xsl:variable name="n1"       select="substring-before($n1_n2,':')"/>
   <xsl:variable name="n2"       select="substring-after($n1_n2,':')"/>
   <xsl:variable name="s1_s2"    select="substring-before(substring-after($id1,'/'),'/')"/>
   <xsl:variable name="s1"       select="substring-before($s1_s2,':')"/>
   <xsl:variable name="s2"       select="substring-after($s1_s2,':')"/>
   <xsl:variable name="l1_l2"    select="substring-before(substring-after(substring-after($id1,'/'),'/'),'/')"/>
   <xsl:variable name="l1"       select="substring-before($l1_l2,':')"/>
   <xsl:variable name="l2"       select="substring-after($l1_l2,':')"/>
   <xsl:variable name="m1_m2"    select="substring-after(substring-after(substring-after($id1,'/'),'/'),'/')"/>
   <xsl:variable name="m1"       select="substring-before($m1_m2,':')"/>
   <xsl:variable name="m2"       select="substring-after($m1_m2,':')"/>
   <xsl:choose>
      <xsl:when test="$sym.atom=''">
         <xsl:call-template name="total">
            <xsl:with-param name="s1_s2"> <xsl:value-of select="$s1_s2"/> </xsl:with-param>
            <xsl:with-param name="s1">    <xsl:value-of select="$s1"/>    </xsl:with-param>
            <xsl:with-param name="s2">    <xsl:value-of select="$s2"/>    </xsl:with-param>
         </xsl:call-template>
         <xsl:call-template name="interstitial">
            <xsl:with-param name="s1_s2"> <xsl:value-of select="$s1_s2"/> </xsl:with-param>
            <xsl:with-param name="s1">    <xsl:value-of select="$s1"/>    </xsl:with-param>
            <xsl:with-param name="s2">    <xsl:value-of select="$s2"/>    </xsl:with-param>
         </xsl:call-template>
         <xsl:call-template name="partial">
            <xsl:with-param name="sym.atom"> <xsl:value-of select="$sym.atom"/> </xsl:with-param>
            <xsl:with-param name="sym">      <xsl:value-of select="$sym"/>      </xsl:with-param>
            <xsl:with-param name="n1_n2">    <xsl:value-of select="$n1_n2"/>    </xsl:with-param>
            <xsl:with-param name="n1">       <xsl:value-of select="$n1"/>       </xsl:with-param>
            <xsl:with-param name="n2">       <xsl:value-of select="$n2"/>       </xsl:with-param>
            <xsl:with-param name="s1_s2">    <xsl:value-of select="$s1_s2"/>    </xsl:with-param>
            <xsl:with-param name="s1">       <xsl:value-of select="$s1"/>       </xsl:with-param>
            <xsl:with-param name="s2">       <xsl:value-of select="$s2"/>       </xsl:with-param>
            <xsl:with-param name="l1_l2">    <xsl:value-of select="$l1_l2"/>    </xsl:with-param>
            <xsl:with-param name="l1">       <xsl:value-of select="$l1"/>       </xsl:with-param>
            <xsl:with-param name="l2">       <xsl:value-of select="$l2"/>       </xsl:with-param>
            <xsl:with-param name="m1_m2">    <xsl:value-of select="$m1_m2"/>    </xsl:with-param>
            <xsl:with-param name="m1">       <xsl:value-of select="$m1"/>       </xsl:with-param>
            <xsl:with-param name="m2">       <xsl:value-of select="$m2"/>       </xsl:with-param>
         </xsl:call-template>
      </xsl:when>
      <xsl:when test="$sym.atom='t'">
         <xsl:call-template name="total">
            <xsl:with-param name="s1_s2"> <xsl:value-of select="$s1_s2"/> </xsl:with-param>
            <xsl:with-param name="s1">    <xsl:value-of select="$s1"/>    </xsl:with-param>
            <xsl:with-param name="s2">    <xsl:value-of select="$s2"/>    </xsl:with-param>
         </xsl:call-template>
      </xsl:when>
      <xsl:when test="$sym.atom='i'">
         <xsl:call-template name="interstitial">
            <xsl:with-param name="s1_s2"> <xsl:value-of select="$s1_s2"/> </xsl:with-param>
            <xsl:with-param name="s1">    <xsl:value-of select="$s1"/>    </xsl:with-param>
            <xsl:with-param name="s2">    <xsl:value-of select="$s2"/>    </xsl:with-param>
         </xsl:call-template>
      </xsl:when>
      <xsl:otherwise>
         <xsl:call-template name="partial">
            <xsl:with-param name="sym.atom"> <xsl:value-of select="$sym.atom"/> </xsl:with-param>
            <xsl:with-param name="sym">      <xsl:value-of select="$sym"/>      </xsl:with-param>
            <xsl:with-param name="n1_n2">    <xsl:value-of select="$n1_n2"/>    </xsl:with-param>
            <xsl:with-param name="n1">       <xsl:value-of select="$n1"/>       </xsl:with-param>
            <xsl:with-param name="n2">       <xsl:value-of select="$n2"/>       </xsl:with-param>
            <xsl:with-param name="s1_s2">    <xsl:value-of select="$s1_s2"/>    </xsl:with-param>
            <xsl:with-param name="s1">       <xsl:value-of select="$s1"/>       </xsl:with-param>
            <xsl:with-param name="s2">       <xsl:value-of select="$s2"/>       </xsl:with-param>
            <xsl:with-param name="l1_l2">    <xsl:value-of select="$l1_l2"/>    </xsl:with-param>
            <xsl:with-param name="l1">       <xsl:value-of select="$l1"/>       </xsl:with-param>
            <xsl:with-param name="l2">       <xsl:value-of select="$l2"/>       </xsl:with-param>
            <xsl:with-param name="m1_m2">    <xsl:value-of select="$m1_m2"/>    </xsl:with-param>
            <xsl:with-param name="m1">       <xsl:value-of select="$m1"/>       </xsl:with-param>
            <xsl:with-param name="m2">       <xsl:value-of select="$m2"/>       </xsl:with-param>
         </xsl:call-template>
      </xsl:otherwise>
   </xsl:choose>
</xsl:template>
<!--*************************************************************************************-->
<xsl:template name="total">
   <xsl:param name="s1_s2"/>
   <xsl:param name="s1"/>
   <xsl:param name="s2"/>
<xsl:element name="totaldos">
<xsl:text>
</xsl:text>
   <xsl:for-each select="/dos/totaldos/diagram">
      <xsl:choose>
         <xsl:when test="$s1='' or $s2=''">
            <xsl:if test="$s1_s2=@nspin or $s1_s2=''">
               <xsl:copy-of select="."/>
<xsl:text>
</xsl:text>
            </xsl:if>
         </xsl:when>
         <xsl:otherwise>
            <xsl:if test="@nspin&lt;=$s2 and @nspin&gt;=$s1">
               <xsl:copy-of select="."/>
<xsl:text>
</xsl:text>
            </xsl:if>
         </xsl:otherwise>
      </xsl:choose>
   </xsl:for-each>
</xsl:element>
<xsl:text>
</xsl:text>
</xsl:template>
<!--************************************************************************************-->
<xsl:template name="interstitial">
   <xsl:param name="s1_s2"/>
   <xsl:param name="s1"/>   
   <xsl:param name="s2"/>
<xsl:element name="interstitialdos">
<xsl:text>
</xsl:text>
   <xsl:for-each select="/dos/interstitialdos/diagram">
      <xsl:choose>
         <xsl:when test="$s1='' or $s2=''">
            <xsl:if test="$s1_s2=@nspin or $s1_s2=''">
               <xsl:copy-of select="."/>
<xsl:text>
</xsl:text>
            </xsl:if>
         </xsl:when>
         <xsl:otherwise>
            <xsl:if test="@nspin&lt;=$s2 and @nspin&gt;=$s1">
               <xsl:copy-of select="."/>
<xsl:text>
</xsl:text>
            </xsl:if>
         </xsl:otherwise>
       </xsl:choose>
   </xsl:for-each>
</xsl:element>
<xsl:text>
</xsl:text>
</xsl:template>
<!--************************************************************************************-->
<xsl:template name="partial">
   <xsl:param name="sym.atom"/>
   <xsl:param name="sym"/>
   <xsl:param name="n1_n2"/>
   <xsl:param name="n1"/>   
   <xsl:param name="n2"/>
   <xsl:param name="s1_s2"/>
   <xsl:param name="s1"/>   
   <xsl:param name="s2"/>
   <xsl:param name="l1_l2"/>
   <xsl:param name="l1"/>   
   <xsl:param name="l2"/>
   <xsl:param name="m1_m2"/>
   <xsl:param name="m1"/>   
   <xsl:param name="m2"/>
   <xsl:for-each select="/dos/partialdos">
   <xsl:choose>
      <xsl:when test="$sym=''">
         <xsl:choose>
            <xsl:when test="$sym.atom=''">
               <xsl:call-template name="condition_n">
                  <xsl:with-param name="n1_n2"> <xsl:value-of select="$n1_n2"/> </xsl:with-param>
                  <xsl:with-param name="n1">    <xsl:value-of select="$n1"/>    </xsl:with-param>
                  <xsl:with-param name="n2">    <xsl:value-of select="$n2"/>    </xsl:with-param>
                  <xsl:with-param name="s1_s2"> <xsl:value-of select="$s1_s2"/> </xsl:with-param>
                  <xsl:with-param name="s1">    <xsl:value-of select="$s1"/>    </xsl:with-param>
                  <xsl:with-param name="s2">    <xsl:value-of select="$s2"/>    </xsl:with-param>
                  <xsl:with-param name="l1_l2"> <xsl:value-of select="$l1_l2"/> </xsl:with-param>
                  <xsl:with-param name="l1">    <xsl:value-of select="$l1"/>    </xsl:with-param>
                  <xsl:with-param name="l2">    <xsl:value-of select="$l2"/>    </xsl:with-param>
                  <xsl:with-param name="m1_m2"> <xsl:value-of select="$m1_m2"/> </xsl:with-param>
                  <xsl:with-param name="m1">    <xsl:value-of select="$m1"/>    </xsl:with-param>
                  <xsl:with-param name="m2">    <xsl:value-of select="$m2"/>    </xsl:with-param>
               </xsl:call-template>
            </xsl:when>
            <xsl:otherwise>    <!--"$sym.atom!=''"-->
               <xsl:if test="$sym.atom=@speciessym">
                  <xsl:call-template name="condition_n">
                     <xsl:with-param name="n1_n2"> <xsl:value-of select="$n1_n2"/> </xsl:with-param>
                     <xsl:with-param name="n1">    <xsl:value-of select="$n1"/>    </xsl:with-param>
                     <xsl:with-param name="n2">    <xsl:value-of select="$n2"/>    </xsl:with-param>
                     <xsl:with-param name="s1_s2"> <xsl:value-of select="$s1_s2"/> </xsl:with-param>
                     <xsl:with-param name="s1">    <xsl:value-of select="$s1"/>    </xsl:with-param>
                     <xsl:with-param name="s2">    <xsl:value-of select="$s2"/>    </xsl:with-param>
                     <xsl:with-param name="l1_l2"> <xsl:value-of select="$l1_l2"/> </xsl:with-param>
                     <xsl:with-param name="l1">    <xsl:value-of select="$l1"/>    </xsl:with-param>
                     <xsl:with-param name="l2">    <xsl:value-of select="$l2"/>    </xsl:with-param>
                     <xsl:with-param name="m1_m2"> <xsl:value-of select="$m1_m2"/> </xsl:with-param>
                     <xsl:with-param name="m1">    <xsl:value-of select="$m1"/>    </xsl:with-param>
                     <xsl:with-param name="m2">    <xsl:value-of select="$m2"/>    </xsl:with-param>
                  </xsl:call-template>
               </xsl:if>
            </xsl:otherwise>
         </xsl:choose>
      </xsl:when>
      <xsl:otherwise>     <!-- $sym!='' -->
         <xsl:if test="$sym=@speciessym">
            <xsl:call-template name="condition_n">
               <xsl:with-param name="n1_n2"> <xsl:value-of select="$n1_n2"/> </xsl:with-param>
               <xsl:with-param name="n1">    <xsl:value-of select="$n1"/>    </xsl:with-param>
               <xsl:with-param name="n2">    <xsl:value-of select="$n2"/>    </xsl:with-param>
               <xsl:with-param name="s1_s2"> <xsl:value-of select="$s1_s2"/> </xsl:with-param>
               <xsl:with-param name="s1">    <xsl:value-of select="$s1"/>    </xsl:with-param>
               <xsl:with-param name="s2">    <xsl:value-of select="$s2"/>    </xsl:with-param>
               <xsl:with-param name="l1_l2"> <xsl:value-of select="$l1_l2"/> </xsl:with-param>
               <xsl:with-param name="l1">    <xsl:value-of select="$l1"/>    </xsl:with-param>
               <xsl:with-param name="l2">    <xsl:value-of select="$l2"/>    </xsl:with-param>
               <xsl:with-param name="m1_m2"> <xsl:value-of select="$m1_m2"/> </xsl:with-param>
               <xsl:with-param name="m1">    <xsl:value-of select="$m1"/>    </xsl:with-param>
               <xsl:with-param name="m2">    <xsl:value-of select="$m2"/>    </xsl:with-param>
            </xsl:call-template>
         </xsl:if>
      </xsl:otherwise>
   </xsl:choose>
   </xsl:for-each>
</xsl:template>
<!--************************************************************************************-->
<xsl:template name="condition_n">
   <xsl:param name="n1_n2"/>
   <xsl:param name="n1"/>
   <xsl:param name="n2"/>
   <xsl:param name="s1_s2"/>
   <xsl:param name="s1"/>
   <xsl:param name="s2"/>
   <xsl:param name="l1_l2"/>
   <xsl:param name="l1"/>
   <xsl:param name="l2"/>
   <xsl:param name="m1_m2"/>
   <xsl:param name="m1"/>
   <xsl:param name="m2"/>
   <xsl:choose>
      <xsl:when test="$n1!='' and $n2!=''">
         <xsl:if test="@atom&gt;=$n1 and @atom&lt;=$n2">
         <xsl:element name="partialdos">
            <xsl:attribute name="type">        <xsl:value-of select="@type"/>       </xsl:attribute>
            <xsl:attribute name="speciessym">  <xsl:value-of select="@speciessym"/> </xsl:attribute>
            <xsl:attribute name="speciesrn">   <xsl:value-of select="@speciesrn"/>  </xsl:attribute>
            <xsl:attribute name="atom">        <xsl:value-of select="@atom"/>       </xsl:attribute>
<xsl:text>
</xsl:text>
            <xsl:call-template name="condition_s">
               <xsl:with-param name="s1_s2"> <xsl:value-of select="$s1_s2"/> </xsl:with-param>
               <xsl:with-param name="s1">    <xsl:value-of select="$s1"/>    </xsl:with-param>
               <xsl:with-param name="s2">    <xsl:value-of select="$s2"/>    </xsl:with-param>
               <xsl:with-param name="l1_l2"> <xsl:value-of select="$l1_l2"/> </xsl:with-param>
               <xsl:with-param name="l1">    <xsl:value-of select="$l1"/>    </xsl:with-param>
               <xsl:with-param name="l2">    <xsl:value-of select="$l2"/>    </xsl:with-param>
               <xsl:with-param name="m1_m2"> <xsl:value-of select="$m1_m2"/> </xsl:with-param>
               <xsl:with-param name="m1">    <xsl:value-of select="$m1"/>    </xsl:with-param>
               <xsl:with-param name="m2">    <xsl:value-of select="$m2"/>    </xsl:with-param>
            </xsl:call-template>
         </xsl:element>
<xsl:text>
</xsl:text>
         </xsl:if>
      </xsl:when>
      <xsl:otherwise>                           <!--"$n1='' or $n2=''"-->
         <xsl:if test="$n1_n2=@atom or $n1_n2=''">
         <xsl:element name="partialdos">
            <xsl:attribute name="type">        <xsl:value-of select="@type"/>       </xsl:attribute>
            <xsl:attribute name="speciessym">  <xsl:value-of select="@speciessym"/> </xsl:attribute>
            <xsl:attribute name="speciesrn">   <xsl:value-of select="@speciesrn"/>  </xsl:attribute>
            <xsl:attribute name="atom">        <xsl:value-of select="@atom"/>       </xsl:attribute>
<xsl:text>
</xsl:text>
            <xsl:call-template name="condition_s">
               <xsl:with-param name="s1_s2"> <xsl:value-of select="$s1_s2"/> </xsl:with-param>
               <xsl:with-param name="s1">    <xsl:value-of select="$s1"/>    </xsl:with-param>
               <xsl:with-param name="s2">    <xsl:value-of select="$s2"/>    </xsl:with-param>
               <xsl:with-param name="l1_l2"> <xsl:value-of select="$l1_l2"/> </xsl:with-param>
               <xsl:with-param name="l1">    <xsl:value-of select="$l1"/>    </xsl:with-param>
               <xsl:with-param name="l2">    <xsl:value-of select="$l2"/>    </xsl:with-param>
               <xsl:with-param name="m1_m2"> <xsl:value-of select="$m1_m2"/> </xsl:with-param>
               <xsl:with-param name="m1">    <xsl:value-of select="$m1"/>    </xsl:with-param>
               <xsl:with-param name="m2">    <xsl:value-of select="$m2"/>    </xsl:with-param>
            </xsl:call-template>
         </xsl:element>
<xsl:text>
</xsl:text>
         </xsl:if>
      </xsl:otherwise>
   </xsl:choose>
</xsl:template>
<!--************************************************************************************-->
<xsl:template name="condition_s">
   <xsl:param name="s1_s2"/>
   <xsl:param name="s1"/>
   <xsl:param name="s2"/>
   <xsl:param name="l1_l2"/>
   <xsl:param name="l1"/>
   <xsl:param name="l2"/>
   <xsl:param name="m1_m2"/>
   <xsl:param name="m1"/>
   <xsl:param name="m2"/>
   <xsl:for-each select="./diagram">
      <xsl:choose>
         <xsl:when test="$s1!='' or $s2!=''">
            <xsl:if test="@nspin&gt;=$s1 and @nspin&lt;=$s2">
               <xsl:call-template name="condition_l">
                  <xsl:with-param name="l1_l2"> <xsl:value-of select="$l1_l2"/> </xsl:with-param>
                  <xsl:with-param name="l1">    <xsl:value-of select="$l1"/>    </xsl:with-param>
                  <xsl:with-param name="l2">    <xsl:value-of select="$l2"/>    </xsl:with-param>
                  <xsl:with-param name="m1_m2"> <xsl:value-of select="$m1_m2"/> </xsl:with-param>
                  <xsl:with-param name="m1">    <xsl:value-of select="$m1"/>    </xsl:with-param>
                  <xsl:with-param name="m2">    <xsl:value-of select="$m2"/>    </xsl:with-param>
               </xsl:call-template> 
            </xsl:if>
         </xsl:when>
         <xsl:otherwise>                               <!-- "$s1!='' and s2!=''" -->
            <xsl:if  test="$s1_s2='' or $s1_s2=@nspin">
               <xsl:call-template name="condition_l">
                  <xsl:with-param name="l1_l2"> <xsl:value-of select="$l1_l2"/> </xsl:with-param>
                  <xsl:with-param name="l1">    <xsl:value-of select="$l1"/>    </xsl:with-param>
                  <xsl:with-param name="l2">    <xsl:value-of select="$l2"/>    </xsl:with-param>
                  <xsl:with-param name="m1_m2"> <xsl:value-of select="$m1_m2"/> </xsl:with-param>
                  <xsl:with-param name="m1">    <xsl:value-of select="$m1"/>    </xsl:with-param>
                  <xsl:with-param name="m2">    <xsl:value-of select="$m2"/>    </xsl:with-param>
               </xsl:call-template>
            </xsl:if>
         </xsl:otherwise>
      </xsl:choose>
   </xsl:for-each>
</xsl:template>
<!--************************************************************************************-->
<xsl:template name="condition_l">
   <xsl:param name="l1_l2"/>
   <xsl:param name="l1"/>
   <xsl:param name="l2"/>
   <xsl:param name="m1_m2"/>
   <xsl:param name="m1"/>
   <xsl:param name="m2"/>
   <xsl:choose>
      <xsl:when test="$l1!='' or $l2!=''">
         <xsl:if test="@l&gt;=$l1 and @l&lt;=$l2">
            <xsl:call-template name="condition_m">
               <xsl:with-param name="m1_m2"> <xsl:value-of select="$m1_m2"/> </xsl:with-param>
               <xsl:with-param name="m1">    <xsl:value-of select="$m1"/>    </xsl:with-param>
               <xsl:with-param name="m2">    <xsl:value-of select="$m2"/>    </xsl:with-param>
            </xsl:call-template>
         </xsl:if>
      </xsl:when>
      <xsl:otherwise>                               <!-- "$l1!='' and l2!=''" -->
         <xsl:if  test="$l1_l2='' or $l1_l2=@l">
            <xsl:call-template name="condition_m">
               <xsl:with-param name="m1_m2"> <xsl:value-of select="$m1_m2"/> </xsl:with-param>
               <xsl:with-param name="m1">    <xsl:value-of select="$m1"/>    </xsl:with-param>
               <xsl:with-param name="m2">    <xsl:value-of select="$m2"/>    </xsl:with-param>
            </xsl:call-template>
         </xsl:if>
      </xsl:otherwise>
   </xsl:choose>
</xsl:template>
<!--************************************************************************************-->
<xsl:template name="condition_m">
   <xsl:param name="m1_m2"/>
   <xsl:param name="m1"/>
   <xsl:param name="m2"/>
   <xsl:choose>
      <xsl:when test="$m1!='' or $m2!=''">
         <xsl:if test="@m&gt;=$m1 and @m&lt;=$m2">
            <xsl:copy-of select="."/>
<xsl:text>
</xsl:text>
         </xsl:if>
      </xsl:when>
      <xsl:otherwise>                               <!-- "$m1!='' and m2!=''" -->
         <xsl:if test="$m1_m2='' or $m1_m2=@m">
            <xsl:copy-of select="."/>
<xsl:text>
</xsl:text>
         </xsl:if>
      </xsl:otherwise>
   </xsl:choose>
</xsl:template>
<!--************************************************************************************-->
</xsl:stylesheet>
