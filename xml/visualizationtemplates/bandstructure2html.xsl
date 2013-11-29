<?xml version="1.0" encoding="UTF-8" ?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
  <xsl:output method="html" />
  <xsl:template match="/">
<html> 
 <head> 
    <meta http-equiv="Content-Type" content="text/html; charset=utf-8" /> 
    <title>Bandstructure <xsl:value-of select="//title"/> </title> 
    <!--[if IE]><script language="javascript" type="text/javascript" src="../excanvas.min.js"></script><![endif]--> 
 
    <script language="javascript" type="text/javascript" src="http://static.exciting-code.org/flot/jquery.js"></script> 
    <script language="javascript" type="text/javascript" src="http://static.exciting-code.org/flot/jquery.flot.js"></script> 
 </head> 
    <body> 
    <h1>Bandstructure <xsl:value-of select="//title"/></h1> 
 
    <div id="placeholder" style="width:500px;height:600px"></div> 
 
    
<script id="source" language="javascript" type="text/javascript"> 

$(function () {

 <xsl:for-each select="//band">
 <xsl:variable name="bandnr" select="position()"></xsl:variable>
  var d<xsl:value-of select="$bandnr"/>= [];
 <xsl:for-each select="point">
   d<xsl:value-of select="$bandnr"/>.push([<xsl:value-of select="@distance"/>,<xsl:value-of select="@eval"/>]);
</xsl:for-each>  
</xsl:for-each> 
   
 
    
    $.plot($("#placeholder"), [
    <xsl:for-each select="//band">
        {  data: d<xsl:value-of select="position()"/>}<xsl:if test="not(position()=count(//band))">,</xsl:if>
        </xsl:for-each>
    ], {
        series: {
            lines: { show: true },
            points: { show: false }
        },
        xaxis: {
            ticks: [   <xsl:for-each select="//vertex">[<xsl:value-of select="@distance"/>, "<xsl:value-of select="@label"/>"]<xsl:if test="not(position()=count(//vertex))">,</xsl:if>
            </xsl:for-each>]
        },
        yaxis: {
            ticks:[[10,10],[0,"Fermi energy"]]
        },
        grid: { grid: { hoverable: true, clickable: true },
            backgroundColor: { colors: ["#fff", "#eee"] }
        }
    });
});


</script> 
 
 </body> 
</html> 
</xsl:template>
</xsl:stylesheet>
