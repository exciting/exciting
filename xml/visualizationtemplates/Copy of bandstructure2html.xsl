<?xml version="1.0" encoding="UTF-8" ?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
  <xsl:output method="xml" />
  <xsl:template match="/">
<html>
  <head>
    <script type="text/javascript" src="http://www.google.com/jsapi"></script>
    <xsl:text>
    </xsl:text>
    <script type="text/javascript"> <xsl:text>
      google.load("visualization", "1", {packages:["linechart","scatterchart"]});
      google.setOnLoadCallback(drawChart);
      function drawChart() {
        var data = new google.visualization.DataTable();
        data.addColumn('number', 'distance');
       </xsl:text>
        <xsl:for-each select="//band">
        <xsl:text>data.addColumn('number', 'band </xsl:text><xsl:value-of select="position()"/> 
        <xsl:text>');
</xsl:text>
        </xsl:for-each>
        
         <xsl:text>
        data.addRows( </xsl:text><xsl:value-of select="count(//band[1]/point)"/> <xsl:text>); 
</xsl:text>
   
          <xsl:for-each select="//band">
               <xsl:variable name="bandnr" select="position()"/>
               <xsl:for-each select="point">
               <xsl:if test="$bandnr=1">
               <xsl:text>data.setValue(</xsl:text> <xsl:value-of select="position()-1"></xsl:value-of>, <xsl:value-of select="0"/>,<xsl:value-of select="@distance"/>);
</xsl:if>
            <xsl:text>data.setValue(</xsl:text>
             <xsl:value-of select="position()-1"/>,<xsl:value-of select="$bandnr"/>, <xsl:value-of select="@eval"/> 
            <xsl:text>);
</xsl:text>
          </xsl:for-each>
          </xsl:for-each>
        <xsl:text>
     
      var chart = new google.visualization.ScatterChart(document.getElementById('chart_div'));
        chart.draw(data, {width: 600, height: 800, titleX: 'distance', titleY: 'Energy', legend: 'none', pointSize: 0,lineSize:2}); }
</xsl:text>
</script>
  </head>

  <body>
    <div id="chart_div"></div>
  
  </body>
</html>
</xsl:template>
</xsl:stylesheet> 