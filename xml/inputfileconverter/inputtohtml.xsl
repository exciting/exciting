<?xml version="1.0" encoding="UTF-8" ?>
<xsl:stylesheet version="1.0"
  xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
  <xsl:template match="/">
    <html>
      <head>
        <title>
          <xsl:value-of select="//title" />
        </title>
      <script type="text/javascript" src="http://ajax.googleapis.com/ajax/libs/jquery/1.2.6/jquery.min.js"/>
        <script src="http://static.exciting-code.org/Jmol.js" />
      </head>
      <body>
      
        <table border="0">

          <tr>
            <td>
              <h2>
                <xsl:value-of select="//title" />
              </h2>
              <script>
               <xsl:text>
        
         
          jmolInitialize("http://jmol.sourceforge.net/jmol");
          jmolCheckBrowser("popup", "http://jmol.sourceforge.net/browsercheck", "onClick");
                   
    
function trim(str, chars) {
	return ltrim(rtrim(str, chars), chars);
}
 
function ltrim(str, chars) {
	chars = chars || "\\s";
	return str.replace(new RegExp("^[" + chars + "]+", "g"), "");
}
 
function rtrim(str, chars) {
	chars = chars || "\\s";
	return str.replace(new RegExp("[" + chars + "]+$", "g"), "");
}

function transformlaticecoordxyz(coordstring){

var basevect=new Array(3);
var tmp;
var sep= new RegExp(" +");

</xsl:text>
                <xsl:for-each select="//basevect">
                   <xsl:text>tmp="</xsl:text>
                  <xsl:value-of select="." />
                   <xsl:text>";
</xsl:text>
                   <xsl:text>basevect[</xsl:text>
                  <xsl:value-of select="position()-1" />
                  <xsl:text>]=new Array(3);
</xsl:text>
                   <xsl:text>basevect[</xsl:text>
                  <xsl:value-of select="position()-1" />
                   <xsl:text>]=trim(tmp).split(sep);
</xsl:text>
                </xsl:for-each>
                <xsl:text>
var coord=trim(coordstring).split(sep);
var output = new Array(3);
output[0]=(coord[0] * basevect[0][0] + coord[1] * basevect[0][1]+ coord[2] * basevect[0][2])*0.529177249;
output[1]=(coord[0] * basevect[1][0] + coord[1] * basevect[1][1] + coord[2] * basevect[1][2])*0.529177249;
output[2]=(coord[0] * basevect[2][0] + coord[1] * basevect[2][1] + coord[2] * basevect[2][2])*0.529177249;
 return output.join(" ");
}
var structure="</xsl:text>
                <xsl:value-of select="count(//atom)" />
                <xsl:text>";
structure=structure+"|testing|";
</xsl:text>
                <xsl:for-each select="//species">

                  <xsl:variable name="symbol">
                    <xsl:choose>

                      <xsl:when test="@chemicalSymbol">
                        <xsl:value-of select="@chemicalSymbol" />
                      </xsl:when>
                      <xsl:otherwise>
                        <xsl:value-of select="'X'" />
                      </xsl:otherwise>
                    </xsl:choose>
                  </xsl:variable>
                  <xsl:for-each select="atom">
                    <xsl:text>structure=structure+"</xsl:text>
                    <xsl:value-of select="$symbol" />
                    <xsl:text> "+transformlaticecoordxyz("</xsl:text>
                    <xsl:value-of select="@coord" />
                    <xsl:text>")+"|";
</xsl:text>
                  </xsl:for-each>
                </xsl:for-each>

                <xsl:text>
    
  jmolAppletInline(400, structure, "background white;", "");


</xsl:text>

              </script>

            </td>
          </tr>
          <xsl:for-each
            select="/input/*[name()!='structure' and name()!='title']">
            <tr>
              <td>

                <xsl:call-template name="attributestotable" />

              </td>
            </tr>

          </xsl:for-each>
        </table>
      </body>
    </html>
  </xsl:template>
  <xsl:template name="attributestotable">
    <h2>
    <xsl:choose>
    <xsl:when test="name(..)='properties'">
    <xsl:element name="a">
    <xsl:attribute name="href">
    <xsl:text>./</xsl:text>      <xsl:value-of select="name(.)" />
    <xsl:text>.xml</xsl:text>
    </xsl:attribute>
      <xsl:value-of select="name(.)" />
    </xsl:element>
    </xsl:when>
    <xsl:otherwise>
          <xsl:value-of select="name(.)" />
    </xsl:otherwise>
    </xsl:choose>
    </h2>
    <table>
      <xsl:for-each select="@*">
        <tr>
          <td>
            <xsl:element name="a">
              <xsl:attribute name="href">
http://exciting-code.org/input-reference-essential-expert#att<xsl:value-of
                select="name(.)" />
	</xsl:attribute>
              <xsl:value-of select="name(.)" />
            </xsl:element>
          </td>
          <td>
            <xsl:value-of select="." />
          </td>
        </tr>
      </xsl:for-each>
      <xsl:for-each select="*">
        <tr>
          <td>
            <xsl:call-template name="attributestotable" />
          </td>
        </tr>
      </xsl:for-each>
    </table>
  </xsl:template>

</xsl:stylesheet>