<?xml version="1.0" encoding="UTF-8" ?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
  <xsl:template name="table">
    <xsl:param name="data" />
    <xsl:element name="table">
      <xsl:for-each select="$data">
        <xsl:element name="tr">
          <xsl:element name="td">
            <xsl:element name="a">
              <xsl:attribute name="href">
				    <xsl:variable name="attname" select="name(.)" />
				    <xsl:call-template name="agraph">
				<xsl:with-param name="data" select="//@*[name()=$attname]" />
				<xsl:with-param name="title" select="name(.)" />
				<xsl:with-param name="size" select="'600x300'" />
				</xsl:call-template>
				    </xsl:attribute>
				    <xsl:attribute name="target">
           				<xsl:text>_blank</xsl:text>
           			</xsl:attribute>
              <xsl:value-of select="name(.)" />
            </xsl:element>
          </xsl:element>
          <xsl:element name="td">
            <xsl:value-of select="." />
          </xsl:element>
        </xsl:element>
      </xsl:for-each>
    </xsl:element>
  </xsl:template>
  <xsl:template name="agraph">
    <xsl:param name="data" />
    <xsl:param name="title" />
    <xsl:param name="size" select="'300x120'" />
    <xsl:text>http://chart.apis.google.com/chart?cht=lc&amp;chs=</xsl:text> 
		<xsl:value-of select="$size"/>
		<xsl:text>&amp;chd=t:</xsl:text>
		<xsl:for-each select="$data">
			<xsl:value-of select="." />
			<xsl:if test="not(position()=count($data))">
				<xsl:text>,</xsl:text>
			</xsl:if>
		</xsl:for-each>
		<xsl:variable name="ymin">
			<xsl:for-each select="$data">
				<xsl:sort select="." order="ascending" data-type="number" />
				<xsl:if test="position()=1">
					<xsl:value-of select="." />
				</xsl:if>
			</xsl:for-each>
		</xsl:variable>
		<xsl:variable name="ymax">
			<xsl:for-each select="$data">
				<xsl:sort select="." order="descending" data-type="number" />
				<xsl:if test="position()=1">
					<xsl:value-of select="." />
				</xsl:if>
			</xsl:for-each>
		</xsl:variable>
		<xsl:text>&amp;chds=</xsl:text>
<xsl:value-of select="$ymin"/><xsl:text>,</xsl:text><xsl:value-of select="$ymax"/>
<xsl:text>&amp;chtt=</xsl:text><xsl:value-of select="$title"/>
<xsl:text>&amp;chxt=x,y&amp;chxr=0,1,</xsl:text>
<xsl:value-of select="count($data)"/><xsl:text>|1,</xsl:text>
<xsl:value-of select="$ymin"/><xsl:text>,</xsl:text><xsl:value-of select="$ymax"/>
</xsl:template>
	<xsl:template match="/">
		<html>
		<head>
		<xsl:element name="title">
		<xsl:value-of select="groundstate/scl/iter/charges/atom/@species"/> scl i=<xsl:value-of 
        select="count(groundstate/scl[last()]/iter)"/>
		<xsl:if test="groundstate/structure"> relax i=<xsl:value-of select="count(groundstate/structure)"/> </xsl:if>
		</xsl:element>
		
        </head>
			<body>
			    <xsl:element name="img">
			    <xsl:attribute name="src">
				<xsl:call-template name="agraph">
				<xsl:with-param name="data" select="/info/groundstate/scl[last()]/iter/energies/@totalEnergy" />
				<xsl:with-param name="title" select="'Total Energy'" />
				</xsl:call-template>
				</xsl:attribute>
				</xsl:element>
				
				<xsl:element name="img">
		    	<xsl:attribute name="src">
				<xsl:call-template name="agraph">
				<xsl:with-param name="data" select="/info/groundstate/scl[last()]/iter/@rmslog10" />
				<xsl:with-param name="title" select="'log10RMS+(convergence)'" />
				</xsl:call-template>
				</xsl:attribute>
				</xsl:element>
				
				<xsl:element name="img">
			    <xsl:attribute name="src">
				<xsl:call-template name="agraph">
				<xsl:with-param name="data" select="/info/groundstate/scl[last()]/iter/charges/@core_leakage" />
				<xsl:with-param name="title" select="'core leakage'" />
				</xsl:call-template>
				</xsl:attribute>
				</xsl:element>
				
				<xsl:element name="img">
			    <xsl:attribute name="src">
				<xsl:call-template name="agraph">
				<xsl:with-param name="data" select="/info/groundstate/scl[last()]/iter/energies/@fermiEnergy" />
				<xsl:with-param name="title" select="'Fermi energy'" />
				</xsl:call-template>
				</xsl:attribute>
				</xsl:element>			
				<xsl:if test="/info/groundstate/structure/@forceMax">
				<xsl:element name="img">
			    <xsl:attribute name="src">
				<xsl:call-template name="agraph">
				<xsl:with-param name="data" select="/info/groundstate/structure/@forceMax" />
				<xsl:with-param name="title" select="'maximal Force magnitude'" />
				</xsl:call-template>
				</xsl:attribute>
				</xsl:element>		
				</xsl:if>
				
				
				<h1>Results</h1>
        <table><tr><td valign = "top">
				<xsl:for-each select ="/info/groundstate/scl[last()]/iter[last()]">
				
				<xsl:element name="div">
				<xsl:attribute name="style"/>
				  <xsl:call-template name="table">
					<xsl:with-param name="data" select="@*"/>
					</xsl:call-template>				
					</xsl:element>
					<xsl:element name="div">
					<xsl:attribute name="style" />
					<h2>Energies</h2>
					<xsl:call-template name="table">
					<xsl:with-param name="data" select="energies/@*"/>
					</xsl:call-template>
					</xsl:element>
					<xsl:element name="div">
					 <xsl:attribute name="style"/> 
					<h2>Charges</h2>
					<xsl:call-template name="table">
					<xsl:with-param name="data" select="charges/@*"/>
					</xsl:call-template>
					</xsl:element>
					
					<xsl:element name="div">
					<xsl:attribute name="style"/>
					<h2>Timing</h2>
					<xsl:call-template name="table">
					<xsl:with-param name="data" select="timing/@*"/>
					</xsl:call-template>
					</xsl:element>
				
			
			    </xsl:for-each>
          </td><td valign = "top">
			    <xsl:if test="/info/groundstate/structure[last()]/species/atom/forces">
			    <xsl:element name="div">
					<xsl:attribute name="style"/>
				
					<xsl:for-each select="/info/groundstate/structure[last()]/species">
		
					<xsl:variable name="chemicalSymbol"><xsl:value-of select="@chemicalSymbol"/></xsl:variable> 
					<xsl:for-each select="atom">
					<xsl:variable name="anr"><xsl:value-of select="position()"/></xsl:variable>
					<h3> Species <xsl:value-of select="$chemicalSymbol"/> Atom <xsl:value-of select="position()"/></h3>
					<h4> Forces </h4>
					<table>
					<tr>
					<td>
          <xsl:element name="a">
          <xsl:attribute name="href">
        <xsl:call-template name="agraph">
        <xsl:with-param name="data" 
        select="/info/groundstate/structure/species[@chemicalSymbol=$chemicalSymbol]/atom[position()=$anr]/forces/@Magnitude" />
        <xsl:with-param name="title" ><xsl:text>species+</xsl:text>
        <xsl:value-of select="@chemicalSymbol"/><xsl:text>+atom+</xsl:text>
         <xsl:value-of select="$anr"/><xsl:text>+Force+Magnitude</xsl:text>
				</xsl:with-param>
				<xsl:with-param name="size" select="'600x300'" />
				</xsl:call-template>
           </xsl:attribute>
           <xsl:attribute name="target">
           	<xsl:text>_blank</xsl:text>
           </xsl:attribute>
					force magnitude  
          </xsl:element>
          </td>
          <td><xsl:value-of select="forces/@Magnitude"/></td>
					
					</tr>
					<xsl:for-each select="forces/*">
					<xsl:variable name="forcename"><xsl:value-of select="name(.)"/></xsl:variable>
					<tr>
					<td><xsl:value-of select="$forcename"/></td>
					<xsl:for-each select="@*">
					<xsl:variable name="attname"><xsl:value-of select="name(.)"/></xsl:variable>
					<td>
					<xsl:element name="a">
					<xsl:attribute name="href">
				 <xsl:call-template name="agraph">
				<xsl:with-param name="data" 
                select="/info/groundstate/structure/species[@chemicalSymbol=$chemicalSymbol]/atom[position()=$anr]/forces/*[name()=$forcename]/@*[name()=$attname]" />
				<xsl:with-param name="title" ><xsl:text>species+</xsl:text>
				<xsl:value-of select="$chemicalSymbol"/><xsl:text>+atom+</xsl:text>
<xsl:value-of select="$anr"/><xsl:text>+</xsl:text>
<xsl:value-of select="$forcename"/><xsl:text>+component+</xsl:text>
				<xsl:value-of select="name(.)" /> 
				</xsl:with-param>
				<xsl:with-param name="size" select="'600x300'" />
				</xsl:call-template>
					</xsl:attribute>
					<xsl:attribute name="target">
           				<xsl:text>_blank</xsl:text>
           			</xsl:attribute>
					<xsl:value-of select="name(.)"/>
					</xsl:element>
					<xsl:text>=</xsl:text>
					<xsl:value-of select="."/>
					</td>
					</xsl:for-each>
					</tr>
					</xsl:for-each>
					<tr>
					<td><b>Position</b></td>
					<xsl:for-each select="@*">
					<xsl:variable name="attname"><xsl:value-of select="name(.)"/></xsl:variable>
					<td>
					<xsl:element name="a">
					<xsl:attribute name="href">
					<xsl:call-template name="agraph">
				<xsl:with-param name="data" select="/info/groundstate/structure/species[@chemicalSymbol=$chemicalSymbol]/atom[position()=$anr]/@*[name()=$attname]" />
				<xsl:with-param name="title" ><xsl:text>species+</xsl:text>
				<xsl:value-of select="$chemicalSymbol"/><xsl:text>+atom+</xsl:text>
<xsl:value-of select="$anr"/><xsl:text>+</xsl:text>
<xsl:text>position+component+</xsl:text>
				<xsl:value-of select="name(.)" /> 
				</xsl:with-param>
				<xsl:with-param name="size" select="'600x300'" />
				</xsl:call-template>
					</xsl:attribute>
					<xsl:attribute name="target">
           				<xsl:text>_blank</xsl:text>
           			</xsl:attribute>
					<xsl:value-of select="name(.)"/>
					</xsl:element>
					<xsl:text>=</xsl:text>
					<xsl:value-of select="."/>
					</td>
					</xsl:for-each>
					</tr>
					</table>
					</xsl:for-each>
					</xsl:for-each>
					</xsl:element>
					</xsl:if>
					<xsl:if test="/info/groundstate/scl/iter/moments">
          <xsl:element name="div">
          <h2>Moments</h2>
          <table>
          <xsl:for-each select="/info/groundstate/scl/iter[last()]/moments/*[name()!='atom']">
           <xsl:variable name="momentname"><xsl:value-of select="name(.)"/> </xsl:variable>
          <tr>
          <td>
          <xsl:value-of select="name(.)"/>
          </td>
          <xsl:for-each select="@*">
          <xsl:variable name="attname"><xsl:value-of select="name(.)"/> </xsl:variable>
           <td>
           <xsl:element name="a">
           <xsl:attribute name="href">
           <xsl:call-template name="agraph">
        <xsl:with-param name="data" select="/info/groundstate/scl[last()]/iter/moments/*[name()=$momentname]/@*[name()=$attname]" />
        <xsl:with-param name="title" ><xsl:text>+</xsl:text>
 <xsl:value-of select="$momentname" /> 
<xsl:text>+component+</xsl:text>
        <xsl:value-of select="name(.)" /> 
        </xsl:with-param>
        <xsl:with-param name="size" select="'600x300'" />
        </xsl:call-template>
           </xsl:attribute>
           <xsl:attribute name="target">
           	<xsl:text>_blank</xsl:text>
           </xsl:attribute>
          <xsl:value-of select="name(.)"/>
          </xsl:element>
          <xsl:text>=</xsl:text>
           <xsl:value-of select="."/>
          </td>
          </xsl:for-each>
          </tr>
          </xsl:for-each>
          </table>
          <xsl:for-each select="/info/groundstate/scl/iter[last()]/moments/atom">
          <xsl:variable name="anr"> <xsl:value-of select="count(.)"/></xsl:variable>
           <h4> <xsl:value-of select="@species"/>  Atom   <xsl:value-of select="count(.)"/> </h4>
          <table>
          <xsl:for-each select="/info/groundstate/scl/iter[last()]/moments/atom/*">
          <xsl:variable name="momentname"><xsl:value-of select="name(.)"/> </xsl:variable>
          <tr>
          <td>
          <xsl:value-of select="name(.)"/>
          </td>
          <xsl:for-each select="@*">
           <xsl:variable name="attname"><xsl:value-of select="name(.)"/> </xsl:variable>
           <td>
           <xsl:element name="a">
           <xsl:attribute name="href">
           <xsl:call-template name="agraph">
        <xsl:with-param name="data" select="/info/groundstate/scl[last()]/iter/moments/atom[$anr]/*[name()=$momentname]/@*[name()=$attname]" />
        <xsl:with-param name="title" ><xsl:text>species+</xsl:text>
        <xsl:value-of select="@chemicalSymbol"/><xsl:text>+atom+</xsl:text>
<xsl:value-of select="$anr"/><xsl:text>+</xsl:text>
 <xsl:value-of select="$momentname" /> 
<xsl:text>+component+</xsl:text>
        <xsl:value-of select="name(.)" /> 
        </xsl:with-param>
        <xsl:with-param name="size" select="'600x300'" />
        </xsl:call-template>
           </xsl:attribute>
           <xsl:attribute name="target">
           	<xsl:text>_blank</xsl:text>
           </xsl:attribute>
          <xsl:value-of select="name(.)"/>
          </xsl:element>
          <xsl:text>=</xsl:text>
           <xsl:value-of select="."/>
          </td>
          </xsl:for-each>
          </tr>
          </xsl:for-each>
          </table>
          </xsl:for-each>
          </xsl:element>
                   </xsl:if>
                   </td></tr></table>
	            </body>
        	</html>	
	</xsl:template>
</xsl:stylesheet>