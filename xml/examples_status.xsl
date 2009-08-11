<?xml version="1.0" encoding="UTF-8" ?>
<xsl:stylesheet version="1.0"
	xmlns:xsl="http://www.w3.org/1999/XSL/Transform">	
<xsl:output method="html" />
	<xsl:template match="/">
<html>
		<head>
			<script>			
						function loadXMLDoc(fname)
						{
							var xmlDoc;
							
  							xmlDoc=document.implementation.createDocument("","",null);
  							
							xmlDoc.async=false;
							xmlDoc.load(fname);
							
							return(xmlDoc);		
						}
						
						function displayResult(filepath_xml,filepath_xsl)
						{
						
							var fpath_xml = filepath_xml;
							var fpath_xsl = filepath_xsl;
							xml=loadXMLDoc(fpath_xml);
							
							xsl=loadXMLDoc(fpath_xsl);
							
							
							  xsltProcessor=new XSLTProcessor();
							  xsltProcessor.importStylesheet(xsl);
							  resultDocument = xsltProcessor.transformToFragment(xml,document);
							    
							  var ref = window.open("");
							  ref.document.write('<html>' + '<head>' + '</head>' + '<body>' + '<div id="test" />' + '</body>' + '</html>');
							  ref.document.getElementById("test").appendChild(resultDocument);
							  ref.document.close();				  
						}		
			</script>
		</head>
		<body>	
		<div id="test" />		
		<table border="2">
			<tr bgcolor="#0000CD">
				<th>n</th>
				<th>status</th>
				<th>name</th>
				<th>modification time</th>	
				<th>iterations</th>	
			</tr>
			<xsl:for-each select="//file[@name='info.xml']">
			<tr>	
				<xsl:variable name="PATH" select="/dirtree/head/path"></xsl:variable>
				<xsl:variable name="DIR"><xsl:value-of select="../@name"/></xsl:variable>
				<xsl:variable name="FILEPATH"><xsl:value-of select="/dirtree/head/path"/>/<xsl:value-of select="../@name"/>/info.xml</xsl:variable>
				<xsl:variable name="FILEPATH_in"><xsl:value-of select="/dirtree/head/path"/>/<xsl:value-of select="../@name"/>/input.xml</xsl:variable>
				<xsl:variable name="SHOW_INFO">displayResult("file://<xsl:value-of select="$FILEPATH"/>","file:///home/tom/git/exciting/examples/info.xsl")</xsl:variable>
				<xsl:variable name="SHOW_INPUT">displayResult("file://<xsl:value-of select="$FILEPATH_in"/>","file:///home/tom/git/exciting/examples/inputtohtml.xsl")</xsl:variable>
				<!-- <xsl:value-of select="$FILEPATH"></xsl:value-of> -->
				<!--  <xsl:value-of select="document($FILEPATH,info/groundstate)"></xsl:value-of> -->
				<th>
					<xsl:value-of select="position()"/>
				</th>
				<xsl:choose>
				<xsl:when test="document($FILEPATH,info/groundstate[@status='finished'])">
				<th bgcolor="#228B22">
					<Status>: calculation finished </Status>
				</th>
				</xsl:when>
				<xsl:otherwise>
				<th bgcolor = "#FF0000">
					<Status>: failed </Status>
				</th>
				</xsl:otherwise>
				</xsl:choose>
				<th>
					<Name> <xsl:value-of select="$DIR"></xsl:value-of> </Name>
				</th>
				<th>
					<Time> <xsl:value-of select="../modify-time"></xsl:value-of> </Time>
				</th>
				<th>
					<It_nr><xsl:value-of select="count(document($FILEPATH)/info/groundstate/scl/iter)"/></It_nr>
				</th>
				<th>
					<xsl:element name="a">
					<xsl:attribute name="href">#</xsl:attribute>
					<xsl:attribute name="onclick"><xsl:value-of select="$SHOW_INFO"/></xsl:attribute>
					<xsl:text>info</xsl:text>
					</xsl:element>
				</th>
				<th>
					<xsl:element name="a">
					<xsl:attribute name="href">#</xsl:attribute>
					<xsl:attribute name="onclick"><xsl:value-of select="$SHOW_INPUT"/></xsl:attribute>
					<xsl:text>input</xsl:text>
					</xsl:element>
				</th>
			</tr>
			</xsl:for-each>
		</table>
		<block>
			<th><xsl:value-of select="count(//file[@name='info.xml'])"/></th>
			<th> of </th>
			<th><xsl:value-of select="count(//file[@name='input.xml'])"/></th>
			<th> input files </th>
		</block>
		
		</body>	
</html>
	</xsl:template>
	
	
</xsl:stylesheet>